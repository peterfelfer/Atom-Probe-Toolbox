function metrics = dataQualityMetrics(pos, options)
% DATAQUALITYMETRICS Assess reconstruction quality and data integrity
%
% metrics = dataQualityMetrics(pos)
% metrics = dataQualityMetrics(pos, 'mass', mass)
%
% Computes various metrics to assess the quality of APT data and
% reconstruction, including spatial resolution estimates, density
% variations, and detection artifacts.
%
% INPUT:
%   pos - Nx3 array of atom positions [x, y, z] in nm
%         OR position table with x, y, z columns (and optionally mass, tof)
%
% OPTIONS:
%   'mass'           - Mass-to-charge values (for mass spectrum quality)
%   'tof'            - Time-of-flight values
%   'detectionEfficiency' - Known detection efficiency (default: estimate)
%   'voxelSize'      - Voxel size for density analysis in nm (default: 1)
%   'showPlots'      - Generate diagnostic plots (default: false)
%   'atomicVolume'   - Expected atomic volume in nm^3 (default: 0.012)
%
% OUTPUT:
%   metrics - Structure containing:
%       .summary - Overall quality assessment string
%       .resolution - Spatial resolution estimates
%           .lateral    - Estimated lateral resolution (nm)
%           .depth      - Estimated depth resolution (nm)
%           .method     - Method used for estimation
%       .density - Density analysis
%           .mean       - Mean atomic density (atoms/nm^3)
%           .std        - Standard deviation of local density
%           .cv         - Coefficient of variation
%           .expected   - Expected density based on atomic volume
%           .efficiency - Estimated detection efficiency
%       .geometry - Specimen geometry
%           .boundingBox - [minX minY minZ; maxX maxY maxZ]
%           .dimensions  - [lengthX lengthY lengthZ]
%           .volume      - Analyzed volume (nm^3)
%           .aspectRatio - Aspect ratios
%       .artifacts - Artifact detection
%           .densityHoles    - Suspected low-density regions
%           .densitySpikes   - Suspected high-density artifacts
%           .trajectoryOverlap - Detected trajectory overlap
%       .massSpectrum - Mass spectrum quality (if mass provided)
%           .peakToBackground - Peak to background ratio
%           .massResolution   - Estimated mass resolution
%           .thermalTails     - Thermal tail assessment
%
% EXAMPLES:
%   % Basic quality check
%   metrics = dataQualityMetrics(pos);
%   disp(metrics.summary);
%
%   % With mass spectrum analysis
%   metrics = dataQualityMetrics(posTable, 'mass', posTable.mc, 'showPlots', true);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos
    options.mass = []
    options.tof = []
    options.detectionEfficiency = []
    options.voxelSize (1,1) double {mustBePositive} = 1
    options.showPlots (1,1) logical = false
    options.atomicVolume (1,1) double {mustBePositive} = 0.012
end

% Handle table input
if istable(pos)
    posArray = [pos.x, pos.y, pos.z];
    if isempty(options.mass) && ismember('mc', pos.Properties.VariableNames)
        options.mass = pos.mc;
    end
else
    posArray = pos;
end

nAtoms = size(posArray, 1);

% Initialize metrics structure
metrics = struct();
metrics.nAtoms = nAtoms;

%% Geometry Analysis
fprintf('Analyzing geometry...\n');
metrics.geometry = analyzeGeometry(posArray);

%% Density Analysis
fprintf('Analyzing density distribution...\n');
metrics.density = analyzeDensity(posArray, options);

%% Resolution Estimation
fprintf('Estimating spatial resolution...\n');
metrics.resolution = estimateResolution(posArray, metrics.density);

%% Artifact Detection
fprintf('Checking for artifacts...\n');
metrics.artifacts = detectArtifacts(posArray, metrics.density, options);

%% Mass Spectrum Quality (if available)
if ~isempty(options.mass)
    fprintf('Analyzing mass spectrum quality...\n');
    metrics.massSpectrum = analyzeMassSpectrum(options.mass);
end

%% Generate Summary
metrics.summary = generateSummary(metrics);

%% Optional Plots
if options.showPlots
    generateDiagnosticPlots(posArray, metrics, options);
end

end

%% Geometry Analysis
function geom = analyzeGeometry(pos)
    geom = struct();

    % Bounding box
    minPos = min(pos, [], 1);
    maxPos = max(pos, [], 1);
    geom.boundingBox = [minPos; maxPos];
    geom.dimensions = maxPos - minPos;
    geom.center = (minPos + maxPos) / 2;

    % Volume (convex hull for better estimate)
    try
        [~, geom.convexHullVolume] = convhull(pos);
    catch
        geom.convexHullVolume = prod(geom.dimensions);
    end

    geom.boundingBoxVolume = prod(geom.dimensions);

    % Aspect ratios
    sortedDims = sort(geom.dimensions, 'descend');
    geom.aspectRatio = sortedDims(1) / sortedDims(3);
    geom.xyAspectRatio = geom.dimensions(1) / geom.dimensions(2);

    % Estimate tip radius (assuming conical specimen)
    % Typical APT specimen is conical, widening with depth (z)
    zLevels = linspace(minPos(3), maxPos(3), 20);
    radii = zeros(length(zLevels)-1, 1);
    for i = 1:length(zLevels)-1
        mask = pos(:,3) >= zLevels(i) & pos(:,3) < zLevels(i+1);
        if sum(mask) > 10
            xyPos = pos(mask, 1:2);
            radii(i) = sqrt(mean(sum((xyPos - mean(xyPos)).^2, 2)));
        end
    end
    geom.estimatedTipRadius = mean(radii(radii > 0));
    geom.radiusVariation = std(radii(radii > 0)) / mean(radii(radii > 0));
end

%% Density Analysis
function dens = analyzeDensity(pos, options)
    dens = struct();

    % Calculate local density using 3D histogram
    voxelSize = options.voxelSize;

    minPos = min(pos, [], 1);
    maxPos = max(pos, [], 1);

    % Create voxel grid
    nVoxels = ceil((maxPos - minPos) / voxelSize);
    nVoxels = max(nVoxels, 1);

    edges = cell(3, 1);
    for d = 1:3
        edges{d} = linspace(minPos(d), maxPos(d), nVoxels(d)+1);
    end

    % Count atoms per voxel
    [counts, ~] = histcounts2(pos(:,1), pos(:,2), edges{1}, edges{2});
    voxelVolume = voxelSize^3;

    % 3D density map
    density3D = zeros(nVoxels);
    for i = 1:nVoxels(1)
        for j = 1:nVoxels(2)
            zMask = pos(:,1) >= edges{1}(i) & pos(:,1) < edges{1}(i+1) & ...
                    pos(:,2) >= edges{2}(j) & pos(:,2) < edges{2}(j+1);
            for k = 1:nVoxels(3)
                mask = zMask & pos(:,3) >= edges{3}(k) & pos(:,3) < edges{3}(k+1);
                density3D(i,j,k) = sum(mask) / voxelVolume;
            end
        end
    end

    % Statistics (excluding empty edge voxels)
    validDensity = density3D(density3D > 0);
    dens.map3D = density3D;
    dens.voxelSize = voxelSize;

    dens.mean = mean(validDensity);
    dens.std = std(validDensity);
    dens.median = median(validDensity);
    dens.cv = dens.std / dens.mean;  % Coefficient of variation

    % Expected density from atomic volume
    dens.expected = 1 / options.atomicVolume;

    % Estimate detection efficiency
    if isempty(options.detectionEfficiency)
        dens.estimatedEfficiency = dens.mean / dens.expected;
    else
        dens.estimatedEfficiency = options.detectionEfficiency;
    end

    % Density profile along z (depth)
    zBins = linspace(minPos(3), maxPos(3), 50);
    dens.depthProfile = struct();
    dens.depthProfile.z = zBins(1:end-1) + diff(zBins)/2;
    dens.depthProfile.density = zeros(length(zBins)-1, 1);

    for i = 1:length(zBins)-1
        mask = pos(:,3) >= zBins(i) & pos(:,3) < zBins(i+1);
        sliceVolume = (maxPos(1)-minPos(1)) * (maxPos(2)-minPos(2)) * (zBins(i+1)-zBins(i));
        dens.depthProfile.density(i) = sum(mask) / sliceVolume;
    end
end

%% Resolution Estimation
function res = estimateResolution(pos, density)
    res = struct();

    % Method 1: From density and expected NN distance
    % In a perfect crystal, atoms are at regular spacing
    % Resolution ~ deviation from ideal positions

    % Nearest neighbor analysis
    nSample = min(size(pos, 1), 10000);
    sampleIdx = randperm(size(pos, 1), nSample);
    posSample = pos(sampleIdx, :);

    tree = KDTreeSearcher(posSample);
    [~, D] = knnsearch(tree, posSample, 'K', 2);
    nnDist = D(:, 2);

    % Expected NN distance for random distribution
    expectedNN = 0.554 * density.mean^(-1/3);

    % Lateral resolution estimate (from NN distribution width)
    res.lateral = std(nnDist);
    res.depth = res.lateral * 0.3;  % Typically better depth resolution

    % Method 2: From density fluctuations
    res.fromDensityCV = density.cv * expectedNN;

    res.method = 'nearest_neighbor_statistics';
    res.expectedNN = expectedNN;
    res.meanNN = mean(nnDist);
    res.stdNN = std(nnDist);

    % Quality indicator
    res.qualityIndex = res.meanNN / expectedNN;  % ~1 for good reconstruction
end

%% Artifact Detection
function art = detectArtifacts(pos, density, options)
    art = struct();

    % Low density regions (potential holes/voids)
    lowThreshold = density.mean * 0.3;
    art.densityHoles = struct();
    art.densityHoles.count = sum(density.map3D(:) < lowThreshold & density.map3D(:) > 0);
    art.densityHoles.fraction = art.densityHoles.count / numel(density.map3D);
    art.densityHoles.threshold = lowThreshold;

    % High density spikes (potential trajectory overlap)
    highThreshold = density.mean * 3;
    art.densitySpikes = struct();
    art.densitySpikes.count = sum(density.map3D(:) > highThreshold);
    art.densitySpikes.fraction = art.densitySpikes.count / numel(density.map3D);
    art.densitySpikes.threshold = highThreshold;

    % Trajectory overlap detection (multiple hits)
    % Check for atoms with very small separation
    nSample = min(size(pos, 1), 20000);
    sampleIdx = randperm(size(pos, 1), nSample);
    posSample = pos(sampleIdx, :);

    tree = KDTreeSearcher(posSample);
    [~, D] = knnsearch(tree, posSample, 'K', 2);
    nnDist = D(:, 2);

    overlapThreshold = 0.05;  % nm, very close atoms suggest overlap
    art.trajectoryOverlap = struct();
    art.trajectoryOverlap.count = sum(nnDist < overlapThreshold);
    art.trajectoryOverlap.fraction = art.trajectoryOverlap.count / nSample;
    art.trajectoryOverlap.threshold = overlapThreshold;

    % Overall artifact score (0-1, lower is better)
    art.overallScore = (art.densityHoles.fraction * 0.3 + ...
                        art.densitySpikes.fraction * 0.3 + ...
                        art.trajectoryOverlap.fraction * 0.4);
end

%% Mass Spectrum Quality
function ms = analyzeMassSpectrum(mass)
    ms = struct();

    % Basic histogram
    edges = 0:0.01:max(mass)*1.1;
    counts = histcounts(mass, edges);
    centers = edges(1:end-1) + 0.005;

    % Find peaks
    [pks, locs] = findpeaks(counts, 'MinPeakHeight', max(counts)*0.01, ...
                            'MinPeakDistance', 10);

    ms.nPeaks = length(pks);
    ms.peakPositions = centers(locs);
    ms.peakHeights = pks;

    % Background estimation (moving minimum)
    windowSize = 100;
    background = movmin(counts, windowSize);

    % Peak to background ratio
    if ~isempty(pks)
        bgAtPeaks = background(locs);
        ms.peakToBackground = mean(pks ./ max(bgAtPeaks, 1));
    else
        ms.peakToBackground = 0;
    end

    % Mass resolution estimate (FWHM of major peak)
    if ~isempty(pks)
        [~, maxPeakIdx] = max(pks);
        peakPos = locs(maxPeakIdx);

        % Find FWHM
        halfMax = pks(maxPeakIdx) / 2;
        leftIdx = find(counts(1:peakPos) < halfMax, 1, 'last');
        rightIdx = peakPos + find(counts(peakPos:end) < halfMax, 1, 'first') - 1;

        if ~isempty(leftIdx) && ~isempty(rightIdx)
            fwhm = centers(rightIdx) - centers(leftIdx);
            ms.massResolution = centers(peakPos) / fwhm;
            ms.fwhm = fwhm;
        else
            ms.massResolution = NaN;
            ms.fwhm = NaN;
        end
    else
        ms.massResolution = NaN;
        ms.fwhm = NaN;
    end

    % Thermal tail assessment (asymmetry of peaks)
    ms.tailAssessment = 'Not implemented';
end

%% Summary Generation
function summary = generateSummary(metrics)
    issues = {};
    quality = 'Good';

    % Check density
    if metrics.density.cv > 0.5
        issues{end+1} = 'High density variation (CV > 0.5)';
        quality = 'Fair';
    end

    % Check artifacts
    if metrics.artifacts.overallScore > 0.1
        issues{end+1} = sprintf('Artifacts detected (score: %.2f)', metrics.artifacts.overallScore);
        quality = 'Fair';
    end

    if metrics.artifacts.trajectoryOverlap.fraction > 0.01
        issues{end+1} = 'Significant trajectory overlap detected';
        quality = 'Poor';
    end

    % Check resolution
    if metrics.resolution.lateral > 0.5
        issues{end+1} = sprintf('Limited lateral resolution (~%.2f nm)', metrics.resolution.lateral);
    end

    % Check geometry
    if metrics.geometry.aspectRatio > 10
        issues{end+1} = 'Highly elongated specimen geometry';
    end

    % Build summary string
    summary = sprintf('=== Data Quality Assessment ===\n');
    summary = [summary, sprintf('Overall Quality: %s\n', quality)];
    summary = [summary, sprintf('Atoms: %d\n', metrics.nAtoms)];
    summary = [summary, sprintf('Volume: %.0f nm^3\n', metrics.geometry.convexHullVolume)];
    summary = [summary, sprintf('Mean density: %.1f atoms/nm^3\n', metrics.density.mean)];
    summary = [summary, sprintf('Est. detection efficiency: %.0f%%\n', metrics.density.estimatedEfficiency * 100)];
    summary = [summary, sprintf('Est. lateral resolution: %.2f nm\n', metrics.resolution.lateral)];

    if ~isempty(issues)
        summary = [summary, sprintf('\nIssues detected:\n')];
        for i = 1:length(issues)
            summary = [summary, sprintf('  - %s\n', issues{i})];
        end
    else
        summary = [summary, sprintf('\nNo significant issues detected.\n')];
    end
end

%% Diagnostic Plots
function generateDiagnosticPlots(pos, metrics, options)
    figure('Position', [100 100 1200 800], 'Name', 'Data Quality Diagnostics');

    % 1. Density histogram
    subplot(2,3,1);
    validDensity = metrics.density.map3D(metrics.density.map3D > 0);
    histogram(validDensity, 50, 'Normalization', 'pdf');
    hold on;
    xline(metrics.density.mean, 'r-', 'LineWidth', 2);
    xline(metrics.density.expected, 'g--', 'LineWidth', 2);
    xlabel('Local density (atoms/nm^3)');
    ylabel('Probability');
    title('Density Distribution');
    legend('Data', 'Mean', 'Expected', 'Location', 'best');

    % 2. Density along depth
    subplot(2,3,2);
    plot(metrics.density.depthProfile.z, metrics.density.depthProfile.density, 'b-', 'LineWidth', 1.5);
    hold on;
    yline(metrics.density.mean, 'r--');
    xlabel('Depth z (nm)');
    ylabel('Density (atoms/nm^3)');
    title('Density vs Depth');

    % 3. XY projection with density
    subplot(2,3,3);
    nSample = min(size(pos,1), 50000);
    sampleIdx = randperm(size(pos,1), nSample);
    scatter(pos(sampleIdx,1), pos(sampleIdx,2), 1, pos(sampleIdx,3), '.');
    axis equal;
    colorbar;
    xlabel('X (nm)');
    ylabel('Y (nm)');
    title('XY Projection (colored by Z)');

    % 4. NN distance distribution
    subplot(2,3,4);
    nSample = min(size(pos,1), 10000);
    sampleIdx = randperm(size(pos,1), nSample);
    tree = KDTreeSearcher(pos(sampleIdx,:));
    [~, D] = knnsearch(tree, pos(sampleIdx,:), 'K', 2);
    histogram(D(:,2), 50, 'Normalization', 'pdf');
    hold on;
    xline(metrics.resolution.expectedNN, 'r--', 'LineWidth', 2);
    xlabel('Nearest neighbor distance (nm)');
    ylabel('Probability');
    title('NN Distance Distribution');
    legend('Data', 'Expected (random)', 'Location', 'best');

    % 5. XZ slice
    subplot(2,3,5);
    yMid = (max(pos(:,2)) + min(pos(:,2))) / 2;
    yWidth = 5;  % nm
    mask = abs(pos(:,2) - yMid) < yWidth;
    scatter(pos(mask,1), pos(mask,3), 1, 'b.');
    axis equal;
    xlabel('X (nm)');
    ylabel('Z (nm)');
    title(sprintf('XZ Slice (Y = %.0f +/- %.0f nm)', yMid, yWidth));

    % 6. Summary text
    subplot(2,3,6);
    axis off;
    text(0.1, 0.9, 'Quality Metrics Summary', 'FontSize', 12, 'FontWeight', 'bold');
    text(0.1, 0.75, sprintf('Atoms: %d', metrics.nAtoms), 'FontSize', 10);
    text(0.1, 0.65, sprintf('Mean density: %.1f atoms/nm^3', metrics.density.mean), 'FontSize', 10);
    text(0.1, 0.55, sprintf('Density CV: %.2f', metrics.density.cv), 'FontSize', 10);
    text(0.1, 0.45, sprintf('Est. efficiency: %.0f%%', metrics.density.estimatedEfficiency*100), 'FontSize', 10);
    text(0.1, 0.35, sprintf('Lateral resolution: %.2f nm', metrics.resolution.lateral), 'FontSize', 10);
    text(0.1, 0.25, sprintf('Artifact score: %.2f', metrics.artifacts.overallScore), 'FontSize', 10);

    sgtitle('APT Data Quality Diagnostics');
end
