function [bgResult, info] = tofSpecBackgroundDetermination(pos, bgRegions, varargin)
% TOFSPECBACKGROUNDDETERMINATION Determine TOF background and track over experiment.
%
% Analyzes time-of-flight (TOF) background in specified regions and tracks
% how it changes over the course of the experiment. The background rate in
% TOF space can be converted to mass-to-charge space using the 1/sqrt(m/n)
% relationship.
%
% SYNTAX
%   [bgResult, info] = tofSpecBackgroundDetermination(pos, bgRegions)
%   [bgResult, info] = tofSpecBackgroundDetermination(pos, bgRegions, Name, Value, ...)
%
% INPUT
%   pos         - Table with atom probe data. Must contain either:
%                 - 'tof' column: time-of-flight values [ns], OR
%                 - 'mc' column: mass-to-charge values [Da] (requires 'tofConversion')
%                 Should contain 'ionIdx' for proper event ordering.
%
%   bgRegions   - Background regions specified as:
%                 - Nx2 array of [tofStart, tofEnd] pairs in ns, OR
%                 - 'auto': automatically detect valleys (requires 'mc' and 'ranges')
%                 - struct with fields 'tof' or 'mc' for explicit specification
%
% NAME-VALUE PAIRS
%   'blockSize'       - Number of ions per block (default: 1e6)
%   'numBlocks'       - Alternative: specify number of blocks (overrides blockSize)
%   'tofConversion'   - Conversion factor k where m/n = k * t^2 [Da/ns^2]
%                       If not provided and only 'mc' exists, will be estimated
%   'tofBinWidth'     - Bin width for TOF histogram [ns] (default: 1)
%   'mcBinWidth'      - Bin width for m/c histogram [Da] (default: 0.01)
%   'ranges'          - Range table (for 'auto' background detection)
%   'excludeRanged'   - Logical, exclude ranged regions from background (default: true)
%   'minPeakDistance' - Minimum distance from isolated peaks [Da] (default: 0.3)
%   'clusterThreshold' - Ranges closer than this are considered a cluster [Da] (default: 2.0)
%   'clusterExtraDistance' - Extra distance to add around peak clusters [Da] (default: 0.5)
%   'normalization'   - 'auto' | 'ion' | 'pulse' (default: 'auto')
%   'useUnrangedOnly' - Logical, estimate background from unranged ions only (default: true)
%   'tofRateMethod'   - 'regions' | 'histMean' (default: 'regions')
%   'plotResult'      - Logical, create diagnostic plots (default: false)
%   'mcRange'         - [mcMin, mcMax] range for output background curve (default: [0, 200])
%                       Note: mcRange is automatically extended to cover actual data range
%
% OUTPUT
%   bgResult    - Structure containing:
%                 .tofRate      - Background rate per block [counts/ns/normalization]
%                 .tofRateError - Standard error of tofRate per block
%                 .blockIdx     - Block indices [start, end] into pos table
%                 .blockCenter  - Center ion index of each block
%                 .mcCurve      - Table with m/c values and background per block
%                 .meanRate     - Mean background rate across all blocks
%                 .rateStd      - Standard deviation of rate across blocks
%                 .rateTrend    - Linear trend coefficient [rate change per event index]
%                 .normalization - 'ion' or 'pulse'
%                 .totalNormalization - Total normalization count
%
%   info        - Structure containing:
%                 .tofConversion   - Conversion factor k used
%                 .bgRegionsTof    - Background regions in TOF space
%                 .totalBgCounts   - Total counts in background regions per block
%                 .totalBgTime     - Total TOF range width used [ns]
%                 .numBlocks       - Number of blocks analyzed
%                 .method          - 'direct' or 'fromMC' depending on input
%                 .normalization   - 'ion' or 'pulse'
%                 .totalNormalization - Total normalization count
%                 .useUnrangedOnly - true if unranged-only estimation is used
%
% PHYSICS
%   The mass-to-charge ratio relates to time-of-flight as:
%       m/n = k * t^2
%
%   where k is an instrument constant [Da/ns^2].
%
%   If background is constant in TOF space at rate B_tof [counts/ns/normalization],
%   the background in mass spectrum space is:
%       B_mc(m/n) = B_tof / (2 * sqrt(k * m/n))
%
%   This gives the characteristic 1/sqrt(m/n) decay in mass spectra.
%
% EXAMPLE
%   % Define background regions in TOF space
%   bgRegions = [100, 120; 250, 280; 400, 450];  % ns
%
%   % Analyze background in 10 blocks
%   [bgResult, info] = tofSpecBackgroundDetermination(pos, bgRegions, ...
%       'numBlocks', 10, 'plotResult', true);
%
%   % Get background at m/c = 56 Da for block 5
%   mc56 = 56;
%   bg56 = bgResult.tofRate(5) / (2 * sqrt(info.tofConversion * mc56));
%
%   % Auto-detect background regions from ranges
%   [bgResult, info] = tofSpecBackgroundDetermination(pos, 'auto', ...
%       'ranges', ranges, 'numBlocks', 10);
%
% SEE ALSO
%   backgroundEstimate, massSpecPlot, posCalculateConcentrationDeconvolved
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% Parse inputs
p = inputParser;
addRequired(p, 'pos', @istable);
addRequired(p, 'bgRegions');
addParameter(p, 'blockSize', 1e6, @(x) isnumeric(x) && x > 0);
addParameter(p, 'numBlocks', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'tofConversion', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'tofBinWidth', 1, @(x) isnumeric(x) && x > 0);
addParameter(p, 'mcBinWidth', 0.01, @(x) isnumeric(x) && x > 0);
addParameter(p, 'ranges', [], @(x) isempty(x) || istable(x));
addParameter(p, 'excludeRanged', true, @islogical);
addParameter(p, 'minPeakDistance', 0.3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'clusterThreshold', 2.0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'clusterExtraDistance', 0.5, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'normalization', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'useUnrangedOnly', true, @islogical);
addParameter(p, 'tofRateMethod', 'regions', @(x) ischar(x) || isstring(x));
addParameter(p, 'plotResult', false, @islogical);
addParameter(p, 'mcRange', [0, 200], @(x) isnumeric(x) && numel(x) == 2);

parse(p, pos, bgRegions, varargin{:});
opts = p.Results;

%% Validate input data
hasTof = ismember('tof', pos.Properties.VariableNames);
hasMc = ismember('mc', pos.Properties.VariableNames);
hasIonIdx = ismember('ionIdx', pos.Properties.VariableNames);
hasDeltaP = ismember('deltaP', pos.Properties.VariableNames);

if ~hasTof && ~hasMc
    error('tofSpecBackgroundDetermination:noData', ...
        'pos table must contain either ''tof'' or ''mc'' column.');
end

% Determine method
if hasTof
    info.method = 'direct';
    tofData = pos.tof;
else
    info.method = 'fromMC';
    if isempty(opts.tofConversion)
        % Estimate k from the data - use median m/c and assume typical TOF
        % This is approximate; user should provide tofConversion for accuracy
        warning('tofSpecBackgroundDetermination:estimatingK', ...
            'Estimating TOF conversion factor. For accurate results, provide ''tofConversion'' parameter.');
        % Typical APT: k ≈ 0.0005 Da/ns^2 for ~100mm flight path
        opts.tofConversion = 5e-4;
    end
    % Convert m/c to TOF: t = sqrt(m/n / k)
    tofData = sqrt(pos.mc / opts.tofConversion);
end

info.tofConversion = opts.tofConversion;

% Get event ordering
if hasIonIdx
    [~, sortOrder] = sort(pos.ionIdx);
else
    warning('tofSpecBackgroundDetermination:noIonIdx', ...
        'No ionIdx column found. Using row order as event sequence.');
    sortOrder = (1:height(pos))';
end

numIons = height(pos);

% Determine normalization basis
normMode = lower(string(opts.normalization));
if hasDeltaP
    if normMode ~= "pulse"
        warning('tofSpecBackgroundDetermination:forcePulseNormalization', ...
            'deltaP found. Enforcing pulse normalization.');
    end
    normMode = "pulse";
else
    if normMode == "" || normMode == "auto"
        normMode = "ion";
    end
    if normMode == "pulse"
        warning('tofSpecBackgroundDetermination:noDeltaP', ...
            'Normalization by pulse requested but deltaP column not found. Using ion normalization.');
        normMode = "ion";
    end
end

if normMode == "pulse"
    totalNorm = sum(pos.deltaP, 'omitnan');
else
    totalNorm = numIons;
end
if ~isfinite(totalNorm) || totalNorm <= 0
    totalNorm = numIons;
    normMode = "ion";
end
info.normalization = char(normMode);
info.totalNormalization = totalNorm;
info.useUnrangedOnly = opts.useUnrangedOnly;
info.tofRateMethod = char(lower(string(opts.tofRateMethod)));

%% Process background regions
useUnrangedOnly = opts.useUnrangedOnly;
rateMethod = lower(string(opts.tofRateMethod));
bgRegionsTof = [];
bgRegionsTofPlot = [];
if useUnrangedOnly
    if isempty(opts.ranges)
        error('tofSpecBackgroundDetermination:noRanges', ...
            'useUnrangedOnly requires a range table (''ranges'').');
    end
end

if ischar(bgRegions) && strcmpi(bgRegions, 'auto')
    % Auto-detect background regions
    if isempty(opts.ranges)
        error('tofSpecBackgroundDetermination:noRanges', ...
            'Automatic background detection requires ''ranges'' parameter.');
    end
    [bgRegionsAuto, opts.mcRange] = autoDetectBgRegions(pos.mc, opts.ranges, opts.minPeakDistance, ...
        opts.tofConversion, opts.mcRange, opts.clusterThreshold, opts.clusterExtraDistance);
    bgRegionsTofPlot = bgRegionsAuto;
    if ~useUnrangedOnly && rateMethod == "regions"
        bgRegionsTof = bgRegionsAuto;
    end
elseif isstruct(bgRegions)
    if isfield(bgRegions, 'tof')
        bgRegionsTof = bgRegions.tof;
        bgRegionsTofPlot = bgRegionsTof;
    elseif isfield(bgRegions, 'mc')
        % Convert m/c regions to TOF
        bgRegionsTof = sqrt(bgRegions.mc / opts.tofConversion);
        bgRegionsTofPlot = bgRegionsTof;
    else
        error('tofSpecBackgroundDetermination:invalidRegions', ...
            'bgRegions struct must have ''tof'' or ''mc'' field.');
    end
else
    % Assume direct TOF regions
    bgRegionsTof = bgRegions;
    bgRegionsTofPlot = bgRegionsTof;
end

info.bgRegionsTof = bgRegionsTof;
info.bgRegionsTofPlot = bgRegionsTofPlot;
if isempty(bgRegionsTofPlot)
    bgRegionsTofPlot = bgRegionsTof;
    info.bgRegionsTofPlot = bgRegionsTofPlot;
end

% Calculate total TOF width of background regions
if isempty(bgRegionsTof)
    totalBgTofWidth = NaN;
    info.totalBgTime = NaN;
else
    totalBgTofWidth = sum(bgRegionsTof(:,2) - bgRegionsTof(:,1));
    info.totalBgTime = totalBgTofWidth;
end

%% Determine block boundaries
if ~isempty(opts.numBlocks)
    numBlocks = opts.numBlocks;
    blockSize = ceil(numIons / numBlocks);
else
    blockSize = opts.blockSize;
    numBlocks = ceil(numIons / blockSize);
end

info.numBlocks = numBlocks;

% Create block boundaries
blockIdx = zeros(numBlocks, 2);
for b = 1:numBlocks
    blockIdx(b, 1) = (b-1) * blockSize + 1;
    blockIdx(b, 2) = min(b * blockSize, numIons);
end

%% Analyze background per block
tofRate = zeros(numBlocks, 1);
tofRateError = zeros(numBlocks, 1);
totalBgCounts = zeros(numBlocks, 1);
bgWidthPerBlock = zeros(numBlocks, 1);

if useUnrangedOnly
    % Build unranged mask from range table using m/c values
    inAnyRange = false(height(pos), 1);
    for r = 1:height(opts.ranges)
        inAnyRange = inAnyRange | (pos.mc >= opts.ranges.mcbegin(r) & pos.mc <= opts.ranges.mcend(r));
    end
    isUnranged = ~inAnyRange;
end

for b = 1:numBlocks
    % Get ions in this block (in event order)
    blockIons = sortOrder(blockIdx(b,1):blockIdx(b,2));
    blockTof = tofData(blockIons);
    numBlockIons = length(blockIons);
    if normMode == "pulse"
        numBlockNorm = sum(pos.deltaP(blockIons), 'omitnan');
    else
        numBlockNorm = numBlockIons;
    end
    if ~isfinite(numBlockNorm) || numBlockNorm <= 0
        numBlockNorm = numBlockIons;
    end

    % Count ions in background regions or estimate from histogram mean
    bgCounts = 0;
    if rateMethod == "histmean"
        binCounts = [];
        if numel(blockTof) > 1
            tofEdges = min(blockTof):opts.tofBinWidth:max(blockTof);
            if numel(tofEdges) < 2
                tofEdges = [min(blockTof), min(blockTof) + opts.tofBinWidth];
            end
            [binCounts, ~] = histcounts(blockTof, tofEdges);
            binRate = binCounts / (opts.tofBinWidth * numBlockNorm);
            tofRate(b) = mean(binRate);
            tofRateError(b) = std(binRate) / max(1, sqrt(numel(binRate)));
            bgWidthPerBlock(b) = (numel(binCounts)) * opts.tofBinWidth;
        else
            tofRate(b) = 0;
            tofRateError(b) = 0;
            bgWidthPerBlock(b) = 0;
        end
        totalBgCounts(b) = sum(binCounts);
        continue;
    end

    if useUnrangedOnly
        blockUnranged = isUnranged(blockIons);
        bgCounts = sum(blockUnranged);
        if numel(blockTof) > 1
            bgWidthPerBlock(b) = max(blockTof) - min(blockTof);
        else
            bgWidthPerBlock(b) = 0;
        end
    else
        for r = 1:size(bgRegionsTof, 1)
            inRegion = blockTof >= bgRegionsTof(r,1) & blockTof < bgRegionsTof(r,2);
            bgCounts = bgCounts + sum(inRegion);
        end
        bgWidthPerBlock(b) = totalBgTofWidth;
    end

    totalBgCounts(b) = bgCounts;

    % Calculate background rate [counts/ns/normalization]
    % Rate = (counts in bg region) / (tof width) / (normalization in block)
    bgWidth = bgWidthPerBlock(b);
    if ~isfinite(bgWidth) || bgWidth <= 0
        bgWidth = eps;
    end
    tofRate(b) = bgCounts / bgWidth / numBlockNorm;

    % Standard error assuming Poisson statistics
    if bgCounts > 0
        tofRateError(b) = sqrt(bgCounts) / bgWidth / numBlockNorm;
    else
        tofRateError(b) = 1 / bgWidth / numBlockNorm; % Upper limit
    end
end

info.totalBgCounts = totalBgCounts;
info.totalBgTimePerBlock = bgWidthPerBlock;

%% Calculate statistics
bgResult.tofRate = tofRate;
bgResult.tofRateError = tofRateError;
bgResult.blockIdx = blockIdx;
bgResult.blockCenter = mean(blockIdx, 2);
bgResult.meanRate = mean(tofRate);
bgResult.rateStd = std(tofRate);
bgResult.normalization = char(normMode);
bgResult.totalNormalization = totalNorm;

% Linear trend
blockCenters = bgResult.blockCenter;
if numBlocks > 1
    p = polyfit(blockCenters, tofRate, 1);
    bgResult.rateTrend = p(1);  % Rate change per event index
else
    bgResult.rateTrend = 0;
end

%% Generate m/c background curve for each block
mcValues = (opts.mcRange(1) + opts.mcBinWidth/2 : opts.mcBinWidth : opts.mcRange(2))';
mcBg = zeros(length(mcValues), numBlocks);

k = opts.tofConversion;
if isempty(k)
    k = 5e-4;  % Default estimate
end

for b = 1:numBlocks
    % B_mc = B_tof / (2 * sqrt(k * m/n))
    mcBg(:, b) = tofRate(b) ./ (2 * sqrt(k * mcValues));
end

% Create table
mcCurveTable = array2table(mcValues, 'VariableNames', {'mc'});
for b = 1:numBlocks
    mcCurveTable.(['block' num2str(b)]) = mcBg(:, b);
end
mcCurveTable.mean = mean(mcBg, 2);

bgResult.mcCurve = mcCurveTable;

%% Validate background against actual spectrum data
% Compare model prediction to observed counts in background regions
if ismember('mc', pos.Properties.VariableNames) && ~isempty(bgRegionsTofPlot)
    [validationStats, correctionFactor] = validateBackgroundFit(pos.mc, bgRegionsTofPlot, ...
        bgResult.meanRate, k, opts.mcBinWidth, opts.mcRange, totalNorm);
    bgResult.validationStats = validationStats;
    bgResult.correctionFactor = correctionFactor;
    bgResult.suggestedRate = bgResult.meanRate * correctionFactor;
else
    bgResult.validationStats = struct();
    bgResult.correctionFactor = 1.0;
    bgResult.suggestedRate = bgResult.meanRate;
end

%% Plotting
if opts.plotResult
    plotBackgroundAnalysis(bgResult, info, pos, tofData, opts);
end

end


function [stats, correctionFactor] = validateBackgroundFit(mc, bgRegionsTof, tofRate, k, binWidth, mcRange, totalNorm)
% Validate background fit by comparing model to observed counts in background regions
%
% Returns statistics per region and an overall correction factor

    % Convert TOF regions to m/c
    bgRegionsMc = bgRegionsTof.^2 * k;

    % Build histogram
    mcEdges = mcRange(1):binWidth:mcRange(2);
    [counts, ~] = histcounts(mc, mcEdges);
    mcCenters = mcEdges(1:end-1) + binWidth/2;
    if nargin < 7 || isempty(totalNorm)
        totalNorm = numel(mc);
    end
    countsNorm = counts / totalNorm;  % counts per bin per normalization

    % Background model: B_mc = tofRate / (2 * sqrt(k * m/c)) * binWidth
    bgModel = tofRate ./ (2 * sqrt(k * mcCenters)) * binWidth;

    % Analyze each background region
    stats = struct('mcStart', {}, 'mcEnd', {}, 'meanObs', {}, 'meanModel', {}, ...
        'ratio', {}, 'numBins', {});
    allRatios = [];

    for r = 1:size(bgRegionsMc, 1)
        mcStart = max(bgRegionsMc(r, 1), mcRange(1));
        mcEnd = min(bgRegionsMc(r, 2), mcRange(2));
        if mcEnd <= mcStart + binWidth
            continue;
        end

        inRegion = mcCenters >= mcStart & mcCenters <= mcEnd;
        if sum(inRegion) < 2
            continue;
        end

        regionObs = countsNorm(inRegion);
        regionModel = bgModel(inRegion);

        meanObs = mean(regionObs);
        meanModel = mean(regionModel);

        if meanModel > 0
            ratio = meanObs / meanModel;
        else
            ratio = NaN;
        end

        stats(end+1) = struct(...
            'mcStart', mcStart, ...
            'mcEnd', mcEnd, ...
            'meanObs', meanObs, ...
            'meanModel', meanModel, ...
            'ratio', ratio, ...
            'numBins', sum(inRegion));

        if ~isnan(ratio)
            allRatios = [allRatios; ratio];
        end
    end

    % Overall correction factor (use median to be robust to outliers)
    if ~isempty(allRatios)
        correctionFactor = median(allRatios);
    else
        correctionFactor = 1.0;
    end
end


%% Helper functions

function [bgRegionsTof, mcRange] = autoDetectBgRegions(mc, ranges, minDist, k, mcRange, clusterThreshold, clusterExtraDistance)
% Auto-detect background regions between peaks
% Handles touching/overlapping ranges by merging them first
% Keeps extra distance from peak clusters (groups of closely-spaced peaks)

if isempty(k)
    k = 5e-4;
end

if nargin < 6 || isempty(clusterThreshold)
    clusterThreshold = 2.0;  % Da - ranges closer than this are considered a cluster
end
if nargin < 7 || isempty(clusterExtraDistance)
    clusterExtraDistance = 0.5;  % Extra distance to add for clustered peaks
end

% Extend mcRange to cover actual data if needed
actualMax = max(mc);
actualMin = min(mc);
if mcRange(2) < actualMax
    mcRange(2) = actualMax + 5;  % Add 5 Da beyond last data point
end
if mcRange(1) > actualMin
    mcRange(1) = max(0, actualMin - 1);
end

% Sort ranges by start position (keeping begin/end paired)
rangeMatrix = [ranges.mcbegin, ranges.mcend];
rangeMatrix = sortrows(rangeMatrix, 1);
numRanges = size(rangeMatrix, 1);

leftDist = inf(numRanges, 1);
rightDist = inf(numRanges, 1);

for i = 1:numRanges
    if i > 1
        leftDist(i) = rangeMatrix(i, 1) - rangeMatrix(i-1, 2);
    end
    if i < numRanges
        rightDist(i) = rangeMatrix(i+1, 1) - rangeMatrix(i, 2);
    end
end

% Determine exclusion distance for each range edge
% Use larger distance if the range is part of a cluster
leftExclusion = zeros(numRanges, 1);  % no padding before peaks
rightExclusion = ones(numRanges, 1) * minDist;  % padding after peaks

for i = 1:numRanges
    % Check if this range is in a cluster (close neighbors on either side)
    inClusterLeft = leftDist(i) < clusterThreshold;
    inClusterRight = rightDist(i) < clusterThreshold;

    if inClusterLeft || inClusterRight
        rightExclusion(i) = minDist + clusterExtraDistance;
    end
end

% Merge overlapping/touching ranges (considering exclusion zones)
mergedRanges = [];
mergedExclusion = [];  % [leftExcl, rightExcl] for each merged range
currentStart = rangeMatrix(1, 1);
currentEnd = rangeMatrix(1, 2);
currentLeftExcl = leftExclusion(1);
currentRightExcl = rightExclusion(1);

for i = 2:numRanges
    nextStart = rangeMatrix(i, 1);
    nextEnd = rangeMatrix(i, 2);

    % Check if ranges touch or overlap (with tolerance)
    if nextStart <= currentEnd + 0.01
        % Merge: extend current range, keep max exclusion
        currentEnd = max(currentEnd, nextEnd);
        currentRightExcl = rightExclusion(i);
    else
        % Save current range and start new one
        mergedRanges = [mergedRanges; currentStart, currentEnd];
        mergedExclusion = [mergedExclusion; currentLeftExcl, currentRightExcl];
        currentStart = nextStart;
        currentEnd = nextEnd;
        currentLeftExcl = leftExclusion(i);
        currentRightExcl = rightExclusion(i);
    end
end
% Don't forget the last range
mergedRanges = [mergedRanges; currentStart, currentEnd];
mergedExclusion = [mergedExclusion; currentLeftExcl, currentRightExcl];

% Now find gaps between merged ranges
bgRegionsMc = [];
numMerged = size(mergedRanges, 1);

% Before first range
firstExcl = mergedExclusion(1, 1);
if mergedRanges(1, 1) - firstExcl > mcRange(1)
    regionStart = max(mcRange(1), 0.1);  % Avoid zero for 1/sqrt model
    regionEnd = mergedRanges(1, 1) - firstExcl;
    if regionEnd > regionStart + 0.1  % Allow smaller edge region
        bgRegionsMc = [bgRegionsMc; regionStart, regionEnd];
    end
end

% Between merged ranges
for i = 1:numMerged - 1
    gapStart = mergedRanges(i, 2) + mergedExclusion(i, 2);
    gapEnd = mergedRanges(i + 1, 1) - mergedExclusion(i + 1, 1);
    if gapEnd > gapStart + 0.5  % Minimum 0.5 Da width
        bgRegionsMc = [bgRegionsMc; gapStart, gapEnd];
    end
end

% After last range - this is often a good clean region
lastExcl = mergedExclusion(end, 2);
if mergedRanges(end, 2) + lastExcl < mcRange(2)
    regionStart = mergedRanges(end, 2) + lastExcl;
    regionEnd = mcRange(2);
    if regionEnd > regionStart + 0.1
        bgRegionsMc = [bgRegionsMc; regionStart, regionEnd];
    end
end

% Warn if very few background regions found
if isempty(bgRegionsMc)
    warning('tofSpecBackgroundDetermination:noBgRegions', 'No background regions found between ranges. Try reducing minPeakDistance.');
    % Fallback: use regions at very low and very high m/c
    bgRegionsMc = [max(1, mcRange(1)), min(mcRange(1) + 3, mergedRanges(1,1) - 0.5); mergedRanges(end, 2) + 1, mcRange(2)];
elseif size(bgRegionsMc, 1) < 3
    warning('tofSpecBackgroundDetermination:fewBgRegions', 'Only %d background region(s) found. Results may be less reliable.', size(bgRegionsMc, 1));
end

% Print background regions for debugging
fprintf('Background regions (m/c):\n');
for i = 1:size(bgRegionsMc, 1)
    fprintf('  [%.2f - %.2f] Da (width: %.2f)\n', bgRegionsMc(i,1), bgRegionsMc(i,2), bgRegionsMc(i,2)-bgRegionsMc(i,1));
end

% Convert to TOF: t = sqrt(m/c / k)
bgRegionsTof = sqrt(bgRegionsMc / k);

end


function plotBackgroundAnalysis(bgResult, info, pos, tofData, opts)
% Create diagnostic plots for TOF background analysis
%
% Creates a comprehensive figure with 6 subplots for troubleshooting:
% 1. Background rate over experiment (temporal stability)
% 2. TOF histogram with background regions highlighted
% 3. Mass spectrum with background regions and fit comparison
% 4. Zoomed background regions with local statistics
% 5. Residual analysis (observed - model)
% 6. Background level comparison (fitted vs alternatives)

figure('Name', 'TOF Background Analysis', 'Position', [100, 100, 1400, 900]);
regionStats = [];

% Get mass spectrum data
hasMc = ismember('mc', pos.Properties.VariableNames);
if hasMc
    mcEdges = opts.mcRange(1):opts.mcBinWidth:opts.mcRange(2);
    [counts, ~] = histcounts(pos.mc, mcEdges);
    mcCenters = mcEdges(1:end-1) + opts.mcBinWidth/2;
    if ismember('deltaP', pos.Properties.VariableNames) && strcmpi(info.normalization, 'pulse')
        totalNorm = sum(pos.deltaP, 'omitnan');
    else
        totalNorm = height(pos);
    end
    countsNorm = counts / (opts.mcBinWidth * totalNorm);  % counts per Da per normalization
end

mcValues = bgResult.mcCurve.mc;
bgPerBin = bgResult.mcCurve.mean;  % background per Da per normalization

% Convert background regions from TOF to m/c for plotting
k = info.tofConversion;
if isempty(k)
    k = 5e-4;
end
bgRegionsTofPlot = info.bgRegionsTofPlot;
if isempty(bgRegionsTofPlot)
    bgRegionsTofPlot = info.bgRegionsTof;
end
if isempty(bgRegionsTofPlot) && ismember('mc', pos.Properties.VariableNames) && ~isempty(opts.ranges)
    % Fallback: recompute plot regions from ranges for visualization
    [bgRegionsTofPlot, ~] = autoDetectBgRegions(pos.mc, opts.ranges, opts.minPeakDistance, ...
        k, opts.mcRange, opts.clusterThreshold, opts.clusterExtraDistance);
end
bgRegionsMc = bgRegionsTofPlot.^2 * k;  % m/c = k * t^2

%% Subplot 1: Background rate over experiment
subplot(2, 3, 1);
errorbar(bgResult.blockCenter / 1e6, bgResult.tofRate, ...
    bgResult.tofRateError, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
hold on;
if info.numBlocks > 1
    xFit = [bgResult.blockCenter(1), bgResult.blockCenter(end)] / 1e6;
    yFit = bgResult.meanRate + bgResult.rateTrend * ...
        ([bgResult.blockCenter(1), bgResult.blockCenter(end)] - mean(bgResult.blockCenter));
    plot(xFit, yFit, 'r--', 'LineWidth', 1.5);
end
yline(bgResult.meanRate, 'k:', 'LineWidth', 1);
hold off;
xlabel('Ion index [10^6 ions]');
if strcmpi(info.normalization, 'pulse')
    ylabel('Background rate [cts/ns/pulse]');
else
    ylabel('Background rate [cts/ns/totalIons]');
end
title('1. Rate Over Experiment');
legend('Measured', 'Trend', 'Mean', 'Location', 'best');
grid on;

%% Subplot 2: TOF histogram with background regions
subplot(2, 3, 2);
tofEdges = min(tofData):opts.tofBinWidth:max(tofData);
if ismember('deltaP', pos.Properties.VariableNames) && strcmpi(info.normalization, 'pulse')
    totalNorm = sum(pos.deltaP, 'omitnan');
else
    totalNorm = height(pos);
end
[counts, edges] = histcounts(tofData, tofEdges);
tofCenters = edges(1:end-1) + opts.tofBinWidth/2;
scale = opts.tofBinWidth * totalNorm;
if isfinite(scale) && scale > 0
    counts = counts / scale;
end
barHandle = bar(tofCenters, counts, 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
barHandle.DisplayName = 'TOF spectrum';
hold on;
% Use log scale with a positive lower bound so patches render correctly
set(gca, 'YScale', 'log');
minPositive = min(counts(counts > 0));
if isempty(minPositive) || ~isfinite(minPositive)
    minPositive = 1e-6;
end
ylim([minPositive, max(counts) * 1.2 + minPositive]);
ylims = ylim;
    if ~isempty(bgRegionsTofPlot)
        for r = 1:size(bgRegionsTofPlot, 1)
            tStart = bgRegionsTofPlot(r,1);
            tEnd = bgRegionsTofPlot(r,2);
            if tEnd <= tStart
                continue;
            end
            patchHandle = patch([tStart, tEnd, tEnd, tStart], ...
                [ylims(1), ylims(1), ylims(2), ylims(2)], 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'g');
            if r == 1
                patchHandle.DisplayName = 'BG regions';
            else
                patchHandle.HandleVisibility = 'off';
            end
        end
    else
        plot(nan, nan, 'g', 'LineWidth', 2, 'DisplayName', 'BG regions (n/a)');
    end
% Draw horizontal line for average background level in TOF space
% Background per TOF bin = tofRate * binWidth * total normalization
if ismember('deltaP', pos.Properties.VariableNames) && strcmpi(info.normalization, 'pulse')
    totalNorm = sum(pos.deltaP, 'omitnan');
else
    totalNorm = height(pos);
end
bgCountsPerTofBin = bgResult.meanRate * opts.tofBinWidth * totalNorm;
if isfinite(opts.tofBinWidth) && opts.tofBinWidth > 0 && isfinite(totalNorm) && totalNorm > 0
    bgRateLine = bgCountsPerTofBin / (opts.tofBinWidth * totalNorm);
else
    bgRateLine = bgResult.meanRate;
end
bgLine = yline(bgRateLine, 'r-', 'LineWidth', 2, 'Label', sprintf('BG = %.3g cts/ns', bgRateLine));
bgLine.DisplayName = 'BG rate';
% Also show suggested (corrected) background if different
if isfield(bgResult, 'correctionFactor') && abs(bgResult.correctionFactor - 1) > 0.05
    bgCorrected = bgRateLine * bgResult.correctionFactor;
    bgSuggested = yline(bgCorrected, 'b--', 'LineWidth', 2, 'Label', sprintf('Suggested = %.3g', bgCorrected));
    bgSuggested.DisplayName = 'Suggested';
end
hold off;
xlabel('Time-of-flight [ns]');
if strcmpi(info.normalization, 'pulse')
    ylabel('Counts/ns/pulse');
else
    ylabel('Counts/ns/totalIons');
end
title('2. TOF Spectrum + BG Regions');
legend('show', 'Location', 'best');

%% Subplot 3: Mass spectrum with background overlay and regions
subplot(2, 3, 3);
if hasMc
    spectrumLine = semilogy(mcCenters, countsNorm, 'k-', 'LineWidth', 0.5);
    spectrumLine.DisplayName = 'Spectrum';
    hold on;

    % Background model (fitted)
    semilogy(mcValues, bgPerBin, 'r-', 'LineWidth', 2, 'DisplayName', 'Model (fitted)');

    % Suggested correction based on validation
    if isfield(bgResult, 'correctionFactor') && isfinite(bgResult.correctionFactor) ...
            && abs(bgResult.correctionFactor - 1) > 0.02
        semilogy(mcValues, bgPerBin * bgResult.correctionFactor, 'b--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Model x%.2f (suggested)', bgResult.correctionFactor));
    end

    % Highlight background regions
    ylims = ylim;
    for r = 1:size(bgRegionsMc, 1)
        mcStart = max(bgRegionsMc(r,1), opts.mcRange(1));
        mcEnd = min(bgRegionsMc(r,2), opts.mcRange(2));
        if mcEnd > mcStart
            patch([mcStart, mcEnd, mcEnd, mcStart], ...
                [ylims(1), ylims(1), ylims(2), ylims(2)], 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'g', ...
                'HandleVisibility', 'off');
        end
    end
    hold off;

    xlabel('Mass-to-charge [Da]');
    if strcmpi(info.normalization, 'pulse')
        ylabel('Counts/Da/pulse');
    else
        ylabel('Counts/Da/totalIons');
    end
    title('3. Mass Spectrum + Background');
    xlim(opts.mcRange);
    legend('show', 'Location', 'northeast');
end

%% Subplot 4: Zoomed view of background regions with statistics
subplot(2, 3, 4);
if hasMc && ~isempty(bgRegionsMc)
    hold on;

    regionStats = [];
    colors = lines(size(bgRegionsMc, 1));
    plot(nan, nan, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 8, 'DisplayName', 'Region data');
    plot(nan, nan, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2, 'DisplayName', 'Model');

    for r = 1:size(bgRegionsMc, 1)
        mcStart = max(bgRegionsMc(r, 1), opts.mcRange(1));
        mcEnd = min(bgRegionsMc(r, 2), opts.mcRange(2));
        if mcEnd <= mcStart + opts.mcBinWidth
            continue;
        end

        % Get data in this region
        inRegion = mcCenters >= mcStart & mcCenters <= mcEnd;
        if sum(inRegion) < 3
            continue;
        end

        regionMc = mcCenters(inRegion);
        regionCounts = countsNorm(inRegion);

        % Get model prediction for this region
        modelInRegion = mcValues >= mcStart & mcValues <= mcEnd;
        modelMc = mcValues(modelInRegion);
        modelBg = bgPerBin(modelInRegion);

        % Plot data points
        plot(regionMc, regionCounts, '.', 'Color', colors(r,:), 'MarkerSize', 8);

        % Plot model
        plot(modelMc, modelBg, '-', 'Color', colors(r,:), 'LineWidth', 2);

        % Calculate statistics for this region
        if ~isempty(regionCounts) && ~isempty(modelBg)
            % Interpolate model to data points
            modelAtData = interp1(modelMc, modelBg, regionMc, 'linear', 'extrap');

            meanObs = mean(regionCounts);
            meanModel = mean(modelAtData);
            ratio = meanObs / meanModel;

            regionStats = [regionStats; struct(...
                'mcStart', mcStart, 'mcEnd', mcEnd, ...
                'meanObs', meanObs, 'meanModel', meanModel, ...
                'ratio', ratio, 'numPoints', sum(inRegion))];
        end
    end
    hold off;

    xlabel('Mass-to-charge [Da]');
    if strcmpi(info.normalization, 'pulse')
        ylabel('Counts/Da/pulse');
    else
        ylabel('Counts/Da/totalIons');
    end
    title('4. Background Regions (zoom)');
    set(gca, 'YScale', 'log');
    grid on;

    % Add statistics text
    if ~isempty(regionStats)
        ratios = [regionStats.ratio];
        meanRatio = mean(ratios);
        text(0.02, 0.98, sprintf('Mean obs/model: %.2f', meanRatio), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'BackgroundColor', 'w');
    end
else
    xlabel('Mass-to-charge [Da]');
    if strcmpi(info.normalization, 'pulse')
        ylabel('Counts/Da/pulse');
    else
        ylabel('Counts/Da/totalIons');
    end
    title('4. Background Regions (zoom)');
    text(0.5, 0.5, 'No BG regions (unranged-only)', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'BackgroundColor', 'w');
end

%% Subplot 5: Residual analysis
subplot(2, 3, 5);
if hasMc && ~isempty(bgRegionsMc)
    hold on;

    allResiduals = [];
    allMc = [];
    plot(nan, nan, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 6, 'DisplayName', 'Obs / Model');

    for r = 1:size(bgRegionsMc, 1)
        mcStart = max(bgRegionsMc(r, 1), opts.mcRange(1));
        mcEnd = min(bgRegionsMc(r, 2), opts.mcRange(2));
        if mcEnd <= mcStart + opts.mcBinWidth
            continue;
        end

        inRegion = mcCenters >= mcStart & mcCenters <= mcEnd;
        if sum(inRegion) < 2
            continue;
        end

        regionMc = mcCenters(inRegion);
        regionCounts = countsNorm(inRegion);

        % Interpolate model to data
        modelAtData = interp1(mcValues, bgPerBin, regionMc, 'linear', 'extrap');

        % Residual as ratio (observed / model)
        residualRatio = regionCounts ./ modelAtData;

        allResiduals = [allResiduals; residualRatio(:)];
        allMc = [allMc; regionMc(:)];

        plot(regionMc, residualRatio, '.', 'MarkerSize', 6);
    end

    unityLine = yline(1, 'k-', 'LineWidth', 1.5);
    unityLine.DisplayName = 'Unity';
    hiLine = yline(1.5, 'r--', 'LineWidth', 1, 'Label', '1.5x');
    hiLine.DisplayName = '1.5x';
    loLine = yline(0.67, 'r--', 'LineWidth', 1, 'Label', '0.67x');
    loLine.DisplayName = '0.67x';
    hold off;

    xlabel('Mass-to-charge [Da]');
    ylabel('Observed / Model');
    title('5. Residual Ratio in BG Regions');
    ylim([0, max(3, max(allResiduals) * 1.1)]);
    grid on;

    % Add statistics
    if ~isempty(allResiduals)
        medianRatio = median(allResiduals);
        text(0.02, 0.98, sprintf('Median ratio: %.2f\nSuggested correction: x%.2f', ...
            medianRatio, medianRatio), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontSize', 9, 'BackgroundColor', 'w');
    end
else
    xlabel('Mass-to-charge [Da]');
    ylabel('Observed / Model');
    title('5. Residual Ratio in BG Regions');
    text(0.5, 0.5, 'No BG regions (unranged-only)', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'BackgroundColor', 'w');
end

%% Subplot 6: Background level comparison
subplot(2, 3, 6);
if hasMc
    % Show spectrum with multiple background levels for comparison
    specLine = semilogy(mcCenters, countsNorm, 'k-', 'LineWidth', 0.5);
    specLine.DisplayName = 'Spectrum';
    hold on;

    % Fitted background
    semilogy(mcValues, bgPerBin, 'r-', 'LineWidth', 2, 'DisplayName', sprintf('Fitted (%.4f)', bgResult.meanRate * 1e4));

    % Alternative levels for comparison
    scalingFactors = [1.25, 1.5, 1.75, 2.0];
    colors = [0.8 0.4 0; 0 0.5 0.8; 0.5 0 0.8; 0.8 0 0.4];

    for i = 1:length(scalingFactors)
        sf = scalingFactors(i);
        semilogy(mcValues, bgPerBin * sf, '--', 'Color', colors(i,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('x%.2f', sf));
    end

    if isfield(bgResult, 'correctionFactor') && isfinite(bgResult.correctionFactor) ...
            && abs(bgResult.correctionFactor - 1) > 0.02
        semilogy(mcValues, bgPerBin * bgResult.correctionFactor, 'b-', 'LineWidth', 2, ...
            'DisplayName', sprintf('Suggested x%.2f', bgResult.correctionFactor));
    end
    hold off;

    xlabel('Mass-to-charge [Da]');
    if strcmpi(info.normalization, 'pulse')
        ylabel('Counts/Da/pulse');
    else
        ylabel('Counts/Da/totalIons');
    end
    title('6. Background Level Comparison');
    xlim(opts.mcRange);
    legend('show', 'Location', 'northeast', 'FontSize', 8);

    % Add note about how to use
    text(0.02, 0.02, 'Compare lines to valleys to find best match', ...
        'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
        'FontSize', 8, 'FontAngle', 'italic');
end

%% Summary title
if ~isempty(regionStats)
    meanRatio = mean([regionStats.ratio]);
    suggestedRate = bgResult.meanRate * meanRatio;
    if strcmpi(info.normalization, 'pulse')
        unitLabel = 'counts/ns/pulse';
    else
        unitLabel = 'counts/ns/totalIons';
    end
    sgtitle(sprintf(['TOF Background Analysis\n' ...
        'Fitted rate: %.4e %s | Obs/Model ratio: %.2f | Suggested rate: %.4e'], ...
        bgResult.meanRate, unitLabel, meanRatio, suggestedRate), 'FontSize', 11);
else
    if strcmpi(info.normalization, 'pulse')
        sgtitle(sprintf('TOF Background Analysis - Mean rate: %.4e counts/ns/pulse', bgResult.meanRate));
    else
        sgtitle(sprintf('TOF Background Analysis - Mean rate: %.4e counts/ns/totalIons', bgResult.meanRate));
    end
end

%% Print diagnostic summary to command window
fprintf('\n=== TOF Background Diagnostic Summary ===\n');
if strcmpi(info.normalization, 'pulse')
fprintf('Fitted background rate: %.4e counts/ns/pulse\n', bgResult.meanRate);
else
    fprintf('Fitted background rate: %.4e counts/ns/totalIons\n', bgResult.meanRate);
end
fprintf('Number of background regions: %d\n', size(bgRegionsMc, 1));
if isfield(info, 'totalBgTime') && isfinite(info.totalBgTime)
    fprintf('Total TOF width used: %.1f ns\n', info.totalBgTime);
elseif isfield(info, 'totalBgTimePerBlock') && ~isempty(info.totalBgTimePerBlock)
    fprintf('Mean TOF width used: %.1f ns\n', mean(info.totalBgTimePerBlock, 'omitnan'));
else
    fprintf('Total TOF width used: n/a\n');
end

if ~isempty(regionStats)
    fprintf('\nPer-region statistics:\n');
    fprintf('  Region [Da]      Obs/Model   Points\n');
    for i = 1:length(regionStats)
        rs = regionStats(i);
        fprintf('  [%5.1f - %5.1f]    %.2f       %d\n', ...
            rs.mcStart, rs.mcEnd, rs.ratio, rs.numPoints);
    end

    meanRatio = mean([regionStats.ratio]);
    fprintf('\nOverall obs/model ratio: %.2f\n', meanRatio);
    fprintf('>>> Suggested correction: multiply background by %.2f <<<\n', meanRatio);
    fprintf('>>> Suggested rate: %.4e (vs fitted %.4e) <<<\n', ...
        bgResult.meanRate * meanRatio, bgResult.meanRate);
end
fprintf('==========================================\n\n');

end
