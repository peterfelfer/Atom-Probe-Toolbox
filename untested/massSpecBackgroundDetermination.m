function [bgResult, info] = massSpecBackgroundDetermination(pos, bgRegions, varargin)
% massSpecBackgroundDetermination Determine background in mass spectrum and track over experiment.
%
% Analyzes background in mass-to-charge (m/c) space and tracks how it
% changes over the course of the experiment using block-based statistics.
%
% SYNTAX
%   [bgResult, info] = massSpecBackgroundDetermination(pos, bgRegions)
%   [bgResult, info] = massSpecBackgroundDetermination(pos, bgRegions, Name, Value, ...)
%
% INPUT
%   pos         - Table with atom probe data. Must contain 'mc' column.
%                 Should contain 'ionIdx' for proper event ordering.
%
%   bgRegions   - Background regions specified as:
%                 - Nx2 array of [mcStart, mcEnd] pairs in Da, OR
%                 - 'auto': automatically detect valleys (requires 'ranges')
%                 - struct with field 'mc' for explicit specification
%
% NAME-VALUE PAIRS
%   'blockSize'       - Number of ions per block (default: [])
%   'numBlocks'       - Alternative: specify number of blocks (overrides blockSize)
%   'mcBinWidth'      - Bin width for m/c histogram [Da] (default: 0.01)
%   'ranges'          - Range table (for 'auto' background detection)
%   'minPeakDistance' - Minimum distance from isolated peaks [Da] (default: 0.3)
%   'clusterThreshold' - Ranges closer than this are considered a cluster [Da] (default: 2.0)
%   'clusterExtraDistance' - Extra distance to add around peak clusters [Da] (default: 0.5)
%   'fitType'         - 'invSqrt' (default: 'invSqrt')
%   'normalization'   - 'auto' | 'ion' | 'pulse' (default: 'auto')
%   'plotResult'      - Logical, create diagnostic plots (default: false)
%   'mcRange'         - [mcMin, mcMax] range for output background curve (default: [0, 200])
%
% OUTPUT
%   bgResult    - Table containing:
%                 .ionIdx      - Center ion index of each block
%                 .coeff       - invSqrt coefficient per block
%
%   info        - Structure containing:
%                 .bgRegionsMc   - Background regions in m/c space
%                 .totalBgWidth  - Total m/c width used [Da]
%                 .numBlocks     - Number of blocks analyzed
%                 .normalization - 'ion' or 'pulse'
%                 .fitType       - 'invSqrt'
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Parse inputs
p = inputParser;
addRequired(p, 'pos', @istable);
addRequired(p, 'bgRegions');
addParameter(p, 'blockSize', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'numBlocks', [], @(x) isempty(x) || (isnumeric(x) && x > 0));
addParameter(p, 'mcBinWidth', 0.01, @(x) isnumeric(x) && x > 0);
addParameter(p, 'ranges', [], @(x) isempty(x) || istable(x));
addParameter(p, 'minPeakDistance', 0.3, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'clusterThreshold', 2.0, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'clusterExtraDistance', 0.5, @(x) isnumeric(x) && x >= 0);
addParameter(p, 'fitType', 'invSqrt', @(x) ischar(x) || isstring(x));
addParameter(p, 'normalization', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'plotResult', false, @islogical);
addParameter(p, 'mcRange', [0, 200], @(x) isnumeric(x) && numel(x) == 2);

parse(p, pos, bgRegions, varargin{:});
opts = p.Results;

if ~ismember('mc', pos.Properties.VariableNames)
    error('massSpecBackgroundDetermination:missingMc', ...
        'pos table must contain mc column.');
end

hasIonIdx = ismember('ionIdx', pos.Properties.VariableNames);
hasDeltaP = ismember('deltaP', pos.Properties.VariableNames);

% Determine normalization basis (enforce pulse if available)
normMode = lower(string(opts.normalization));
if hasDeltaP
    if normMode ~= "pulse"
        warning('massSpecBackgroundDetermination:forcePulseNormalization', ...
            'deltaP found. Enforcing pulse normalization.');
    end
    normMode = "pulse";
else
    if normMode == "" || normMode == "auto"
        normMode = "ion";
    end
    if normMode == "pulse"
        warning('massSpecBackgroundDetermination:noDeltaP', ...
            'Normalization by pulse requested but deltaP column not found. Using ion normalization.');
        normMode = "ion";
    end
end

numIons = height(pos);
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

% Event ordering
if hasIonIdx
    [~, sortOrder] = sort(pos.ionIdx);
else
    warning('massSpecBackgroundDetermination:noIonIdx', ...
        'No ionIdx column found. Using row order as event sequence.');
    sortOrder = (1:height(pos))';
end

%% Process background regions
if ischar(bgRegions) && strcmpi(bgRegions, 'auto')
    if isempty(opts.ranges)
        error('massSpecBackgroundDetermination:noRanges', ...
            'Automatic background detection requires ''ranges'' parameter.');
    end
    [bgRegionsMc, opts.mcRange] = autoDetectBgRegionsMc(pos.mc, opts.ranges, opts.minPeakDistance, ...
        opts.mcRange, opts.clusterThreshold, opts.clusterExtraDistance);
elseif isstruct(bgRegions)
    if isfield(bgRegions, 'mc')
        bgRegionsMc = bgRegions.mc;
    else
        error('massSpecBackgroundDetermination:invalidRegions', ...
            'bgRegions struct must have ''mc'' field.');
    end
else
    bgRegionsMc = bgRegions;
end

info.bgRegionsMc = bgRegionsMc;
if isempty(bgRegionsMc)
    error('massSpecBackgroundDetermination:noBgRegions', ...
        'No background regions specified or detected.');
end

totalBgWidth = sum(bgRegionsMc(:,2) - bgRegionsMc(:,1));
info.totalBgWidth = totalBgWidth;

%% Determine block boundaries
if ~isempty(opts.numBlocks)
    numBlocks = opts.numBlocks;
    blockSize = ceil(numIons / numBlocks);
elseif ~isempty(opts.blockSize)
    blockSize = opts.blockSize;
    numBlocks = ceil(numIons / blockSize);
else
    numBlocks = 1;
    blockSize = numIons;
end
info.numBlocks = numBlocks;

blockIdx = zeros(numBlocks, 2);
for b = 1:numBlocks
    blockIdx(b, 1) = (b-1) * blockSize + 1;
    blockIdx(b, 2) = min(b * blockSize, numIons);
end

%% Analyze background per block
fitType = lower(string(opts.fitType));
if fitType == "1/sqrt" || fitType == "onesqrt" || fitType == "inv"
    fitType = "invsqrt";
end
if fitType ~= "invsqrt"
    warning('massSpecBackgroundDetermination:fitTypeOverride', ...
        'Only invSqrt fit is supported. Using invSqrt.');
    fitType = "invsqrt";
end
bgLevel = zeros(numBlocks, 1);
bgCoeff = zeros(numBlocks, 1);
fitError = zeros(numBlocks, 1);

mcValues = (opts.mcRange(1) + opts.mcBinWidth/2 : opts.mcBinWidth : opts.mcRange(2))';
mcCurve = zeros(numel(mcValues), numBlocks);

for b = 1:numBlocks
    blockIons = sortOrder(blockIdx(b,1):blockIdx(b,2));
    blockMc = pos.mc(blockIons);
    numBlockIons = numel(blockIons);

    if normMode == "pulse"
        numBlockNorm = sum(pos.deltaP(blockIons), 'omitnan');
    else
        numBlockNorm = numBlockIons;
    end
    if ~isfinite(numBlockNorm) || numBlockNorm <= 0
        numBlockNorm = numBlockIons;
    end

    mcEdges = opts.mcRange(1):opts.mcBinWidth:opts.mcRange(2);
    [counts, ~] = histcounts(blockMc, mcEdges);
    mcCenters = mcEdges(1:end-1) + opts.mcBinWidth/2;
    countsNorm = counts / (opts.mcBinWidth * numBlockNorm); % counts/Da/normalization

    inBg = false(size(mcCenters));
    for r = 1:size(bgRegionsMc, 1)
        inBg = inBg | (mcCenters >= bgRegionsMc(r,1) & mcCenters <= bgRegionsMc(r,2));
    end

    if ~any(inBg)
    bgLevel(b) = NaN;
    bgCoeff(b, 1) = NaN;
    mcCurve(:, b) = NaN;
    fitError(b) = NaN;
        continue;
    end

    if fitType == "linear"
        pfit = polyfit(mcCenters(inBg), countsNorm(inBg), 1);
        bgCoeff(b, 1) = pfit(1);
        bgLevel(b) = mean(countsNorm(inBg));
        mcCurve(:, b) = polyval(pfit, mcValues);
        fitError(b) = std(countsNorm(inBg) - polyval(pfit, mcCenters(inBg)));
    elseif fitType == "invsqrt"
        xfit = 1 ./ sqrt(mcCenters(inBg));
        yfit = countsNorm(inBg);
        denom = sum(xfit.^2);
        if denom > 0
            a = sum(xfit .* yfit) / denom;
        else
            a = 0;
        end
        bgCoeff(b, 1) = a;
        bgLevel(b) = mean(yfit);
        mcCurve(:, b) = a ./ sqrt(mcValues);
        fitError(b) = std(yfit - a * xfit);
    else
        bgLevel(b) = mean(countsNorm(inBg));
        bgCoeff(b, 1) = bgLevel(b);
        mcCurve(:, b) = bgLevel(b);
        fitError(b) = std(countsNorm(inBg) - bgLevel(b));
    end
end

blockCenter = mean(blockIdx, 2);
bgResult = table(blockCenter, bgCoeff, 'VariableNames', {'ionIdx', 'coeff'});
if strcmpi(normMode, 'pulse')
    coeffUnit = 'counts/Da/pulse * sqrt(Da)';
else
    coeffUnit = 'counts/Da/totalIons * sqrt(Da)';
end
bgResult.Properties.VariableUnits = {'ion index', coeffUnit};

validLevels = bgLevel(isfinite(bgLevel));
meanLevel = mean(validLevels, 'omitnan');
levelStd = std(validLevels, 'omitnan');
if numBlocks > 1
    p = polyfit(blockCenter, bgLevel, 1);
    levelTrend = p(1);
else
    levelTrend = 0;
end

mcCurveTable = array2table(mcValues, 'VariableNames', {'mc'});
for b = 1:numBlocks
    mcCurveTable.(['block' num2str(b)]) = mcCurve(:, b);
end
mcCurveTable.mean = mean(mcCurve, 2, 'omitnan');

info.fitType = char(fitType);
info.bgLevel = bgLevel;
info.bgCoeff = bgCoeff;
info.blockIdx = blockIdx;
info.blockCenter = blockCenter;
info.fitError = fitError;
info.normalization = char(normMode);
info.totalNormalization = totalNorm;
info.meanLevel = meanLevel;
info.levelStd = levelStd;
info.levelTrend = levelTrend;
info.mcCurve = mcCurveTable;

%% Plotting
if opts.plotResult
    plotMassSpecBackground(info, pos, bgRegionsMc, opts);
end

end

function [bgRegionsMc, mcRange] = autoDetectBgRegionsMc(mc, ranges, minDist, mcRange, clusterThreshold, clusterExtraDistance)
% Auto-detect background regions in m/c space from ranges

if nargin < 5 || isempty(clusterThreshold)
    clusterThreshold = 2.0;
end
if nargin < 6 || isempty(clusterExtraDistance)
    clusterExtraDistance = 0.5;
end

actualMax = max(mc);
actualMin = min(mc);
if mcRange(2) < actualMax
    mcRange(2) = actualMax + 5;
end
if mcRange(1) > actualMin
    mcRange(1) = max(0, actualMin - 1);
end

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

leftExclusion = zeros(numRanges, 1);
rightExclusion = ones(numRanges, 1) * minDist;
for i = 1:numRanges
    inCluster = (leftDist(i) < clusterThreshold) || (rightDist(i) < clusterThreshold);
    if inCluster
        rightExclusion(i) = minDist + clusterExtraDistance;
    end
end

% Merge overlapping/touching ranges
mergedRanges = [];
mergedExclusion = [];
currentStart = rangeMatrix(1, 1);
currentEnd = rangeMatrix(1, 2);
currentLeftExcl = leftExclusion(1);
currentRightExcl = rightExclusion(1);

for i = 2:numRanges
    nextStart = rangeMatrix(i, 1);
    nextEnd = rangeMatrix(i, 2);
    if nextStart <= currentEnd + 0.01
        currentEnd = max(currentEnd, nextEnd);
        currentRightExcl = rightExclusion(i);
    else
        mergedRanges = [mergedRanges; currentStart, currentEnd];
        mergedExclusion = [mergedExclusion; currentLeftExcl, currentRightExcl];
        currentStart = nextStart;
        currentEnd = nextEnd;
        currentLeftExcl = leftExclusion(i);
        currentRightExcl = rightExclusion(i);
    end
end
mergedRanges = [mergedRanges; currentStart, currentEnd];
mergedExclusion = [mergedExclusion; currentLeftExcl, currentRightExcl];

bgRegionsMc = [];
numMerged = size(mergedRanges, 1);

% Before first range
firstExcl = mergedExclusion(1, 1);
if mergedRanges(1, 1) - firstExcl > mcRange(1)
    regionStart = max(mcRange(1), 0.1);
    regionEnd = mergedRanges(1, 1) - firstExcl;
    if regionEnd > regionStart + 0.1
        bgRegionsMc = [bgRegionsMc; regionStart, regionEnd];
    end
end

% Between merged ranges
for i = 1:numMerged - 1
    gapStart = mergedRanges(i, 2) + mergedExclusion(i, 2);
    gapEnd = mergedRanges(i + 1, 1) - mergedExclusion(i + 1, 1);
    if gapEnd > gapStart + 0.5
        bgRegionsMc = [bgRegionsMc; gapStart, gapEnd];
    end
end

% After last range
lastExcl = mergedExclusion(end, 2);
if mergedRanges(end, 2) + lastExcl < mcRange(2)
    regionStart = mergedRanges(end, 2) + lastExcl;
    regionEnd = mcRange(2);
    if regionEnd > regionStart + 0.1
        bgRegionsMc = [bgRegionsMc; regionStart, regionEnd];
    end
end

if isempty(bgRegionsMc)
    warning('massSpecBackgroundDetermination:noBgRegions', ...
        'No background regions found between ranges. Try reducing minPeakDistance.');
end
end

function plotMassSpecBackground(info, pos, bgRegionsMc, opts)
figure('Name', 'Mass Spectrum Background Analysis', 'Position', [120, 120, 1200, 700]);

hasMc = ismember('mc', pos.Properties.VariableNames);
if ~hasMc
    return;
end

mcEdges = opts.mcRange(1):opts.mcBinWidth:opts.mcRange(2);
[counts, ~] = histcounts(pos.mc, mcEdges);
mcCenters = mcEdges(1:end-1) + opts.mcBinWidth/2;
if strcmpi(info.normalization, 'pulse') && ismember('deltaP', pos.Properties.VariableNames)
    totalNorm = sum(pos.deltaP, 'omitnan');
else
    totalNorm = height(pos);
end
countsNorm = counts / (opts.mcBinWidth * totalNorm);

%% Subplot 1: Background level over experiment
subplot(2, 2, 1);
bgLevelPlot = info.bgLevel;
minPositive = min(bgLevelPlot(bgLevelPlot > 0));
if isempty(minPositive) || ~isfinite(minPositive)
    minPositive = eps;
end
bgLevelPlot(~isfinite(bgLevelPlot) | bgLevelPlot <= 0) = minPositive;
semilogy(info.blockCenter / 1e6, bgLevelPlot, 'o-', 'LineWidth', 1.5);
hold on;
meanLevel = info.meanLevel;
if ~isfinite(meanLevel) || meanLevel <= 0
    meanLevel = minPositive;
end
yline(meanLevel, 'k:', 'LineWidth', 1);
hold off;
xlabel('Ion index [10^6 ions]');
if strcmpi(info.normalization, 'pulse')
    ylabel('Background [cts/Da/pulse]');
else
    ylabel('Background [cts/Da/totalIons]');
end
title('1. Background Level Over Experiment');
legend('Measured', 'Mean', 'Location', 'best');
grid on;

%% Subplot 2: Mass spectrum with background regions
subplot(2, 2, 2);
semilogy(mcCenters, countsNorm, 'k-', 'LineWidth', 0.5, 'DisplayName', 'Spectrum');
hold on;
bgLine = info.mcCurve.mean;
semilogy(info.mcCurve.mc, bgLine, 'r-', 'LineWidth', 2, 'DisplayName', 'Background (mean)');
ylims = ylim;
for r = 1:size(bgRegionsMc, 1)
    patch([bgRegionsMc(r,1), bgRegionsMc(r,2), bgRegionsMc(r,2), bgRegionsMc(r,1)], ...
        [ylims(1), ylims(1), ylims(2), ylims(2)], 'g', 'FaceAlpha', 0.15, 'EdgeColor', 'g', ...
        'HandleVisibility', 'off');
end
hold off;
xlabel('Mass-to-charge [Da]');
if strcmpi(info.normalization, 'pulse')
    ylabel('Counts/Da/pulse');
else
    ylabel('Counts/Da/totalIons');
end
title('2. Mass Spectrum + Background');
legend('show', 'Location', 'northeast');
grid on;

%% Subplot 3: Background regions zoom + stats
subplot(2, 2, 3);
hold on;
colors = lines(size(bgRegionsMc, 1));
regionStats = [];
for r = 1:size(bgRegionsMc, 1)
    inRegion = mcCenters >= bgRegionsMc(r,1) & mcCenters <= bgRegionsMc(r,2);
    if sum(inRegion) < 3
        continue;
    end
    plot(mcCenters(inRegion), countsNorm(inRegion), '.', 'Color', colors(r,:), 'MarkerSize', 8);
    modelAtData = interp1(info.mcCurve.mc, bgLine, mcCenters(inRegion), 'linear', 'extrap');
    plot(mcCenters(inRegion), modelAtData, '-', 'Color', colors(r,:), 'LineWidth', 1.5);
    meanObs = mean(countsNorm(inRegion));
    meanModel = mean(modelAtData);
    regionStats = [regionStats; meanObs / meanModel]; %#ok<AGROW>
end
hold off;
set(gca, 'YScale', 'log');
xlabel('Mass-to-charge [Da]');
if strcmpi(info.normalization, 'pulse')
    ylabel('Counts/Da/pulse');
else
    ylabel('Counts/Da/totalIons');
end
title('3. Background Regions (zoom)');
grid on;
if ~isempty(regionStats)
    text(0.02, 0.98, sprintf('Mean obs/model: %.2f', mean(regionStats)), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'FontSize', 9, 'BackgroundColor', 'w');
end

%% Subplot 4: Residual ratio
subplot(2, 2, 4);
hold on;
for r = 1:size(bgRegionsMc, 1)
    inRegion = mcCenters >= bgRegionsMc(r,1) & mcCenters <= bgRegionsMc(r,2);
    if sum(inRegion) < 3
        continue;
    end
    modelAtData = interp1(info.mcCurve.mc, bgLine, mcCenters(inRegion), 'linear', 'extrap');
    plot(mcCenters(inRegion), countsNorm(inRegion) ./ modelAtData, '.', 'MarkerSize', 6);
end
yline(1, 'k-', 'LineWidth', 1.5);
hold off;
xlabel('Mass-to-charge [Da]');
ylabel('Observed / Model');
title('4. Residual Ratio in BG Regions');
grid on;
end
