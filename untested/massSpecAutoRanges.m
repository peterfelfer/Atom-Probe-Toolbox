function [rangeTable, info] = massSpecAutoRanges(posOrMc, peaks, assigned, options)
% MASSSPECAUTORANGES Build non-overlapping ranges for assigned peaks.
%
% [rangeTable, info] = massSpecAutoRanges(posOrMc, peaks, assigned)
%
% INPUT
% posOrMc: pos table with mc column or numeric mc vector
% peaks:   table from massSpecFindPeaks (must contain mc)
% assigned: struct array from massSpecAutoAssignIons info.assigned
%
% OPTIONS
%   'binWidth'           - bin width in Da (default: 0.02)
%   'baselineQuantile'   - quantile for baseline (default: 0.1)
%   'baselineWindow'     - baseline window in Da (default: 0.5)
%   'smoothSpan'         - smoothing span in Da (default: 0.08)
%   'noiseWindow'        - noise window in Da (default: 0.5)
%   'noiseK'             - threshold multiplier (default: 3)
%   'spacingTolerance'   - relative tolerance for isotope spacing (default: 0.2)
%
% OUTPUT
% rangeTable: table with mcMin, mcMax, mcMid, ionName, ion, chargeState, peakMc
% info:       struct with spectrum arrays and settings
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    posOrMc
    peaks table
    assigned struct
    options.binWidth (1,1) double {mustBePositive} = 0.02
    options.baselineQuantile (1,1) double {mustBeGreaterThanOrEqual(options.baselineQuantile,0), mustBeLessThanOrEqual(options.baselineQuantile,1)} = 0.1
    options.baselineWindow (1,1) double {mustBePositive} = 0.5
    options.smoothSpan (1,1) double {mustBePositive} = 0.08
    options.noiseWindow (1,1) double {mustBePositive} = 0.5
    options.noiseK (1,1) double {mustBePositive} = 3
    options.spacingTolerance (1,1) double {mustBeNonnegative} = 0.2
end

if ~istable(peaks) || ~ismember('mc', peaks.Properties.VariableNames)
    error('massSpecAutoRanges:invalidPeaks', 'peaks must contain mc.');
end

% Resolve mc vector
if istable(posOrMc)
    if ~ismember('mc', posOrMc.Properties.VariableNames)
        error('massSpecAutoRanges:missingMc', 'pos table must contain mc column.');
    end
    mc = posOrMc.mc;
else
    mc = posOrMc;
end
mc = mc(:);
mc = mc(isfinite(mc));

% Spectrum
binWidth = options.binWidth;
mcMin = min(mc);
mcMax = max(mc);
edges = mcMin:binWidth:(ceil(mcMax / binWidth) * binWidth);
counts = histcounts(mc, edges);
centers = edges(1:end-1) + binWidth/2;

% Baseline + residual + noise
baselineBins = max(3, round(options.baselineWindow / binWidth));
if mod(baselineBins, 2) == 0
    baselineBins = baselineBins + 1;
end
baseline = localMovingQuantile(counts, options.baselineQuantile, baselineBins);
residual = counts - baseline;

smoothBins = max(3, round(options.smoothSpan / binWidth));
if mod(smoothBins, 2) == 0
    smoothBins = smoothBins + 1;
end
smoothResidual = smoothdata(residual, 'sgolay', smoothBins);

noiseBins = max(3, round(options.noiseWindow / binWidth));
if mod(noiseBins, 2) == 0
    noiseBins = noiseBins + 1;
end
med = movmedian(residual, noiseBins, 'omitnan');
noise = movmedian(abs(residual - med), noiseBins, 'omitnan') * 1.4826;
noise(noise == 0) = eps;

% Build peak list from assigned ions (matched peaks only)
rows = [];
for i = 1:numel(assigned)
    ionStr = assigned(i).ion;
    cs = assigned(i).chargeState;
    matchedMc = assigned(i).matchedMc(:);
    if isempty(matchedMc)
        continue;
    end
    for k = 1:numel(matchedMc)
        rows = [rows; {matchedMc(k), ionStr, cs, assigned(i).score}]; %#ok<AGROW>
    end
end

if isempty(rows)
    rangeTable = table();
    info = struct();
    return;
end

peakMc = cell2mat(rows(:,1));
ionStr = rows(:,2);
chargeState = cell2mat(rows(:,3));
score = cell2mat(rows(:,4));

% Map to nearest detected peaks for heights
peakIdx = zeros(size(peakMc));
peakHeight = zeros(size(peakMc));
for i = 1:numel(peakMc)
    [~, idx] = min(abs(peaks.mc - peakMc(i)));
    peakIdx(i) = idx;
    if ismember('prominence', peaks.Properties.VariableNames)
        peakHeight(i) = peaks.prominence(idx);
    elseif ismember('height', peaks.Properties.VariableNames)
        peakHeight(i) = peaks.height(idx);
    else
        peakHeight(i) = 1;
    end
end

% Enforce isotope spacing per charge state (drop weaker peaks if too close)
keep = true(size(peakMc));
for i = 1:numel(unique(ionStr))
    ionName = ionStr{i};
    idxIon = find(strcmp(ionStr, ionName));
    if isempty(idxIon)
        continue;
    end
    cs = chargeState(idxIon(1));
    minSep = 1 / cs;
    minAllowed = minSep * (1 - options.spacingTolerance);
    [mcSorted, order] = sort(peakMc(idxIon));
    idxSorted = idxIon(order);
    for k = 1:(numel(idxSorted)-1)
        if mcSorted(k+1) - mcSorted(k) < minAllowed
            a = idxSorted(k);
            b = idxSorted(k+1);
            if peakHeight(a) >= peakHeight(b)
                keep(b) = false;
            else
                keep(a) = false;
            end
        end
    end
end

peakMc = peakMc(keep);
ionStr = ionStr(keep);
chargeState = chargeState(keep);
score = score(keep);
peakIdx = peakIdx(keep);

% Sort peaks by mc
[peakMc, order] = sort(peakMc);
ionStr = ionStr(order);
chargeState = chargeState(order);
score = score(order);
peakIdx = peakIdx(order);

% Determine boundaries (non-overlapping)
peakBins = zeros(size(peakMc));
for i = 1:numel(peakMc)
    [~, peakBins(i)] = min(abs(centers - peakMc(i)));
end

boundaries = nan(numel(peakMc)+1,1);
% left boundary for first peak
boundaries(1) = findLeftBoundary(peakBins(1), smoothResidual, noise, options.noiseK);
% middle boundaries
for i = 1:(numel(peakMc)-1)
    leftBin = peakBins(i);
    rightBin = peakBins(i+1);
    if rightBin <= leftBin
        midBin = leftBin;
    else
        [~, midRel] = min(smoothResidual(leftBin:rightBin));
        midBin = leftBin + midRel - 1;
    end
    boundaries(i+1) = midBin;
end
% right boundary for last peak
boundaries(end) = findRightBoundary(peakBins(end), smoothResidual, noise, options.noiseK);

mcMinRanges = zeros(numel(peakMc),1);
mcMaxRanges = zeros(numel(peakMc),1);
for i = 1:numel(peakMc)
    lbin = max(1, boundaries(i));
    rbin = min(numel(centers), boundaries(i+1));
    if rbin <= lbin
        lbin = max(1, peakBins(i)-1);
        rbin = min(numel(centers), peakBins(i)+1);
    end
    mcMinRanges(i) = centers(lbin);
    mcMaxRanges(i) = centers(rbin);
end

mcMidRanges = (mcMinRanges + mcMaxRanges) / 2;
ionName = categorical(string(ionStr));
ion = ionStr;

rangeTable = table(mcMinRanges, mcMaxRanges, mcMidRanges, ionName, ion, chargeState, peakMc, score, ...
    'VariableNames', {'mcMin','mcMax','mcMid','ionName','ion','chargeState','peakMc','score'});

info = struct();
info.centers = centers;
info.counts = counts;
info.baseline = baseline;
info.residual = residual;
info.smoothResidual = smoothResidual;
info.noise = noise;
info.settings = options;
end

function q = localMovingQuantile(x, quantile, windowSize)
    n = numel(x);
    half = floor(windowSize / 2);
    q = zeros(size(x));
    for i = 1:n
        i1 = max(1, i - half);
        i2 = min(n, i + half);
        q(i) = prctile(x(i1:i2), quantile * 100);
    end
end

function lbin = findLeftBoundary(peakBin, smoothResidual, noise, noiseK)
    lbin = peakBin;
    for i = peakBin:-1:1
        if smoothResidual(i) <= noiseK * noise(i)
            lbin = i;
            break;
        end
    end
end

function rbin = findRightBoundary(peakBin, smoothResidual, noise, noiseK)
    rbin = peakBin;
    for i = peakBin:numel(smoothResidual)
        if smoothResidual(i) <= noiseK * noise(i)
            rbin = i;
            break;
        end
    end
end
