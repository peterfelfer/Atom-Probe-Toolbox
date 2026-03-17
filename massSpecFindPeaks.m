function [peakTable, info] = massSpecFindPeaks(posOrMc, options)
% MASSSPECFINDPEAKS Data-driven peak detection for APT mass spectra.
%
% [peakTable, info] = massSpecFindPeaks(posOrMc)
% [peakTable, info] = massSpecFindPeaks(posOrMc, 'binWidth', 0.02)
%
% INPUT:
%   posOrMc   - pos table with mc column or numeric mc vector (Da)
%
% OPTIONS:
%   'binWidth'           - Histogram bin width in Da (default: 0.02)
%   'mcRange'            - [min max] range in Da (default: data range)
%   'usePulseNormalization' - Normalize by total pulses if deltaP exists (default: true)
%   'baselineQuantile'   - Quantile for baseline estimate (default: 0.1)
%   'baselineWindow'     - Window in Da for baseline (default: 0.5)
%   'smoothSpan'         - Smoothing span in Da (default: 0.1)
%   'noiseWindow'        - Window in Da for noise estimate (default: 0.5)
%   'minProminenceFactor'- Min prominence = factor * local noise (default: 6)
%   'minWidthBins'       - Minimum peak width in bins (default: 2)
%   'minDistanceDa'      - Minimum peak separation in Da (default: 0.04)
%   'minPeakHeight'      - Minimum peak height (default: 0)
%
% OUTPUT:
%   peakTable - table with detected peak properties
%   info      - struct with spectrum arrays and settings
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    posOrMc
    options.binWidth (1,1) double {mustBePositive} = 0.02
    options.mcRange (1,2) double = [NaN NaN]
    options.usePulseNormalization (1,1) logical = true
    options.baselineQuantile (1,1) double {mustBeGreaterThanOrEqual(options.baselineQuantile,0), mustBeLessThanOrEqual(options.baselineQuantile,1)} = 0.1
    options.baselineWindow (1,1) double {mustBePositive} = 0.5
    options.smoothSpan (1,1) double {mustBePositive} = 0.1
    options.noiseWindow (1,1) double {mustBePositive} = 0.5
    options.noiseFloorQuantile (1,1) double {mustBeGreaterThanOrEqual(options.noiseFloorQuantile,0), mustBeLessThanOrEqual(options.noiseFloorQuantile,1)} = 0.2
    options.minProminenceFactor (1,1) double {mustBePositive} = 6
    options.minHeightFactor (1,1) double {mustBePositive} = 4
    options.minWidthBins (1,1) double {mustBePositive} = 2
    options.minDistanceDa (1,1) double {mustBePositive} = 0.04
    options.minPeakHeight (1,1) double = 0
    options.minProminenceAbsolute (1,1) double {mustBeNonnegative} = 0
end

% Resolve mc vector
if istable(posOrMc)
    if ~ismember('mc', posOrMc.Properties.VariableNames)
        error('massSpecFindPeaks:missingMc', 'pos table must contain mc column.');
    end
    mc = posOrMc.mc;
    totalPulses = NaN;
    if options.usePulseNormalization && ismember('deltaP', posOrMc.Properties.VariableNames)
        totalPulses = sum(posOrMc.deltaP, 'omitnan');
    end
else
    mc = posOrMc;
    totalPulses = NaN;
end

mc = mc(:);
mc = mc(isfinite(mc));
if isempty(mc)
    peakTable = table();
    info = struct();
    return;
end

% Histogram
binWidth = options.binWidth;
if any(isnan(options.mcRange))
    mcMin = min(mc);
    mcMax = max(mc);
else
    mcMin = options.mcRange(1);
    mcMax = options.mcRange(2);
end
mcMax = max(mcMax, mcMin + binWidth);
maxEdge = ceil(mcMax / binWidth) * binWidth;
edges = mcMin:binWidth:maxEdge;
counts = histcounts(mc, edges);
centers = edges(1:end-1) + binWidth/2;

% Normalize to counts/Da[/pulse]
scaleLabel = 'counts/Da';
if options.usePulseNormalization && isfinite(totalPulses) && totalPulses > 0
    counts = counts / (totalPulses * binWidth);
    scaleLabel = 'counts/Da/pulse';
else
    counts = counts / binWidth;
end

% Baseline via moving quantile
baselineBins = max(3, round(options.baselineWindow / binWidth));
if mod(baselineBins, 2) == 0
    baselineBins = baselineBins + 1;
end
if exist('movquantile', 'file')
    baseline = movquantile(counts, options.baselineQuantile, baselineBins, 'omitnan');
elseif exist('movprctile', 'file')
    baseline = movprctile(counts, options.baselineQuantile * 100, baselineBins, 'omitnan');
else
    baseline = localMovingQuantile(counts, options.baselineQuantile, baselineBins);
end

% Residual and smoothing
residual = counts - baseline;
smoothBins = max(3, round(options.smoothSpan / binWidth));
if mod(smoothBins, 2) == 0
    smoothBins = smoothBins + 1;
end
smoothResidual = smoothdata(residual, 'sgolay', smoothBins);

% Noise estimate (robust)
noiseBins = max(3, round(options.noiseWindow / binWidth));
if mod(noiseBins, 2) == 0
    noiseBins = noiseBins + 1;
end
if exist('movmad', 'file')
    noise = movmad(residual, noiseBins, 1, 'omitnan') * 1.4826;
else
    med = movmedian(residual, noiseBins, 'omitnan');
    noise = movmedian(abs(residual - med), noiseBins, 'omitnan') * 1.4826;
end

% Peak detection
minDistBins = max(1, round(options.minDistanceDa / binWidth));
[pks, locs, widths, proms] = findpeaks( ...
    smoothResidual, ...
    'MinPeakDistance', minDistBins, ...
    'MinPeakWidth', options.minWidthBins, ...
    'MinPeakHeight', options.minPeakHeight);

if isempty(locs)
    peakTable = table();
else
    sigma = noise(locs);
    sigma(sigma == 0) = eps;
    noisePos = noise(noise > 0);
    if ~isempty(noisePos)
        sigmaFloor = quantile(noisePos, options.noiseFloorQuantile);
        sigma = max(sigma, sigmaFloor);
    end
    minProm = max(options.minProminenceAbsolute, options.minProminenceFactor * sigma);
    keep = proms >= minProm & pks >= options.minHeightFactor * sigma;

    locs = locs(keep);
    pks = pks(keep);
    widths = widths(keep);
    proms = proms(keep);
    [leftBases, rightBases] = estimatePeakBases(smoothResidual, locs);
    sigma = sigma(keep);

    % Peak area (baseline-subtracted)
    areas = zeros(size(locs));
    for i = 1:numel(locs)
        l = max(1, round(leftBases(i)));
        r = min(numel(residual), round(rightBases(i)));
        areas(i) = sum(max(residual(l:r), 0)) * binWidth;
    end

    peakMc = centers(locs);
    peakHeight = pks;
    peakProm = proms;
    widthBins = widths;
    widthDa = widths * binWidth;
    leftMc = centers(max(1, round(leftBases)));
    rightMc = centers(min(numel(centers), round(rightBases)));
    snr = peakProm ./ sigma;

    peakTable = table(peakMc(:), peakHeight(:), peakProm(:), snr(:), widthBins(:), widthDa(:), ...
        leftMc(:), rightMc(:), areas(:), ...
        'VariableNames', {'mc', 'height', 'prominence', 'snr', 'widthBins', 'widthDa', 'leftMc', 'rightMc', 'area'});
end

info = struct();
info.centers = centers;
info.edges = edges;
info.counts = counts;
info.baseline = baseline;
info.residual = residual;
info.smoothResidual = smoothResidual;
info.noise = noise;
info.binWidth = binWidth;
info.scaleLabel = scaleLabel;
info.totalPulses = totalPulses;
info.settings = options;
end

function [leftBases, rightBases] = estimatePeakBases(signal, locs)
    n = numel(signal);
    leftBases = zeros(size(locs));
    rightBases = zeros(size(locs));
    for i = 1:numel(locs)
        l = locs(i);
        r = locs(i);
        while l > 1 && signal(l) > 0
            l = l - 1;
        end
        while r < n && signal(r) > 0
            r = r + 1;
        end
        leftBases(i) = l;
        rightBases(i) = r;
    end
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
