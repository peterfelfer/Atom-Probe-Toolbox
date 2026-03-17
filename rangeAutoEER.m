function [rangeTable, info] = rangeAutoEER(posOrMc, peakMcPositions, options)
% RANGEAUTOEER Automatic ranging using Equal Error Rate criterion.
%
% [rangeTable, info] = rangeAutoEER(pos, peakMc)
% [rangeTable, info] = rangeAutoEER(pos, peakMc, 'purityFactor', 2.0)
%
% Given a set of peak positions, estimates the background and places range
% boundaries where missed signal atoms equal included background atoms
% (Equal Error Rate). A purityFactor parameter controls the trade-off
% between purity (high purityFactor) and recovery (low purityFactor).
%
% The boundary equation for each half of a peak is:
%   missed_signal = purityFactor * included_background
%
% INPUT
%   posOrMc          - pos table with mc column, or numeric mc vector (Da)
%   peakMcPositions  - numeric vector of peak m/c locations (Da)
%
% OPTIONS
%   'binWidth'       - histogram bin width in Da (default: 0.01)
%   'purityFactor'   - trade-off parameter k (default: 1.0)
%                      k = 1.0  balanced (EER, maximises classification probability)
%                      k > 1.0  conservative / high purity (e.g. 2.0)
%                      k < 1.0  inclusive / high recovery (e.g. 0.5)
%   'mcRange'        - [min max] mass-to-charge range in Da (default: data range)
%   'background'     - pre-computed background counts per bin (default: [])
%                      If empty, background is estimated automatically via
%                      A/sqrt(mc) fitted to the lower quantile of the spectrum.
%   'bgWindow'       - margin around peaks excluded from BG fit in Da (default: 0.3)
%   'minPurity'      - minimum purity to keep a range (default: 0.0)
%   'minCounts'      - minimum background-corrected counts to keep (default: 0)
%   'showPlot'       - display diagnostic plot (default: false)
%
% OUTPUT
%   rangeTable - table with columns:
%       mcbegin     - left range boundary (Da)
%       mcend       - right range boundary (Da)
%       peakMc      - peak maximum location (Da)
%       purity      - signal / (signal + background) in range
%       recovery    - signal in range / total peak signal
%       counts      - total counts in range
%       bgCounts    - estimated background counts in range
%       sigCounts   - background-corrected signal counts
%
%   info - struct with:
%       centers, counts, background  - spectrum arrays
%       peaks           - sorted peak positions used
%       binWidth        - bin width used
%       purityFactor    - purity factor used
%       bgCoefficient   - A from A/sqrt(mc) fit (if auto-estimated)
%       settings        - all options
%
% ALGORITHM
%   1. Build histogram of the m/c data.
%   2. Estimate background via A/sqrt(mc) fitted to the lower quantile of
%      the spectrum, or use user-provided background.
%   3. Partition the spectrum at valleys between adjacent peaks so each
%      peak's EER search is limited to its own region.
%   4. For each peak half (left and right of maximum), walk outward and
%      compute cumulative missed_signal and included_background.
%      Place boundary where missed = purityFactor * included_bg.
%   5. Filter ranges by purity and count thresholds.
%
% EXAMPLE
%   pos = posLoad('dataset.pos');
%
%   % Use with known peak positions (e.g. from ion identification)
%   peakMc = [14.0, 27.97, 28.47, 28.97];
%   [ranges, info] = rangeAutoEER(pos, peakMc, 'purityFactor', 1.0);
%   disp(ranges);
%
%   % Conservative ranging
%   [ranges, ~] = rangeAutoEER(pos, peakMc, 'purityFactor', 2.0);
%
%   % Inclusive ranging
%   [ranges, ~] = rangeAutoEER(pos, peakMc, 'purityFactor', 0.5);
%
%   % With custom background
%   [ranges, ~] = rangeAutoEER(pos, peakMc, 'background', myBg);
%
% SEE ALSO
%   massSpecFindPeaks, massSpecAutoAssignIons, backgroundEstimate,
%   rangeAdd, rangeAddAll
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    posOrMc
    peakMcPositions (:,1) double {mustBeNonempty}
    options.binWidth (1,1) double {mustBePositive} = 0.01
    options.purityFactor (1,1) double {mustBePositive} = 1.0
    options.mcRange (1,2) double = [NaN NaN]
    options.background double = []
    options.bgWindow (1,1) double {mustBePositive} = 0.3
    options.minPurity (1,1) double {mustBeNonnegative} = 0.0
    options.minCounts (1,1) double {mustBeNonnegative} = 0
    options.fitTails (1,1) logical = true
    options.showPlot (1,1) logical = false
end

%% Resolve mc vector
if istable(posOrMc)
    if ~ismember('mc', posOrMc.Properties.VariableNames)
        error('rangeAutoEER:missingMc', 'pos table must contain mc column.');
    end
    mc = posOrMc.mc;
else
    mc = posOrMc;
end
mc = mc(:);
mc = mc(isfinite(mc));
if isempty(mc)
    rangeTable = table();
    info = struct();
    return;
end

%% Build histogram
binWidth = options.binWidth;
if any(isnan(options.mcRange))
    mcMin = min(mc);
    mcMax = max(mc);
else
    mcMin = options.mcRange(1);
    mcMax = options.mcRange(2);
end
edges = mcMin:binWidth:(ceil(mcMax / binWidth) * binWidth + binWidth);
rawCounts = histcounts(mc, edges);
centers = edges(1:end-1) + binWidth / 2;
totalIons = numel(mc);

%% Background
bgCoefficient = NaN;
if ~isempty(options.background)
    % User-provided background
    bgRaw = options.background(:)';
    if numel(bgRaw) ~= numel(centers)
        error('rangeAutoEER:bgSizeMismatch', ...
            'Background length (%d) must match number of bins (%d).', ...
            numel(bgRaw), numel(centers));
    end
else
    % Auto-estimate: global A/sqrt(mc) fit to inter-peak regions.
    % This is used as the info output and as a fallback. The actual EER
    % computation uses per-peak local background (linear interpolation
    % between valley endpoints), which is more accurate for dense spectra.
    [bgRaw, bgCoefficient] = estimateBackgroundInvSqrt(centers, rawCounts, peakMcPositions, options.bgWindow);
end

%% Prepare peaks — sort and find bin indices
peakMcPositions = sort(peakMcPositions(:));

% Remove peaks outside the histogram range
peakMcPositions = peakMcPositions(peakMcPositions >= centers(1) & ...
                                  peakMcPositions <= centers(end));
nPeaks = numel(peakMcPositions);

if nPeaks == 0
    rangeTable = table();
    info = struct('centers', centers, 'counts', rawCounts, ...
        'background', bgRaw, 'peaks', peakMcPositions, ...
        'binWidth', binWidth, 'settings', options);
    return;
end

% For each peak, find the actual local maximum in the raw spectrum near
% the given position. Use a small search radius (±0.05 Da) to avoid
% snapping to a neighboring peak.
pkBins = zeros(nPeaks, 1);
snapRadius = round(0.05 / binWidth);
for p = 1:nPeaks
    nominalBin = findBin(centers, peakMcPositions(p));
    lo = max(1, nominalBin - snapRadius);
    hi = min(numel(centers), nominalBin + snapRadius);
    [~, relIdx] = max(rawCounts(lo:hi));
    pkBins(p) = lo + relIdx - 1;
end

%% Partition spectrum at valleys between adjacent peaks
% Use the raw spectrum (not smoothed) to find the minimum between each
% pair of adjacent peaks. Savitzky-Golay smoothing causes ringing near
% sharp peaks, which misplaces valleys. The raw spectrum minimum is
% robust — noise only shifts it by ±1 bin, which is negligible.
leftLimit = ones(nPeaks, 1);
rightLimit = ones(nPeaks, 1) * numel(centers);

for p = 1:(nPeaks - 1)
    regionBins = pkBins(p):pkBins(p+1);
    if numel(regionBins) > 2
        [~, minIdx] = min(rawCounts(regionBins));
        valleyBin = regionBins(minIdx);
    else
        valleyBin = round((pkBins(p) + pkBins(p+1)) / 2);
    end
    rightLimit(p) = valleyBin;
    leftLimit(p+1) = valleyBin;
end

%% Fit power law tails for each peak (optional)
% The thermal tail of an APT peak follows a power law: signal ~ a * dx^(-c)
% where dx is the distance from the peak maximum. When enabled, the fitted
% tails of neighbouring peaks are added to the effective background, so the
% EER boundary accounts for the classification probability between adjacent
% peaks. When disabled (fitTails=false), only the monotone decrease
% constraint and valley partitioning are used (global background mode).
tailFits = struct('aRight', {}, 'cRight', {}, 'aLeft', {}, 'cLeft', {});
for p = 1:nPeaks
    tailFits(p).aRight = NaN; tailFits(p).cRight = NaN;
    tailFits(p).aLeft = NaN;  tailFits(p).cLeft = NaN;
end

if options.fitTails
for p = 1:nPeaks

    pkMc = centers(pkBins(p));

    % Right tail: fit from peak+0.05 to midpoint toward right neighbor
    if p < nPeaks
        midRight = (pkMc + centers(pkBins(p+1))) / 2;
    else
        midRight = centers(rightLimit(p));
    end
    fitMask = centers > pkMc + 0.05 & centers < midRight;
    if sum(fitMask) >= 5
        dx = centers(fitMask) - pkMc;
        sig = max(1, rawCounts(fitMask) - bgRaw(fitMask));
        pf = polyfit(log(dx), log(sig), 1);
        tailFits(p).cRight = -pf(1);
        tailFits(p).aRight = exp(pf(2));
    end

    % Left tail: fit from peak-0.05 to midpoint toward left neighbor
    if p > 1
        midLeft = (centers(pkBins(p-1)) + pkMc) / 2;
    else
        midLeft = centers(leftLimit(p));
    end
    fitMask = centers < pkMc - 0.05 & centers > midLeft;
    if sum(fitMask) >= 5
        dx = pkMc - centers(fitMask);
        sig = max(1, rawCounts(fitMask) - bgRaw(fitMask));
        pf = polyfit(log(dx), log(sig), 1);
        tailFits(p).cLeft = -pf(1);
        tailFits(p).aLeft = exp(pf(2));
    end
end
end % if options.fitTails

%% Compute EER boundaries with neighbor tail subtraction
k = options.purityFactor;
mcbegin = nan(nPeaks, 1);
mcend = nan(nPeaks, 1);

for p = 1:nPeaks
    lBin = leftLimit(p);
    rBin = rightLimit(p);
    regionBins = lBin:rBin;
    nRegion = numel(regionBins);
    pkInRegion = pkBins(p) - lBin + 1;
    pkInRegion = max(1, min(nRegion, pkInRegion));

    regionCts = rawCounts(regionBins);
    regionBg = bgRaw(regionBins);
    regionMc = centers(regionBins);

    % Subtract the right neighbor's left tail from the right half
    neighborTailRight = zeros(1, nRegion);
    if p < nPeaks && isfinite(tailFits(p+1).aLeft) && tailFits(p+1).cLeft > 0
        nbrMc = centers(pkBins(p+1));
        dx = nbrMc - regionMc;
        valid = dx > 0.02;
        neighborTailRight(valid) = tailFits(p+1).aLeft .* dx(valid).^(-tailFits(p+1).cLeft);
    end

    % Subtract the left neighbor's right tail from the left half
    neighborTailLeft = zeros(1, nRegion);
    if p > 1 && isfinite(tailFits(p-1).aRight) && tailFits(p-1).cRight > 0
        nbrMc = centers(pkBins(p-1));
        dx = regionMc - nbrMc;
        valid = dx > 0.02;
        neighborTailLeft(valid) = tailFits(p-1).aRight .* dx(valid).^(-tailFits(p-1).cRight);
    end

    % Effective background = global bg + neighbor tail contributions
    effectiveBg = regionBg + neighborTailRight + neighborTailLeft;

    % EER on left half
    mcbegin(p) = findEERBoundaryLocal( ...
        regionCts, effectiveBg, regionMc, pkInRegion, k, 'left');

    % EER on right half
    mcend(p) = findEERBoundaryLocal( ...
        regionCts, effectiveBg, regionMc, pkInRegion, k, 'right');
end

%% Compute quality metrics and cross-contamination matrix
purity = nan(nPeaks, 1);
recovery = nan(nPeaks, 1);
counts = nan(nPeaks, 1);
bgCounts = nan(nPeaks, 1);
sigCounts = nan(nPeaks, 1);

% Cross-contamination matrix: contamination(p, q) = estimated number of
% atoms from peak q's tail that fall inside range p. This can be used to
% correct the composition: corrected_counts(p) = counts(p) - sum(contamination(p,:))
contamination = zeros(nPeaks, nPeaks);

for p = 1:nPeaks
    inRange = centers >= mcbegin(p) & centers <= mcend(p);
    counts(p) = sum(rawCounts(inRange));
    bgCounts(p) = sum(bgRaw(inRange));

    % Compute contamination from each neighbour's tail
    if options.fitTails
        for q = 1:nPeaks
            if q == p; continue; end
            tf = tailFits(q);
            [~, qBin] = min(abs(centers - peakMcPositions(q)));
            qMc = centers(qBin);

            if qMc < centers(pkBins(p)) && isfinite(tf.aRight) && tf.cRight > 0
                % Peak q is to the left; its right tail extends into range p
                dx = centers(inRange) - qMc;
                valid = dx > 0.02;
                if any(valid)
                    contamination(p, q) = sum(tf.aRight .* dx(valid).^(-tf.cRight));
                end
            elseif qMc > centers(pkBins(p)) && isfinite(tf.aLeft) && tf.cLeft > 0
                % Peak q is to the right; its left tail extends into range p
                dx = qMc - centers(inRange);
                valid = dx > 0.02;
                if any(valid)
                    contamination(p, q) = sum(tf.aLeft .* dx(valid).^(-tf.cLeft));
                end
            end
        end
    end

    totalContam = sum(contamination(p, :));
    sigCounts(p) = max(0, counts(p) - bgCounts(p) - totalContam);

    if counts(p) > 0
        purity(p) = sigCounts(p) / counts(p);
    else
        purity(p) = 0;
    end

    % Recovery: fraction of total signal in this peak's search region
    inRegion = centers >= centers(leftLimit(p)) & centers <= centers(rightLimit(p));
    totalSig = max(1, sum(rawCounts(inRegion)) - sum(bgRaw(inRegion)));
    recovery(p) = min(1, sigCounts(p) / totalSig);
end

%% Filter ranges
keepRange = purity >= options.minPurity & sigCounts >= options.minCounts;
keepRange = keepRange & isfinite(mcbegin) & isfinite(mcend) & mcend > mcbegin;

%% Build output table
rangeTable = table(mcbegin, mcend, peakMcPositions, purity, recovery, ...
    counts, bgCounts, sigCounts, ...
    'VariableNames', {'mcbegin', 'mcend', 'peakMc', 'purity', 'recovery', ...
    'counts', 'bgCounts', 'sigCounts'});
rangeTable = rangeTable(keepRange, :);
rangeTable = sortrows(rangeTable, 'mcbegin');

%% Build info struct
info = struct();
info.centers = centers;
info.counts = rawCounts;
info.background = bgRaw;
info.bgNorm = bgRaw / (binWidth * totalIons);
info.peaks = peakMcPositions;
info.binWidth = binWidth;
info.purityFactor = k;
info.totalIons = totalIons;
info.bgCoefficient = bgCoefficient;
info.tailFits = tailFits;
info.contamination = contamination;
info.settings = options;

%% Diagnostic plot
if options.showPlot
    plotDiagnostic(centers, rawCounts, bgRaw, rangeTable, peakMcPositions, k);
end

end


%% ========== Core EER boundary finder ==========

function mcBoundary = findEERBoundaryLocal(counts, bg, centers, pkIdx, k, side)
% Find range boundary using the EER criterion within a local region.
%
% counts, bg, centers are arrays for the region around one peak.
% pkIdx is the index of the peak maximum within this region.

nBins = numel(counts);

switch side
    case 'left'
        if pkIdx <= 1
            mcBoundary = centers(1);
            return;
        end
        idx = pkIdx:-1:1;

    case 'right'
        if pkIdx >= nBins
            mcBoundary = centers(end);
            return;
        end
        idx = pkIdx:nBins;
end

cts = counts(idx);
bgCts = bg(idx);
mc = centers(idx);

% Signal above background at each bin
signal = max(0, cts - bgCts);

% Enforce monotone decrease: the signal from THIS peak can only decrease
% as we walk away from its maximum. Any increase is either noise or a
% neighboring peak's contribution — cap it at the previous level.
% This is robust against noise (no smoothing/thresholds needed) and
% naturally prevents the range from extending into a neighbor's tail.
monoSignal = signal;
for i = 2:numel(monoSignal)
    monoSignal(i) = min(monoSignal(i-1), monoSignal(i));
end

totalSignal = sum(monoSignal);

if totalSignal == 0
    mcBoundary = mc(1);
    return;
end

% Cumulative signal and background as we walk outward from peak
cumSignal = cumsum(monoSignal);
cumBg = cumsum(bgCts);

% Missed signal = what we'd lose if range boundary is at position i
missedSignal = totalSignal - cumSignal;

% Included background = BG atoms inside the range
includedBg = cumBg;

% Find crossing: missed_signal = k * included_bg
balance = missedSignal - k * includedBg;

crossIdx = find(balance <= 0, 1, 'first');

if isempty(crossIdx)
    % EER not reached — stop at the last bin with nonzero mono signal
    lastNonzero = find(monoSignal > 0, 1, 'last');
    if isempty(lastNonzero)
        mcBoundary = mc(1);
    else
        mcBoundary = mc(lastNonzero);
    end
else
    % Interpolate for sub-bin precision
    if crossIdx > 1 && balance(crossIdx - 1) > 0
        frac = balance(crossIdx - 1) / (balance(crossIdx - 1) - balance(crossIdx));
        mcBoundary = mc(crossIdx - 1) + frac * (mc(crossIdx) - mc(crossIdx - 1));
    else
        mcBoundary = mc(crossIdx);
    end
end

end


%% ========== Helper functions ==========

function idx = findBin(centers, mcVal)
    [~, idx] = min(abs(centers - mcVal));
end


function [bgRaw, A] = estimateBackgroundInvSqrt(centers, rawCounts, peakMc, marginDa)
% Estimate background as A/sqrt(mc) by fitting to inter-peak regions.
%
% Masks ±marginDa around each peak, fits A/sqrt(mc) to the remaining
% (unmasked) bins, and evaluates the fit at all bin centers.
%
% This is the physics-based background model: a constant background in
% TOF space transforms to A/sqrt(mc) in mass spectrum space via the
% mc ~ t^2 relationship.

    % Mask bins near peaks
    bgMask = true(size(centers));
    bgMask(centers <= 1) = false;  % avoid 1/sqrt divergence

    for i = 1:numel(peakMc)
        bgMask = bgMask & ~(centers >= peakMc(i) - marginDa & ...
                            centers <= peakMc(i) + marginDa);
    end

    xBg = centers(bgMask);
    yBg = rawCounts(bgMask);

    if numel(xBg) < 10
        % Not enough unmasked bins — fall back to moving lower quantile
        windowBins = max(3, round(marginDa * 2 / mean(diff(centers))));
        if mod(windowBins, 2) == 0
            windowBins = windowBins + 1;
        end
        n = numel(rawCounts);
        half = floor(windowBins / 2);
        lowerEnv = zeros(size(rawCounts));
        for i = 1:n
            i1 = max(1, i - half);
            i2 = min(n, i + half);
            lowerEnv(i) = prctile(rawCounts(i1:i2), 10);
        end
        fitMask = centers > 1 & lowerEnv > 0;
        if sum(fitMask) < 3
            bgRaw = lowerEnv;
            A = NaN;
            return;
        end
        xBg = centers(fitMask);
        yBg = lowerEnv(fitMask);
    end

    % The unmasked bins may still contain peak tails, which inflate the
    % background estimate. Use only the lower quartile of unmasked bins
    % (grouped by mc region) to fit through the true valley floor.
    % Transform to y*sqrt(x) space where background is constant = A.
    yTransformed = yBg(:) .* sqrt(xBg(:));
    % Keep only bins below the median in transformed space —
    % this excludes tail-contaminated bins and fits through the
    % true valley floor between peaks.
    q50 = quantile(yTransformed, 0.50);
    lowMask = yTransformed <= q50;
    if sum(lowMask) < 3
        lowMask = true(size(yTransformed));  % fallback: use all
    end
    xLow = xBg(lowMask);
    yLow = yBg(lowMask);

    % Least squares: y = A * x^(-0.5)
    invSqrtX = 1 ./ sqrt(xLow(:));
    A = sum(yLow(:) .* invSqrtX) / sum(invSqrtX.^2);
    A = max(0, A);

    bgRaw = A ./ sqrt(max(centers, 0.01));
end


function plotDiagnostic(centers, counts, bg, rangeTable, peakMcPositions, k)
% Plot diagnostic figure showing spectrum, background, peaks, and ranges

    figure('Name', sprintf('rangeAutoEER (k = %.2f)', k), 'Color', 'w');

    % Spectrum
    semilogy(centers, max(counts, 0.1), 'k-', 'LineWidth', 0.5, 'DisplayName', 'Spectrum');
    hold on;

    % Background
    semilogy(centers, max(bg, 0.01), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Background');

    % Peak markers
    for i = 1:numel(peakMcPositions)
        xline(peakMcPositions(i), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    end

    % Range shading
    yLim = ylim;
    for i = 1:height(rangeTable)
        x1 = rangeTable.mcbegin(i);
        x2 = rangeTable.mcend(i);
        patch([x1 x2 x2 x1], [yLim(1) yLim(1) yLim(2) yLim(2)], ...
            [0.3 0.6 1], 'FaceAlpha', 0.2, 'EdgeColor', [0.2 0.4 0.8], ...
            'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Range boundary lines
    if height(rangeTable) > 0
        xline(rangeTable.mcbegin(1), '-', 'Color', [0.2 0.4 0.8], ...
            'LineWidth', 1, 'DisplayName', 'Range boundaries');
        for i = 2:height(rangeTable)
            xline(rangeTable.mcbegin(i), '-', 'Color', [0.2 0.4 0.8], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
        end
        for i = 1:height(rangeTable)
            xline(rangeTable.mcend(i), '-', 'Color', [0.2 0.4 0.8], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end

    xlabel('mass-to-charge-state ratio [Da]');
    ylabel('counts per bin');
    title(sprintf('Automatic EER Ranging (k = %.2f, %d ranges)', ...
        k, height(rangeTable)));
    legend('Location', 'northeast');
    set(gca, 'YScale', 'log');
    hold off;
end
