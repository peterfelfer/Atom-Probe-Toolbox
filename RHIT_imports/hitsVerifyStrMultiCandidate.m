function report = hitsVerifyStrMultiCandidate(hitsFile, strFile, options)
% HITSVERIFYSTRMULTICANDIDATE Verify the leading HITS/STR multi-pair signal.
%
% report = hitsVerifyStrMultiCandidate(hitsFile, strFile)
%
% Focused verification for the candidate found by hitsSearchStrMultiPairs:
% HITS multi-hit continuation byte h2.b2 against STR detyt2 difference from
% the same event to the next STR event.
%
% OPTIONS:
%   eventOffset  - HITS eventIdx -> STR eventIdx offset (default: 0)
%   pairMode     - "eventPlus1" or "loadedPlus1" (default: "eventPlus1")
%   nShuffle     - number of shuffled controls (default: 200)
%   verbose      - print report (default: true)
%
% See also: hitsSearchStrMultiPairs
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.eventOffset (1,1) double = 0
    options.pairMode (1,1) string {mustBeMember(options.pairMode, ...
        ["eventPlus1","loadedPlus1"])} = "eventPlus1"
    options.nShuffle (1,1) double = 200
    options.verbose (1,1) logical = true
end

[x, y] = candidateVectors(hitsFile, strFile, options);
ok = isfinite(x) & isfinite(y);
x = double(x(ok));
y = double(y(ok));

fit = polyfit(x, y, 1);
yFit = polyval(fit, x);
resid = y - yFit;
rho = corr(x, y, 'Rows', 'complete');
rhoDetrended = corr(detrendLocal(x), detrendLocal(y), 'Rows', 'complete');
robust = robustMasks(x, y);

shuffleRho = nan(options.nShuffle, 1);
rng(1);
for k = 1:options.nShuffle
    shuffleRho(k) = corr(x(randperm(numel(x))), y, 'Rows', 'complete');
end

binTable = binnedTrend(x, y);

report = struct();
report.hitsFile = hitsFile;
report.strFile = strFile;
report.pairMode = char(options.pairMode);
report.n = numel(x);
report.predictor = 'h2_b2';
report.target = 'd_detyt2';
report.correlation = rho;
report.detrendedCorrelation = rhoDetrended;
report.slope = fit(1);
report.intercept = fit(2);
report.residualStd = std(resid);
report.shuffleAbsCorrelation95 = prctile(abs(shuffleRho), 95);
report.shuffleAbsCorrelationMax = max(abs(shuffleRho));
report.robust = robust;
report.binnedTrend = binTable;

if options.verbose
    fprintf('\n--- HITS/STR multi-hit candidate verification ---\n');
    fprintf('Pair mode: %s\n', options.pairMode);
    fprintf('Predictor: %s  Target: %s\n', report.predictor, report.target);
    fprintf('n=%d  r=%.4f  detrended r=%.4f\n', ...
        report.n, report.correlation, report.detrendedCorrelation);
    fprintf('Linear fit: d_detyt2 = %.3f*h2_b2 %+ .3f TDC\n', ...
        report.slope, report.intercept);
    fprintf('Residual std: %.1f TDC\n', report.residualStd);
    fprintf('Shuffle control |r| 95%%=%.4f max=%.4f (%d shuffles)\n', ...
        report.shuffleAbsCorrelation95, report.shuffleAbsCorrelationMax, ...
        options.nShuffle);
    fprintf('Robust subsets:\n');
    disp(struct2table(report.robust));
    fprintf('Binned trend:\n');
    disp(report.binnedTrend);
end
end

function [x, y] = candidateVectors(hitsFile, strFile, options)
[events, ~, raw] = hitsScanEvents(hitsFile);
eventFields = double(events.nFields);
completeBlocks = floor(max(0, eventFields - 1) / 2);
candidate = events.type == 'MULTI' & completeBlocks >= 2;
eventIdx = find(candidate);

n = numel(eventIdx);
x = nan(n, 1);
for k = 1:n
    ei = eventIdx(k);
    payloadOffset = double(events.byteOffset(ei)) + 4;
    x(k) = double(raw.bytes(payloadOffset + 8 + 3)); % hit2 byte 2
end

[strHits, ~] = strLoad(strFile);
strHits = strCalculatePositions(strHits);
rowByEvent = zeros(max(double(strHits.eventIdx)) + abs(options.eventOffset) + 3, ...
    1, 'uint32');
rowByEvent(double(strHits.eventIdx)) = uint32((1:height(strHits))');

baseEvent = double(eventIdx) + options.eventOffset;
rowA = zeros(size(baseEvent));
inBounds = baseEvent >= 1 & baseEvent <= numel(rowByEvent);
rowA(inBounds) = double(rowByEvent(baseEvent(inBounds)));

switch options.pairMode
    case "eventPlus1"
        rowB = zeros(size(baseEvent));
        eventB = baseEvent + 1;
        inBounds = eventB >= 1 & eventB <= numel(rowByEvent);
        rowB(inBounds) = double(rowByEvent(eventB(inBounds)));
    case "loadedPlus1"
        rowB = rowA + 1;
end

valid = rowA >= 1 & rowA <= height(strHits) & ...
        rowB >= 1 & rowB <= height(strHits);
y = nan(n, 1);
y(valid) = double(strHits.detyt2(rowB(valid))) - ...
    double(strHits.detyt2(rowA(valid)));
end

function T = binnedTrend(x, y)
[u, ~, ic] = unique(x);
count = accumarray(ic, 1);
meanY = accumarray(ic, y, [], @mean);
medianY = accumarray(ic, y, [], @median);
stdY = accumarray(ic, y, [], @std);
T = table(u, count, meanY, medianY, stdY, ...
    'VariableNames', {'h2_b2','count','mean_d_detyt2', ...
                      'median_d_detyt2','std_d_detyt2'});
end

function s = robustMasks(x, y)
normalByte = x <= 64;
yLimits = prctile(y, [1 99]);
centralY = y >= yLimits(1) & y <= yLimits(2);
madY = mad(y, 1);
if madY > 0
    madCentralY = abs(y - median(y)) <= 6 * madY;
else
    madCentralY = true(size(y));
end

s = struct();
s.all_n = numel(x);
s.x_le_64_n = nnz(normalByte);
s.x_le_64_r = corrSubset(x, y, normalByte);
s.x_le_64_detrended_r = corrSubset(detrendLocal(x), detrendLocal(y), normalByte);
s.central99_n = nnz(centralY);
s.central99_r = corrSubset(x, y, centralY);
s.central99_detrended_r = corrSubset(detrendLocal(x), detrendLocal(y), centralY);
s.x_le_64_central99_n = nnz(normalByte & centralY);
s.x_le_64_central99_r = corrSubset(x, y, normalByte & centralY);
s.x_le_64_central99_detrended_r = corrSubset(detrendLocal(x), ...
    detrendLocal(y), normalByte & centralY);
s.x_le_64_mad_n = nnz(normalByte & madCentralY);
s.x_le_64_mad_r = corrSubset(x, y, normalByte & madCentralY);
s.x_le_64_mad_detrended_r = corrSubset(detrendLocal(x), detrendLocal(y), ...
    normalByte & madCentralY);
end

function rho = corrSubset(x, y, sel)
ok = sel(:) & isfinite(x(:)) & isfinite(y(:));
if nnz(ok) < 10
    rho = NaN;
else
    rho = corr(x(ok), y(ok), 'Rows', 'complete');
end
end

function y = detrendLocal(x)
x = double(x(:));
w = min(101, max(11, 2*floor(numel(x) / 1000) + 1));
y = x - movmean(x, w, 'omitnan');
end
