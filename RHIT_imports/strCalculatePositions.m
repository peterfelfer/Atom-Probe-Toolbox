function hits = strCalculatePositions(hits)
% STRCALCPOSITIONS Compute detector positions and TOF from DLD timing data.
%
% hits = strCalculatePositions(hits)
%
% Computes detxRaw, detyRaw, detwRaw (delay-line positions in TDC counts)
% and tof (time-of-flight in TDC counts) from the six raw delay-line end
% timings (detxt1/t2, detyt1/t2, detwt1/t2).
%
% For partial hits (some channels missing), the 45-degree geometric
% constraint is used to recover the missing position:
%   w ≈ a*x + b*y + c  (fitted from complete events)
%
% DETECTOR GEOMETRY:
%   DL-x at 0°   → detxRaw = detxt1 - detxt2
%   DL-y at 90°  → detyRaw = detyt1 - detyt2
%   DL-w at 45°  → detwRaw = detwt1 - detwt2 (redundant diagonal)
%
% HIT CLASSIFICATION (added as 'hitType' column):
%   3 = full hit    (all 3 delay lines)
%   2 = partial hit (2 delay lines, 3rd recovered from geometry)
%   0 = incomplete  (fewer than 2 delay lines — no position possible)
%
% TOF CALCULATION:
%   tof = mean of available delay-line sums: (t1+t2)/2 per DL,
%   corrected for per-DL offset (fitted from complete events).
%
% INPUTS/OUTPUTS:
%   hits - Table from strLoad. Columns detxRaw, detyRaw, detwRaw, tof,
%          and hitType are added.
%
% EXAMPLE:
%   [hits, meta] = strLoad('data.STR');
%   hits = strCalculatePositions(hits);
%   histogram(hits.tof(hits.hitType >= 2), 1000);
%
% See also: strLoad
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
end

%% Identify complete delay lines
hasX = ~isnan(hits.detxt1) & ~isnan(hits.detxt2);
hasY = ~isnan(hits.detyt1) & ~isnan(hits.detyt2);
hasW = ~isnan(hits.detwt1) & ~isnan(hits.detwt2);

nFull = sum(hasX & hasY & hasW);
nPartial2 = sum((hasX & hasY & ~hasW) | (hasX & ~hasY & hasW) | (~hasX & hasY & hasW));
nIncomplete = height(hits) - nFull - nPartial2;

fprintf('Events: %d full (3 DL), %d partial (2 DL), %d incomplete (<2 DL)\n', ...
    nFull, nPartial2, nIncomplete);

%% Compute raw positions for complete delay lines
detxRaw = NaN(height(hits), 1);
detyRaw = NaN(height(hits), 1);
detwRaw = NaN(height(hits), 1);

detxRaw(hasX) = hits.detxt1(hasX) - hits.detxt2(hasX);
detyRaw(hasY) = hits.detyt1(hasY) - hits.detyt2(hasY);
detwRaw(hasW) = hits.detwt1(hasW) - hits.detwt2(hasW);

%% Fit geometric constraint from full hits: detwRaw = a*detxRaw + b*detyRaw + c
% Use sum-consistency filter to exclude mismatched events from the fit
full = hasX & hasY & hasW;
sumX_f = hits.detxt1(full) + hits.detxt2(full);
sumY_f = hits.detyt1(full) + hits.detyt2(full);
sumW_f = hits.detwt1(full) + hits.detwt2(full);
dXY = sumX_f - sumY_f;
dXW = sumX_f - sumW_f;
medXY = median(dXY); medXW = median(dXW);
sumGood = abs(dXY - medXY) < 50 & abs(dXW - medXW) < 50;

fitPool = find(full);
fitPool = fitPool(sumGood);
nFit = min(500000, numel(fitPool));
fitIdx = fitPool(1:nFit);

A = [detxRaw(fitIdx), detyRaw(fitIdx), ones(nFit, 1)];
p = A \ detwRaw(fitIdx);
a = p(1); b = p(2); c = p(3);

resid = detwRaw(fitIdx) - A * p;
fprintf('Geometry fit: w = %.4f*x + %.4f*y + %.1f (residual std=%.1f TDC)\n', ...
    a, b, c, std(resid));

%% Recover missing positions for partial hits (2 of 3 DLs present)

% Missing x: recover from y and w
missX = ~hasX & hasY & hasW;
if any(missX)
    detxRaw(missX) = (detwRaw(missX) - b * detyRaw(missX) - c) / a;
end

% Missing y: recover from x and w
missY = hasX & ~hasY & hasW;
if any(missY)
    detyRaw(missY) = (detwRaw(missY) - a * detxRaw(missY) - c) / b;
end

% Missing w: recover from x and y (not strictly needed, but complete the set)
missW = hasX & hasY & ~hasW;
if any(missW)
    detwRaw(missW) = a * detxRaw(missW) + b * detyRaw(missW) + c;
end

nRecovered = sum(missX) + sum(missY) + sum(missW);
fprintf('Recovered %d partial hits (%d missX, %d missY, %d missW)\n', ...
    nRecovered, sum(missX), sum(missY), sum(missW));

%% Compute TOF from delay-line sums
% Each DL sum = t1 + t2 = 2*TOF + DL_offset
% Fit DL offsets from full hits
sumX = hits.detxt1 + hits.detxt2;
sumY = hits.detyt1 + hits.detyt2;
sumW = hits.detwt1 + hits.detwt2;

% DL offsets relative to mean sum (from full hits)
meanSum = (sumX(fitIdx) + sumY(fitIdx) + sumW(fitIdx)) / 3;
offX = median(sumX(fitIdx) - meanSum);
offY = median(sumY(fitIdx) - meanSum);
offW = median(sumW(fitIdx) - meanSum);

fprintf('DL sum offsets: x=%+.1f, y=%+.1f, w=%+.1f\n', offX, offY, offW);

% TOF = average of corrected sums / 2
tof = NaN(height(hits), 1);
nDL = zeros(height(hits), 1);
tofAccum = zeros(height(hits), 1);

if any(hasX)
    corrX = (sumX - offX) / 2;
    tofAccum(hasX) = tofAccum(hasX) + corrX(hasX);
    nDL(hasX) = nDL(hasX) + 1;
end
if any(hasY)
    corrY = (sumY - offY) / 2;
    tofAccum(hasY) = tofAccum(hasY) + corrY(hasY);
    nDL(hasY) = nDL(hasY) + 1;
end
if any(hasW)
    corrW = (sumW - offW) / 2;
    tofAccum(hasW) = tofAccum(hasW) + corrW(hasW);
    nDL(hasW) = nDL(hasW) + 1;
end

hasAnyTof = nDL > 0;
tof(hasAnyTof) = tofAccum(hasAnyTof) ./ nDL(hasAnyTof);

%% Hit type classification: count of complete delay lines (0-3)
hitType = uint8(hasX) + uint8(hasY) + uint8(hasW);

%% Add columns to table
hits.detxRaw = detxRaw;
hits.detyRaw = detyRaw;
hits.detwRaw = detwRaw;
hits.tof     = tof;
hits.hitType = hitType;

fprintf('Result: %d events with position, %d with TOF\n', ...
    sum(~isnan(detxRaw) & ~isnan(detyRaw)), sum(~isnan(tof)));

end
