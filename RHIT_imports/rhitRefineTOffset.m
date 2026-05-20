function [tOffsetRefined, info] = rhitRefineTOffset(hits, refSpectrum, bowl, tOffsetGuess, options)
% RHITREFINETOFFSET Refine t_offset by cross-correlation against a reference.
%
% [t, info] = rhitRefineTOffset(hits, refSpectrum, bowl, tOffsetGuess)
%
% Holds the bowl fixed, then searches t_offset (1-D) to maximise the
% cross-correlation of the computed mc histogram against `refSpectrum`
% (typically the EPOS-derived mc histogram).  Sample-agnostic: no peak
% windows, no peak assignments — works for any composition.
%
% INPUTS:
%   hits         - rhitLoad table (detx, dety mm; tof ns; VDC V)
%   refSpectrum  - struct with fields:
%                    .edges   bin edges (Da)
%                    .values  reference counts per bin
%                  (typically built from EPOS via rhitBuildRefSpectrum)
%   bowl         - struct accepted by rhitEvaluateBowl
%   tOffsetGuess - starting estimate (ns) — usually rhitTOffsetXCorr's value
%
% NAME-VALUE OPTIONS:
%   'sampleSize'  - subsample # ions for fast histogramming (default 1e6)
%   'coarseHalf'  - coarse search half-range in ns (default 20)
%   'coarseStep'  - coarse step ns (default 1)
%   'fineHalf'    - fine search half-range ns (default 2)
%   'fineStep'    - fine step ns (default 0.1)
%   'mcLo'/'mcHi' - mc range used for the cost (default 0.5/80)
%
% OUTPUT:
%   tOffsetRefined - the optimal t_offset (ns)
%   info           - struct with cost trace and best cost
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    refSpectrum struct
    bowl struct
    tOffsetGuess (1,1) double
    options.sampleSize (1,1) double = 1e6
    options.coarseHalf (1,1) double = 20
    options.coarseStep (1,1) double = 1
    options.fineHalf (1,1) double = 2
    options.fineStep (1,1) double = 0.1
    options.mcLo (1,1) double = 0.5
    options.mcHi (1,1) double = 80
    options.detrend (1,1) logical = true   % subtract sliding-percentile baseline
    options.detrendWin (1,1) double = 151  % bins; ~1.5 Da at 0.01 Da binning
end

% Subsample ions
nP = height(hits);
if nP > options.sampleSize
    rng(0);
    idx = randperm(nP, round(options.sampleSize));
else
    idx = (1:nP)';
end
sVDC = double(hits.VDC(idx));
sTof = double(hits.tof(idx));
sX   = double(hits.detx(idx));
sY   = double(hits.dety(idx));

bowlVDC = rhitEvaluateBowl(bowl, sX, sY) .* sVDC;

edges = double(refSpectrum.edges(:));
centres = 0.5 * (edges(1:end-1) + edges(2:end));
hRef = double(refSpectrum.values(:));
mask = centres >= options.mcLo & centres <= options.mcHi;

if options.detrend
    % Detrended cross-correlation: subtract a sliding 25th-percentile
    % baseline so the high pMass baseline (~10^3 counts) doesn't drown
    % out the peak signal.  Required when refSpectrum is the embedded
    % RHIT pMass; harmless when it's a clean EPOS histogram.
    hRefRef = detrendBase(hRef, options.detrendWin);
    nrm = sqrt(sum(hRefRef(mask) .^ 2));
    hRefNorm = hRefRef / max(nrm, 1e-12);
    useDot = true;
else
    hRefLog = log1p(hRef);
    hRefNorm = hRefLog / max(1e-9, sum(hRefLog));
    useDot = false;
end

cost = @(t) computeCost(t, bowlVDC, sTof, edges, hRefNorm, mask, options, useDot);

% Coarse → fine 1D search
bestT = tOffsetGuess;
bestC = cost(bestT);
costTrace = [];
coarseGrid = -options.coarseHalf:options.coarseStep:options.coarseHalf;
for dt = coarseGrid
    tt = tOffsetGuess + dt;
    c_v = cost(tt);
    costTrace(end+1, :) = [tt, c_v]; %#ok<AGROW>
    if c_v < bestC
        bestC = c_v;
        bestT = tt;
    end
end
fineGrid = -options.fineHalf:options.fineStep:options.fineHalf;
for dt = fineGrid
    tt = bestT + dt;
    c_v = cost(tt);
    costTrace(end+1, :) = [tt, c_v]; %#ok<AGROW>
    if c_v < bestC
        bestC = c_v;
        bestT = tt;
    end
end

tOffsetRefined = bestT;
info = struct();
info.bestCost = bestC;
info.costTrace = costTrace;
end


function c = computeCost(t, bowlVDC, sTof, edges, hRefNorm, mask, opt, useDot)
tofCorr = sTof - t;
mcT = bowlVDC .* (tofCorr .* tofCorr);
hT = histcounts(mcT(mcT > 0 & mcT < opt.mcHi + 1), edges);
hT = double(hT(:));
if useDot
    hTRef = detrendBase(hT, opt.detrendWin);
    nrm = sqrt(sum(hTRef(mask) .^ 2));
    hTNorm = hTRef / max(nrm, 1e-12);
else
    hTLog = log1p(hT);
    hTNorm = hTLog / max(1e-9, sum(hTLog));
end
% Maximise cross-correlation at zero lag → minimise its negative
c = -sum(hTNorm(mask) .* hRefNorm(mask));
end


function out = detrendBase(h, win)
%DETRENDBASE  Subtract a sliding-window 25th-percentile baseline.
% Each bin's baseline is approximated by a smoothed map of the local
% lower quartile, isolating the peaks above the continuous noise.
n = numel(h);
half = floor(win / 2);
% Coarse-grid baseline: lower quartile of each anchor's window
step = max(1, floor(win / 4));
anchors = 1:step:n;
baseAnchors = zeros(size(anchors));
for k = 1:numel(anchors)
    a = anchors(k);
    lo = max(1, a - half);
    hi = min(n, a + half);
    s = sort(h(lo:hi));
    baseAnchors(k) = s(max(1, ceil(0.25 * numel(s))));
end
% Interpolate baseline to every bin and smooth with a uniform kernel
base = interp1(anchors, baseAnchors, 1:n, 'linear', 'extrap')';
kernel = ones(win, 1) / win;
base = conv(base, kernel, 'same');
out = max(h - base, 0);
end
