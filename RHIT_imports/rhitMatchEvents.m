function [matchedE, matchedR, info] = rhitMatchEvents(hits, epos, options)
% RHITMATCHEVENTS Drift-tracking RHIT-to-EPOS event matcher.
%
% [matchedE, matchedR] = rhitMatchEvents(hits, epos)
% [matchedE, matchedR, info] = rhitMatchEvents(hits, epos, 'tOffsetGuess', 291.6, ...)
%
% Two-pass matcher:
%   Pass 1 (drift model). Probe ~200 EPOS event indices spread across the
%       run.  At each, do a wide-window position+TOF search to find the
%       best RHIT match.  The resulting (ei, r_offset) pairs are filtered
%       (median-rate outlier rejection) and form a piecewise-linear
%       offset(ei) interpolator.
%   Pass 2 (tight match).  For every probe (every `stride`-th EPOS event)
%       use the drift model to predict the RHIT index, then accept the
%       single closest candidate within `[scan_back, scan_fwd]` that
%       passes (VDC, det_x, det_y, TOF) tolerances.
%
% INPUTS:
%   hits  - table from rhitLoad with detx, dety (mm), VDC (V), tof (ns).
%   epos  - table from posLoad with detx, dety (mm), VDC (V), tof (ns).
%
% NAME-VALUE OPTIONS:
%   'tOffsetGuess'  - initial TOF offset (ns).  If omitted, computed by
%                     rhitTOffsetXCorr.
%   'icf'           - epos.detx / hits.detx ratio (mm/mm).  Auto if [].
%   'stride'        - thin EPOS by this factor (default 200).
%   'scanBack'      - back-search width in events (default 50).
%   'scanFwd'       - forward-search width (default 2000).
%   'fallbackFwd'   - wider one-shot fallback (default 50000).
%   'vdcTol'        - V (default 0.5).
%   'posTol'        - mm (default 0.15).
%   'tofTol'        - ns around tOffsetGuess (default 1.5).
%   'driftAnchors'  - target anchor count (default 200).
%   'verbose'       - logical (default false).
%
% OUTPUTS:
%   matchedE / matchedR  - 1-based row indices into epos / hits.
%   info  - struct with anchor counts, failure counts, ICF, t_offset_used.
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    epos table
    options.tOffsetGuess (1,1) double = NaN
    options.icf (1,1) double = NaN
    options.stride (1,1) double = 200
    options.scanBack (1,1) double = 50
    options.scanFwd (1,1) double = 2000
    options.fallbackFwd (1,1) double = 50000
    options.vdcTol (1,1) double = 0.5
    options.posTol (1,1) double = 0.15
    options.tofTol (1,1) double = 1.5
    options.driftAnchors (1,1) double = 200
    % Cameca-style ROI selectors. Empty means "use all".
    options.vdcMin (1,1) double = -inf       % drop low-V (cap) region
    options.vdcMax (1,1) double = inf        % drop high-V (microfracture) region
    options.detRadiusMaxMm (1,1) double = inf  % XY FOV: keep |det| < r
    options.verbose (1,1) logical = false
end

% --- Pull arrays once (fast) ---
hVDC  = double(hits.VDC);     eVDC  = double(epos.VDC);
hDetx = double(hits.detx);    eDetx = double(epos.detx);
hDety = double(hits.dety);    eDety = double(epos.dety);
hTof  = double(hits.tof);     eTof  = double(epos.tof);

% --- Cameca-style ROI: zero out tolerance for events outside Z range
% or beyond the FOV. We don't actually drop them (would break index
% alignment) — we just refuse to match them, which is equivalent.
if isfinite(options.vdcMin) || isfinite(options.vdcMax) ...
        || isfinite(options.detRadiusMaxMm)
    keepE = eVDC >= options.vdcMin & eVDC <= options.vdcMax ...
        & sqrt(eDetx.^2 + eDety.^2) < options.detRadiusMaxMm;
    if options.verbose
        fprintf('rhitMatchEvents: ROI keeps %d / %d EPOS events (%.1f%%)\n', ...
            sum(keepE), numel(keepE), 100 * sum(keepE) / numel(keepE));
    end
else
    keepE = true(size(eVDC));
end
nR = numel(hVDC); nE = numel(eVDC);

% --- ICF (auto from a probe slice if not supplied) ---
icf = options.icf;
if isnan(icf)
    nP = min(1000, min(nR, nE));
    rxslice = hDetx(1:nP); exslice = eDetx(1:nP);
    valid = abs(rxslice) > 1.0 & abs(exslice ./ rxslice) < 2.0;
    if any(valid)
        icf = median(exslice(valid) ./ rxslice(valid));
    else
        icf = 1.0;
    end
end
if options.verbose
    fprintf('rhitMatchEvents: ICF = %.4f\n', icf);
end

% --- t_offset guess ---
tGuess = options.tOffsetGuess;
if isnan(tGuess)
    tGuess = rhitTOffsetXCorr(hTof, eTof);
    if options.verbose
        fprintf('rhitMatchEvents: t_offset (xcorr) = %.2f ns\n', tGuess);
    end
end

% --- Pass 1: drift anchors ---
if options.verbose
    fprintf('rhitMatchEvents: building drift model...\n');
end
[anchorEi, anchorOff] = buildDriftModel( ...
    hVDC, hDetx, hDety, hTof, eVDC, eDetx, eDety, eTof, ...
    icf, tGuess, options);

if numel(anchorEi) < 2
    warning('rhitMatchEvents:noAnchors', ...
        'Drift-anchor pass found <2 anchors; using r_offset=0 fallback.');
    anchorEi = [1; nE];
    anchorOff = [0; 0];
end

% --- Pass 2: tight matching with predicted offsets ---
posTol2 = options.posTol ^ 2;
maxMatches = floor(nE / options.stride) + 1;
matchedE = zeros(maxMatches, 1, 'int64');
matchedR = zeros(maxMatches, 1, 'int64');
nMatched = 0;
nFailed = 0;

probes = 1:options.stride:nE;
for pIdx = 1:numel(probes)
    ei = probes(pIdx);
    if ~keepE(ei)
        continue
    end
    offset = round(interp1(double(anchorEi), double(anchorOff), ...
        double(ei), 'linear', 'extrap'));
    riPredict = ei + offset;
    rLo = max(1, riPredict - options.scanBack);
    rHi = min(nR, riPredict + options.scanFwd);
    if rHi <= rLo
        continue
    end

    rIdx = rLo:rHi;
    dvdc = abs(hVDC(rIdx) - eVDC(ei));
    ddx = hDetx(rIdx) * icf - eDetx(ei);
    ddy = hDety(rIdx) * icf - eDety(ei);
    d2 = ddx .* ddx + ddy .* ddy;
    dtof = abs(hTof(rIdx) - eTof(ei) - tGuess);
    ok = dvdc < options.vdcTol & d2 < posTol2 & dtof < options.tofTol;
    if ~any(ok)
        nFailed = nFailed + 1;
        continue
    end
    score = d2 + 0.05 * (dvdc .^ 2) + 0.01 * (dtof .^ 2);
    score(~ok) = inf;
    [~, j] = min(score);
    rj = rLo + j - 1;
    nMatched = nMatched + 1;
    matchedE(nMatched) = ei;
    matchedR(nMatched) = rj;
end

matchedE = double(matchedE(1:nMatched));
matchedR = double(matchedR(1:nMatched));

info = struct();
info.icf = icf;
info.tOffsetGuess = tGuess;
info.nAnchors = numel(anchorEi);
info.nMatched = nMatched;
info.nFailed = nFailed;

if options.verbose
    fprintf('rhitMatchEvents: matched=%d / %d probed (failed=%d)\n', ...
        nMatched, numel(probes), nFailed);
end
end


% =========================================================================
% Internal helpers
% =========================================================================

function [anchorEi, anchorOff] = buildDriftModel(hVDC, hDetx, hDety, hTof, ...
    eVDC, eDetx, eDety, eTof, icf, tGuess, opt) %#ok<INUSL>

nE = numel(eVDC); nR = numel(hVDC);

% Anchor schedule: dense early, then linear-spread to nE
nEarly = 20;
nLate = max(20, opt.driftAnchors);
maxEarly = max(2, floor(nE / 100));
early = unique(round(logspace(0, log10(maxEarly), nEarly)));
late  = round(linspace(maxEarly, nE, nLate));
anchorsEi = unique([1, early(:)', late(:)']);
anchorsEi = anchorsEi(anchorsEi >= 1 & anchorsEi <= nE);

posTol2 = opt.posTol ^ 2;
anchorEi  = [];
anchorOff = [];
cur = 0;

for k = 1:numel(anchorsEi)
    ei = anchorsEi(k);
    maxDrift = max(2000, floor(double(ei) * 0.01));
    lo = max(1, ei + cur - 200);
    hi = min(nR, ei + cur + maxDrift);
    rj = anchorSearch(ei, lo, hi, hVDC, hDetx, hDety, hTof, ...
        eVDC, eDetx, eDety, eTof, icf, tGuess, posTol2, opt);
    if isempty(rj)
        % Wide fallback
        lo2 = max(1, ei + cur - 5000);
        hi2 = min(nR, ei + cur + 500000);
        rj = anchorSearch(ei, lo2, hi2, hVDC, hDetx, hDety, hTof, ...
            eVDC, eDetx, eDety, eTof, icf, tGuess, posTol2, opt);
    end
    if ~isempty(rj)
        anchorEi(end+1, 1)  = ei;          %#ok<AGROW>
        anchorOff(end+1, 1) = rj - ei;     %#ok<AGROW>
        cur = rj - ei;
    end
end

if numel(anchorEi) < 2
    return
end

% Outlier rejection: drift accumulates ~linearly with ei.  Drop anchors
% whose normalised drift rate deviates >5x from median.
rate = double(anchorOff) ./ max(double(anchorEi), 1);
if any(anchorEi > 1000)
    anchorRate = median(rate(anchorEi > 1000));
else
    anchorRate = 0;
end
if isnan(anchorRate); anchorRate = 0; end
plausibleMax = max(50000, 5 * abs(anchorRate) * nR);
keep = abs(double(anchorOff)) < plausibleMax;
if anchorRate ~= 0
    keep = keep & abs(rate - anchorRate) < (5 * abs(anchorRate) + 1e-3);
end
anchorEi  = anchorEi(keep);
anchorOff = anchorOff(keep);
end


function rj = anchorSearch(ei, lo, hi, hVDC, hDetx, hDety, hTof, ...
    eVDC, eDetx, eDety, eTof, icf, tGuess, posTol2, opt)
rj = [];
if hi <= lo
    return
end
rIdx = lo:hi;
dvdc = abs(hVDC(rIdx) - eVDC(ei));
ddx = hDetx(rIdx) * icf - eDetx(ei);
ddy = hDety(rIdx) * icf - eDety(ei);
d2 = ddx .* ddx + ddy .* ddy;
dtof = abs(hTof(rIdx) - eTof(ei) - tGuess);
ok = dvdc < opt.vdcTol & d2 < posTol2 & dtof < opt.tofTol;
if ~any(ok)
    return
end
score = d2 + 0.01 * (dtof .^ 2);
score(~ok) = inf;
[~, j] = min(score);
rj = lo + j - 1;
end
