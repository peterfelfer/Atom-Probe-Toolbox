function tOffset = rhitTOffsetXCorr(rhitTof, eposTof, options)
% RHITTOFFSETXCORR Robust t_offset by cross-correlation of TOF histograms.
%
% tOffset = rhitTOffsetXCorr(rhitTof, eposTof)
% tOffset = rhitTOffsetXCorr(rhitTof, eposTof, 'binNs', 0.5, ...)
%
% Builds a TOF histogram for both streams (in ns) and finds the lag that
% maximises the cross-correlation of H_rhit against H_epos. Works even
% when the streams are NOT 1:1 aligned at the start — only the
% *distribution* of TOFs needs to be representative.
%
% NAME-VALUE OPTIONS:
%   'nSample' - subsample size from each stream (default 500000)
%   'binNs'   - histogram bin width in ns (default 0.5; sets lag resolution)
%   'tofMin'  - lower TOF cutoff in ns (default 0)
%   'tofMax'  - upper TOF cutoff in ns (default 8000)
%
% Returns t_offset = tof_rhit - tof_epos, in ns.
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    rhitTof (:,1) double
    eposTof (:,1) double
    options.nSample (1,1) double = 500000
    options.binNs   (1,1) double = 0.5
    options.tofMin  (1,1) double = 0
    options.tofMax  (1,1) double = 8000
end

nR = numel(rhitTof);
nE = numel(eposTof);

if nR > options.nSample
    rng(42);
    rhitTof = rhitTof(randperm(nR, options.nSample));
end
if nE > options.nSample
    rng(42);
    eposTof = eposTof(randperm(nE, options.nSample));
end

edges = options.tofMin:options.binNs:options.tofMax;
hR = histcounts(rhitTof(rhitTof >= options.tofMin & rhitTof < options.tofMax), edges);
hE = histcounts(eposTof(eposTof >= options.tofMin & eposTof < options.tofMax), edges);

% Cross-correlate via FFT (no Signal Processing Toolbox dependency).
% xcorr(x, y)[k] = sum_n x(n+k) y(n).  Implemented as conv(x, flip(y))
% which has length 2N-1 and lag axis -(N-1) .. (N-1).
N = numel(hR);
c = conv(double(hR), flip(double(hE)));
lags = -(N - 1):(N - 1);
[~, k] = max(c);
lagBins = lags(k);
tOffset = double(lagBins) * options.binNs;
end
