function [hits, calib] = rhitCalibrateFromEpos(hits, epos, options)
% RHITCALIBRATEFROMPOS Calibrate RHIT mc using a matching epos file.
%
% [hits, calib] = rhitCalibrateFromEpos(hits, epos)
% [hits, calib] = rhitCalibrateFromEpos(hits, epos, 'bowlKind', 'xy2d', 'bowlDegree', 4)
%
% Voltage + bowl calibration anchored on Cameca's calibrated EPOS mc.
%
% PIPELINE:
%   1. Robust t_offset estimate from cross-correlation of TOF histograms
%      (rhitTOffsetXCorr) — works even when streams aren't 1:1 aligned.
%   2. Drift-tracking event matcher (rhitMatchEvents) — handles cumulative
%      ~0.03%/event drift between EPOS and RHIT.
%   3. 2D xy polynomial bowl fit (rhitFitBowl2D, default degree 4) — the
%      LEAP bowl is asymmetric; a radial-only fit plateaus at ~2.6%
%      residual where 2D drops below 0.1%.
%   4. Apply via rhitApplyCalibration to populate hits.mc.
%
% INPUTS:
%   hits  - Table from rhitLoad (same run as epos)
%   epos  - Table from posLoad of a matching .epos file
%
% NAME-VALUE OPTIONS:
%   'bowlKind'    - 'xy2d' (default) or 'radial'
%   'bowlDegree'  - default 4
%   'stride'      - matcher EPOS-stride (default 200)
%   'verbose'     - print progress (default true)
%   Plus any rhitMatchEvents / rhitFitBowl2D options pass-through.
%
% OUTPUTS:
%   hits  - input table with mc replaced by the calibrated mc
%   calib - calibration struct usable by rhitApplyCalibration:
%             tOffsetNs    - TOF offset (ns)
%             bowl         - struct {kind, degree, coeffs, terms}
%             ICF          - epos.detx / hits.detx (mm/mm)
%             nMatched     - matched-event count
%             fitResidualStd / fitResidualPctC0
%             flightPathMm / detectorHalfsizeMm / t0RhitNs / kf
%             notes
%
% EXAMPLE:
%   [rhitH, ~, meta] = rhitLoad('R56_03063.RHIT');
%   epos = posLoad('R56_03063-v01.epos');
%   [rhitH, calib] = rhitCalibrateFromEpos(rhitH, epos);
%   massSpecPlot(rhitH, 0.01, 'normalised');
%
% See also: rhitLoad, rhitMatchEvents, rhitFitBowl2D, rhitApplyCalibration,
%           rhitTOffsetXCorr, rhitPoolBowl, rhitRefineTOffset
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    epos table
    options.bowlKind (1,:) char = 'xy2d'
    options.bowlDegree (1,1) double = 4
    options.stride (1,1) double = 200
    options.scanBack (1,1) double = 50
    options.scanFwd (1,1) double = 2000
    options.fallbackFwd (1,1) double = 50000
    options.vdcTol (1,1) double = 0.5
    options.posTol (1,1) double = 0.15
    options.tofTol (1,1) double = 1.5
    options.driftAnchors (1,1) double = 200
    options.cOutlierFrac (1,1) double = 0.2
    options.tOffsetGuess (1,1) double = NaN
    options.icf (1,1) double = NaN
    options.vdcMin (1,1) double = -inf
    options.vdcMax (1,1) double = inf
    options.detRadiusMaxMm (1,1) double = inf
    options.verbose (1,1) logical = true
end

nRHIT = height(hits);
nEpos = height(epos);
if options.verbose
    fprintf('Calibrating: %d RHIT events, %d epos events\n', nRHIT, nEpos);
end

% --- 1. Robust t_offset and ICF ---
tGuess = options.tOffsetGuess;
if isnan(tGuess)
    tGuess = rhitTOffsetXCorr(double(hits.tof), double(epos.tof));
    if options.verbose
        fprintf('  t_offset (xcorr) seed: %.2f ns\n', tGuess);
    end
end

% --- 2. Match events ---
[matchedE, matchedR, matchInfo] = rhitMatchEvents(hits, epos, ...
    'tOffsetGuess', tGuess, ...
    'icf', options.icf, ...
    'stride', options.stride, ...
    'scanBack', options.scanBack, ...
    'scanFwd', options.scanFwd, ...
    'fallbackFwd', options.fallbackFwd, ...
    'vdcTol', options.vdcTol, ...
    'posTol', options.posTol, ...
    'tofTol', options.tofTol, ...
    'driftAnchors', options.driftAnchors, ...
    'vdcMin', options.vdcMin, ...
    'vdcMax', options.vdcMax, ...
    'detRadiusMaxMm', options.detRadiusMaxMm, ...
    'verbose', options.verbose);

if options.verbose
    fprintf('  Matched %d events (ICF=%.4f)\n', numel(matchedE), matchInfo.icf);
end

if numel(matchedE) < 100
    warning('rhitCalibrateFromEpos:fewMatches', ...
        'Only %d matches; bowl fit will be unstable.', numel(matchedE));
end

% --- 3. Fit bowl + t_offset on matches ---
[tOffset, bowl, fitInfo] = rhitFitBowl2D(hits, epos, matchedE, matchedR, ...
    'bowlKind', options.bowlKind, ...
    'bowlDegree', options.bowlDegree, ...
    'cOutlierFrac', options.cOutlierFrac);

if options.verbose
    fprintf('  Fit:  t_offset=%.4f ns (std %.4f)\n', ...
        fitInfo.tOffsetNs, fitInfo.tOffsetStdNs);
    fprintf('        bowl=%s deg=%d (%d terms)\n', ...
        bowl.kind, bowl.degree, numel(bowl.coeffs));
    fprintf('        residual std=%.3e (%.3f%% of C0)\n', ...
        fitInfo.fitResidualStd, fitInfo.fitResidualPctC0);
    fnames = fieldnames(fitInfo.LeffMmByRMm);
    for k = 1:numel(fnames)
        fprintf('        L_eff(%s mm) = %.2f mm\n', ...
            erase(fnames{k}, 'r'), fitInfo.LeffMmByRMm.(fnames{k}));
    end
end

% --- 4. Apply ---
calib = struct();
calib.tOffsetNs = tOffset;
calib.bowl = bowl;
calib.ICF = matchInfo.icf;
calib.nMatched = numel(matchedE);
calib.fitResidualStd = fitInfo.fitResidualStd;
calib.fitResidualPctC0 = fitInfo.fitResidualPctC0;
calib.tOffsetStdNs = fitInfo.tOffsetStdNs;

% Legacy field aliases (so old code reading calib.t_offset / calib.C_poly
% keeps working).  C_poly is only populated for radial fits; for an
% xy2d bowl we leave it absent.
calib.t_offset = tOffset;
calib.residual_std = fitInfo.fitResidualStd;
calib.clock_info = 'tof in ns, subtract t_offset before use';
if strcmpi(bowl.kind, 'radial')
    calib.C_poly = bowl.coeffs;
end

% Optional metadata: rhitLoad returns instrument params as a separate
% output (the third), but if the caller stuffed it into UserData we
% capture it here.  Otherwise these fields are simply absent.
if istable(hits) && ~isempty(hits.Properties.UserData) ...
        && isstruct(hits.Properties.UserData) ...
        && isfield(hits.Properties.UserData, 'instrumentParams')
    ip = hits.Properties.UserData.instrumentParams;
    if isfield(ip, 'flight_path_mm')
        calib.flightPathMm = ip.flight_path_mm;
    end
    if isfield(ip, 'detector_halfsize_mm')
        calib.detectorHalfsizeMm = ip.detector_halfsize_mm;
    end
    if isfield(ip, 't0_ns')
        calib.t0RhitNs = ip.t0_ns;
    end
end
if ismember('Vref', hits.Properties.VariableNames) && hits.VDC(1) > 0
    calib.kf = double(hits.Vref(1)) / double(hits.VDC(1));
else
    calib.kf = 1.03;
end
calib.notes = sprintf('nRHIT=%d, nEpos=%d', nRHIT, nEpos);

hits = rhitApplyCalibration(hits, calib);
end
