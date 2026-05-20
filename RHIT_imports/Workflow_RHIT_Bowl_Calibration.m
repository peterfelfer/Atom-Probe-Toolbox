%% Workflow: RHIT 2D bowl + voltage calibration anchored on EPOS
%
% Demonstrates:
%   1. Per-run calibration from a paired RHIT/EPOS file
%   2. Pooling per-run calibrations into one global instrument bowl
%   3. Re-applying the global bowl with cross-correlation t_offset to
%      runs without (or with broken) EPOS matching
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

% Add the toolbox root and RHIT_imports to the MATLAB path. If you've
% already done this, the addpath calls are no-ops.
thisFile = mfilename('fullpath');
toolboxRoot = fileparts(fileparts(thisFile));     % parent of RHIT_imports
addpath(toolboxRoot);                              % posLoad, massSpecPlot, ...
addpath(fullfile(toolboxRoot, 'RHIT_imports'));    % rhitLoad, rhit* helpers

%% ---------- 1. Per-run calibration ----------------------------------

[hits, hist, meta] = rhitLoad('R56_03063.RHIT');
epos = posLoad('R56_03063-v01.epos');

% Stash instrument params on the hits table so the calibration struct
% can pick them up automatically (optional convenience).
hits.Properties.UserData.instrumentParams = meta.instrumentParams;

[hits, calib_03063] = rhitCalibrateFromEpos(hits, epos);
% calib_03063 now contains:
%   .tOffsetNs
%   .bowl              (struct: kind='xy2d', degree=4, coeffs(15), terms(15x2))
%   .ICF
%   .nMatched
%   .fitResidualPctC0  (typically 0.02-0.10% on this LEAP)
%   ...

% Compare against the EPOS-stored Cameca calibration:
massSpecPlot(hits, 0.01, 'normalised');
hold on;
[hE, edE] = histcounts(epos.mc, 'BinWidth', 0.01);
plot(edE(1:end-1), hE / max(hE), 'r-');
legend('our 2D bowl', 'Cameca EPOS mc');
hold off;


%% ---------- 2. Pool many runs into a global bowl --------------------

% Run rhitCalibrateFromEpos on several runs from the same instrument,
% collect their calibration structs, then pool.

runs = {'R56_01274', 'R56_01429', 'R56_01515', 'R56_01554', ...
        'R56_02509', 'R56_03057', 'R56_03063'};
calibs = cell(1, numel(runs));
for k = 1:numel(runs)
    fprintf('=== %s ===\n', runs{k});
    [h, ~, m] = rhitLoad([runs{k}, '.RHIT']);
    h.Properties.UserData.instrumentParams = m.instrumentParams;
    e = posLoad([runs{k}, '-v01.epos']);
    [~, calibs{k}] = rhitCalibrateFromEpos(h, e);
end

% Drop runs with high residual or few matches before pooling
trustworthy = cellfun(@(c) c.fitResidualPctC0 < 0.1 && c.nMatched > 25000, ...
    calibs);
globalBowl = rhitPoolBowl(calibs(trustworthy), 'bowlDegree', 4);
% globalBowl.kind = 'xy2d', .degree = 4, .coeffs (15x1)
% Bowl shape is reproduced from 7 runs to ~0.1% per-cell std.


%% ---------- 3. Apply global bowl + 1-param t_offset to a new RHIT ----

% Useful when:  (a) the EPOS-anchored matcher fails for a run,
%               (b) you want to calibrate a RHIT that has no EPOS at all.
% You provide a reference spectrum to lock t_offset against (typically
% an EPOS mc histogram from any run on the same instrument, or a
% template you build).

[hits2, ~, ~] = rhitLoad('R56_01700.RHIT');

% Build a reference from any matched EPOS — here we re-use one we have:
epos_ref = posLoad('R56_01515-v01.epos');
edges = 0:0.01:80;
refSpectrum = struct('edges', edges, ...
    'values', histcounts(epos_ref.mc, edges));

% Histogram-based starting estimate, then refine against the reference
tGuess = rhitTOffsetXCorr(double(hits2.tof), double(epos_ref.tof));
[tRefined, ~] = rhitRefineTOffset(hits2, refSpectrum, globalBowl, tGuess);

% Build the run-specific calibration struct and apply
calib2 = struct();
calib2.tOffsetNs = tRefined;
calib2.bowl = globalBowl;
hits2 = rhitApplyCalibration(hits2, calib2);

massSpecPlot(hits2, 0.01, 'normalised');


%% ---------- 4. Save / reload calibrations --------------------------

% Save:  save('R56_03063_calib.mat', '-struct', 'calib_03063');
% Save the global bowl:  save('LEAP_globalBowl.mat', '-struct', 'globalBowl');

% Reload + apply to any RHIT from the same instrument:
%   gb = load('LEAP_globalBowl.mat');
%   gb_struct = struct('kind', gb.kind, 'degree', gb.degree, ...
%       'coeffs', gb.coeffs, 'terms', gb.terms);
%   calib = struct('tOffsetNs', tRefined, 'bowl', gb_struct);
%   hits = rhitApplyCalibration(hits, calib);
