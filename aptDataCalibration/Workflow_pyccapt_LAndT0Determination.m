%% L and t0 Determination for pyccapt Data
% Estimate the propagation delay t0 and flight-path length L from
% known isotope peaks using linear regression on time-of-flight data.
%
% The flight-time equation is:
%   tof = t0 + L * sqrt(m_ideal * amu / (2 * e * VDC))
%
% By selecting several peaks with known masses and fitting a line to
% (regression feature, tof), we recover L (slope) and t0 (intercept).

%% 1. Setup
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme.mat');

%% 2. User Parameters
flightPathLength = 110;    % mm (initial estimate)
t0Initial        = 38;     % ns (initial estimate)
pulseMode        = 'voltage';  % 'voltage' or 'laser'
maxMc            = 200;    % Da

%% 3. Load pyccapt Raw HDF5
pos = posLoadPyccapt();

%% 4. Temporal Crop (optional)
% Uncomment to interactively crop:
%   pos = plotExperimentHistory(pos);

%% 5. Fine-Tune t0 (rough)
% Adjust the slider until major peaks are roughly in the right place.
if strcmp(pulseMode, 'voltage')
    VP = pos.VP;
else
    VP = [];
end
[pos, t0] = fineTuneT0(pos, flightPathLength, t0Initial, ...
    'mode', pulseMode, 'VP', VP, 'maxMc', maxMc);

%% 6. Mass Spectrum — Identify Known Peaks
% Use the toolbox mass spectrum and ion identification functions.
spec = massSpecPlot(pos.mc, 0.01, 'normalised');

% Add ions you expect to see in your specimen.
% Adjust these for your material.
ionAdd(spec, 'Al', 1, isotopeTable, colorScheme);
ionAdd(spec, 'Al', 2, isotopeTable, colorScheme);
ionAdd(spec, 'Si', 1, isotopeTable, colorScheme);

%% 7. Build Selected Peaks Table
% For each known peak, define the mass range and the ideal
% mass-to-charge ratio (= atomic weight / charge state).
%
% Edit the values below to match the peaks you identified.
%
% Example for Al and Si:
%   Al-27 (1+): ideal mc = 26.9815 / 1 = 26.9815 Da
%   Al-27 (2+): ideal mc = 26.9815 / 2 = 13.4908 Da
%   Si-28 (1+): ideal mc = 27.9769 / 1 = 27.9769 Da

selectedPeaks = table( ...
    [13.2;  26.7;  27.7], ...   % mcLow
    [13.8;  27.3;  28.3], ...   % mcHigh
    [13.4908; 26.9815; 27.9769], ... % idealMc (weight / charge)
    'VariableNames', {'mcLow', 'mcHigh', 'idealMc'});

disp('Selected peaks for regression:');
disp(selectedPeaks);

%% 8. Determine L and t0
% The function restricts to ions near the detector centre to reduce
% geometric path-length variation, then fits the linear model.
[L, t0Fitted, fitInfo] = determineLAndT0(pos, selectedPeaks, ...
    'showPlot', true, 'detectorWindow', 5);

%% 9. Compare with Fixed-L Model
% Also compute t0 assuming L is known (from instrument geometry):
[~, t0Fixed] = determineLAndT0(pos, selectedPeaks, ...
    'flightPathLength', flightPathLength, 'showPlot', false);
fprintf('Fixed-L model (L = %d mm): t0 = %.2f ns\n', flightPathLength, t0Fixed);

%% 10. Recompute Mass Spectrum with Fitted Parameters
% Use the fitted L and t0 to recalculate mc.
mcArgs = {'mode', pulseMode};
if strcmp(pulseMode, 'voltage')
    mcArgs = [mcArgs, {'VP', pos.VP}];
end
pos.mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, ...
    L, t0Fitted, mcArgs{:});

spec = massSpecPlot(pos.mc, 0.01, 'normalised');
title(sprintf('Calibrated spectrum (L=%.1f mm, t0=%.1f ns)', L, t0Fitted));

fprintf('\nUse these values in your data processing workflow:\n');
fprintf('  flightPathLength = %.2f;  %% mm\n', L);
fprintf('  t0 = %.2f;               %% ns\n', t0Fitted);
