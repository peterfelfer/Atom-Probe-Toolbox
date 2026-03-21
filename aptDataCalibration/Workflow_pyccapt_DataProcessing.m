%% pyccapt Data Processing Workflow
% End-to-end calibration workflow for pyccapt raw HDF5 data.
%
% This workflow combines pyccapt-specific calibration functions with
% standard Atom Probe Toolbox functions for ion identification,
% reconstruction, and visualisation.

%% 1. Setup
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme.mat');

%% 2. User Parameters
% Edit these to match your instrument and experiment.
flightPathLength = 110;    % mm
t0Initial        = 38;     % ns (propagation delay estimate)
pulseMode        = 'voltage';  % 'voltage' or 'laser'
maxMc            = 200;    % Da (maximum mass-to-charge to display)
detEff           = 0.57;   % detection efficiency

%% 3. Load pyccapt Raw HDF5
% A file dialog opens if no path is given.
pos = posLoadPyccapt();
% Alternatively, specify a path:
%   pos = posLoadPyccapt('/path/to/measurement.h5');
% Or use the auto-detecting canonical loader:
%   pos = posLoad('measurement.h5');

%% 4. Temporal Crop (Experiment History)
% A 2-D histogram of ion index vs. TOF is shown with a voltage overlay.
% Draw a rectangle around the stable evaporation region, then
% double-click inside it to confirm.
pos = plotExperimentHistory(pos);

%% 5. Spatial Crop (Field Desorption Map)
% A 2-D histogram of the detector plane is shown.
% Draw a circle around the region of interest, then double-click to confirm.
pos = plotFieldDesorptionMap(pos);

%% 6. Fine-Tune t0
% An interactive slider lets you adjust t0 while watching the mass
% spectrum update in real time. Press "Accept" when a known reference
% peak aligns with its expected position.
if strcmp(pulseMode, 'voltage')
    VP = pos.VP;
else
    VP = [];
end
[pos, t0] = fineTuneT0(pos, flightPathLength, t0Initial, ...
    'mode', pulseMode, 'VP', VP, 'maxMc', maxMc);

%% 7. Select Reference Peak for Calibration
% Visualise the uncalibrated mass spectrum using the standard toolbox
% function, then click two points to define the reference peak range.
spec = massSpecPlot(pos.mc, 0.01, 'normalised');

fprintf('Click two points to define the reference peak range [lo, hi].\n');
[gx, ~] = ginput(2);
peakRange = sort(gx(:))';
fprintf('Reference peak range: [%.3f, %.3f] Da\n', peakRange(1), peakRange(2));
close(gcf);

%% 8. Iterative Voltage + Bowl Correction
maxIterations = 10;
for iter = 1:maxIterations
    fprintf('\n--- Calibration iteration %d / %d ---\n', iter, maxIterations);
    pos = voltageCorrection(pos, peakRange, 'showPlot', false);
    pos = bowlCorrection(pos, peakRange, 'showPlot', false);
end
fprintf('\nCalibration complete.\n');

%% 9. Calibrated Mass Spectrum
% Use the standard toolbox function to inspect the calibrated spectrum.
spec = massSpecPlot(pos.mc, 0.01, 'normalised');

%% 10. Ion Identification (Toolbox)
% Add known ions to the spectrum. Adjust elements and charge states
% for your specimen.
ionAdd(spec, 'Al', 1, isotopeTable, colorScheme);
ionAdd(spec, 'Al', 2, isotopeTable, colorScheme);
ionAdd(spec, 'Si', 1, isotopeTable, colorScheme);
% Add more ions as needed:
%   ionAdd(spec, 'Fe', 2, isotopeTable, colorScheme);

%% 11. Range Definition (Toolbox)
% Automatically add ranges for all identified ions.
rangeAddAll(spec, colorScheme, 3, true);
% Or define ranges manually:
%   rangeAdd(spec, colorScheme);

%% 12. Allocate Ions to Ranges (Toolbox)
rangeTable = rangesExtractFromMassSpec(spec);
pos = posAllocateRange(pos, rangeTable, 'decompose');

%% 13. Composition (Toolbox)
conc = posCalculateConcentrationSimple(pos, detEff, {'unranged'}, '', 'mode', 'atomic');
disp(conc);

%% 14. 3-D Reconstruction (Toolbox)
% Compute tip radius evolution and reconstruct.
kf = 3.3;          % field factor (V/nm)
fieldEvap = 25;     % field evaporation field (V/nm)
ICF = 1.65;         % image compression factor

radiusEvolution = pos.VDC ./ (kf * fieldEvap);
[~, rad] = cart2pol(pos.detx, pos.dety);
effectiveDetArea = max(rad)^2 * pi;

% Atomic volume from isotopeTable (use dominant species)
atomicVolume = 0.0166;  % nm^3 (aluminium default — adjust for your material)
ionVolumes = atomicVolume * ones(height(pos), 1);

pos = posReconstruct3DGeiser(pos, flightPathLength, detEff, ...
    effectiveDetArea, ionVolumes, radiusEvolution, ICF);

%% 15. 3-D Visualisation (Toolbox)
scatterPlotPosWidget(pos, colorScheme);

%% 16. Save Results (Toolbox)
% Save to HDF5 database:
%   posTableAddToHDF5('results.h5', pos);
% Or export to other formats using the toolbox's export functions.
