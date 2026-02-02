%% Step 2: Create mass spectrum and identify ions for A718 (Inconel 718)
% A718 composition: Ni-Cr-Fe-Nb-Mo-Ti-Al superalloy

% Add toolbox to path
addpath(genpath(pwd));

% Output directory
outputDir = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/recons/recon-v01/default/analysis_output';

% Load previously saved pos
load(fullfile(outputDir, 'pos.mat'), 'pos');
fprintf('Loaded pos with %d atoms\n', height(pos));

%% Create mass spectrum figure
fprintf('Creating mass spectrum...\n');
fig1 = figure('Position', [100, 100, 1400, 600], 'Visible', 'off');
spec = massSpecPlot(pos, 'binWidth', 0.01);
xlim([0 100]);  % Focus on main range first
title('Mass Spectrum - A718 (Inconel 718)');
xlabel('Mass-to-Charge (Da)');
ylabel('Counts');

% Save initial mass spectrum
saveas(fig1, fullfile(outputDir, 'mass_spectrum_overview.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_overview.fig'));
fprintf('Saved mass spectrum overview\n');

%% Load isotope table
isoTable = readtable('isotopeTable.csv');

%% Create color scheme
colors = colorSchemeCreate();

%% Add ions for A718 (Inconel 718) - Nickel-based superalloy
% Major elements: Ni, Cr, Fe, Nb, Mo, Ti, Al
% Minor elements: Co, C, Mn, Si, B, H

fprintf('Adding ions to mass spectrum...\n');

% Hydrogen (if present - important for this study on hydrogen)
ionAdd(spec, 'H+', 1, isoTable, colors);
ionAdd(spec, 'H2+', 1, isoTable, colors);

% Carbon
ionAdd(spec, 'C+', 1, isoTable, colors);
ionAdd(spec, 'C+', 2, isoTable, colors);

% Aluminum
ionAdd(spec, 'Al+', 1, isoTable, colors);
ionAdd(spec, 'Al+', 2, isoTable, colors);
ionAdd(spec, 'Al+', 3, isoTable, colors);

% Silicon (minor)
ionAdd(spec, 'Si+', 1, isoTable, colors);
ionAdd(spec, 'Si+', 2, isoTable, colors);

% Titanium
ionAdd(spec, 'Ti+', 1, isoTable, colors);
ionAdd(spec, 'Ti+', 2, isoTable, colors);
ionAdd(spec, 'Ti+', 3, isoTable, colors);

% Chromium
ionAdd(spec, 'Cr+', 1, isoTable, colors);
ionAdd(spec, 'Cr+', 2, isoTable, colors);
ionAdd(spec, 'Cr+', 3, isoTable, colors);

% Manganese (minor)
ionAdd(spec, 'Mn+', 1, isoTable, colors);
ionAdd(spec, 'Mn+', 2, isoTable, colors);

% Iron
ionAdd(spec, 'Fe+', 1, isoTable, colors);
ionAdd(spec, 'Fe+', 2, isoTable, colors);
ionAdd(spec, 'Fe+', 3, isoTable, colors);

% Cobalt (minor)
ionAdd(spec, 'Co+', 1, isoTable, colors);
ionAdd(spec, 'Co+', 2, isoTable, colors);

% Nickel (major element)
ionAdd(spec, 'Ni+', 1, isoTable, colors);
ionAdd(spec, 'Ni+', 2, isoTable, colors);
ionAdd(spec, 'Ni+', 3, isoTable, colors);

% Niobium
ionAdd(spec, 'Nb+', 1, isoTable, colors);
ionAdd(spec, 'Nb+', 2, isoTable, colors);
ionAdd(spec, 'Nb+', 3, isoTable, colors);

% Molybdenum
ionAdd(spec, 'Mo+', 1, isoTable, colors);
ionAdd(spec, 'Mo+', 2, isoTable, colors);
ionAdd(spec, 'Mo+', 3, isoTable, colors);

% Boron (minor, sometimes added)
ionAdd(spec, 'B+', 1, isoTable, colors);
ionAdd(spec, 'B+', 2, isoTable, colors);

%% Save mass spectrum with ions
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ions.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ions.fig'));
fprintf('Saved mass spectrum with ions\n');

%% Save the spec handle for later use
save(fullfile(outputDir, 'spec.mat'), 'spec', 'isoTable', 'colors');

fprintf('Step 2 complete!\n');
