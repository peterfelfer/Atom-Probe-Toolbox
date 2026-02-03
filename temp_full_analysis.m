%% Full Analysis Script for A718 (Inconel 718) EPOS data
% This script loads the EPOS file, creates a mass spectrum, finds peaks,
% assigns ions, adds ranges, and creates a scatter plot visualization.

clearvars; close all;

% Add toolbox to path
addpath(genpath(pwd));

%% File paths
eposFile = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/recons/recon-v01/default/R56_09162-v01.epos';

% Create output folder
[parentDir, ~, ~] = fileparts(eposFile);
outputDir = fullfile(parentDir, 'analysis_output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fprintf('Output directory: %s\n', outputDir);

%% Step 1: Load EPOS file
fprintf('\n=== Step 1: Loading EPOS file ===\n');
pos = posLoad(eposFile);
fprintf('Loaded %d atoms\n', height(pos));
fprintf('Mass range: %.2f - %.2f Da\n', min(pos.mc), max(pos.mc));

%% Step 2: Create mass spectrum
fprintf('\n=== Step 2: Creating mass spectrum ===\n');
fig1 = figure('Position', [100, 100, 1400, 600], 'Visible', 'on');
spec = massSpecPlot(pos, 0.01);
xlim([0 100]);
title('Mass Spectrum - A718 (Inconel 718)');

% Save initial mass spectrum
saveas(fig1, fullfile(outputDir, 'mass_spectrum_overview.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_overview.fig'));
fprintf('Saved mass spectrum overview\n');

%% Step 3: Find peaks using new massSpecFindPeaks function
fprintf('\n=== Step 3: Finding peaks (using new massSpecFindPeaks) ===\n');
[peakTable, peakInfo] = massSpecFindPeaks(pos, 'binWidth', 0.02, 'mcRange', [0 100]);
fprintf('Found %d peaks\n', height(peakTable));

% Display top peaks by prominence
if height(peakTable) > 0
    sortedPeaks = sortrows(peakTable, 'prominence', 'descend');
    fprintf('\nTop 20 peaks by prominence:\n');
    disp(sortedPeaks(1:min(20, height(sortedPeaks)), :));
end

%% Step 4: Auto-assign ions using massSpecAutoAssignIons
fprintf('\n=== Step 4: Auto-assigning ions (using new massSpecAutoAssignIons) ===\n');

% A718 (Inconel 718) expected elements:
% Major: Ni (~50-55%), Cr (~17-21%), Fe (~17%)
% Minor: Nb (~5%), Mo (~3%), Ti (~1%), Al (~0.5%)
% Trace: Co, C, Mn, Si, B, H
elements = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'Co', 'C', 'Mn', 'Si', 'B', 'H', 'O', 'N'};

[ionTable, assignInfo] = massSpecAutoAssignIons(peakTable, elements, ...
    'maxAtoms', 2, ...
    'maxCharge', 3, ...
    'conservative', true);

fprintf('Assigned %d ions\n', height(ionTable));
if height(ionTable) > 0
    disp(ionTable);
end

%% Step 5: Generate ranges using massSpecAutoRanges
fprintf('\n=== Step 5: Generating ranges (using new massSpecAutoRanges) ===\n');
if ~isempty(assignInfo.assigned)
    [rangeTable, rangeInfo] = massSpecAutoRanges(pos, peakTable, assignInfo.assigned, ...
        'binWidth', 0.02);
    fprintf('Generated %d ranges\n', height(rangeTable));

    if height(rangeTable) > 0
        disp(rangeTable(1:min(30, height(rangeTable)), :));
    end
else
    fprintf('No ions assigned, cannot generate ranges.\n');
    rangeTable = table();
end

%% Step 6: Add ions and ranges to mass spectrum manually for major elements
fprintf('\n=== Step 6: Adding ions and ranges to mass spectrum ===\n');

% Load isotope table and color scheme
isoTable = readtable('isotopeTable.csv');
colors = colorSchemeCreate();

% Create new figure for mass spectrum with ions
fig2 = figure('Position', [100, 100, 1600, 800], 'Visible', 'on');
spec2 = massSpecPlot(pos, 0.01);
xlim([0 100]);
title('Mass Spectrum with Ion Assignments - A718');

% Add ions based on expected A718 composition and detected peaks
% Major elements
ionAdd(spec2, 'Ni+', [1 2 3], isoTable, colors);
ionAdd(spec2, 'Cr+', [1 2 3], isoTable, colors);
ionAdd(spec2, 'Fe+', [1 2 3], isoTable, colors);

% Minor elements
ionAdd(spec2, 'Nb+', [1 2 3], isoTable, colors);
ionAdd(spec2, 'Mo+', [1 2 3], isoTable, colors);
ionAdd(spec2, 'Ti+', [1 2 3], isoTable, colors);
ionAdd(spec2, 'Al+', [1 2 3], isoTable, colors);

% Trace elements
ionAdd(spec2, 'Co+', [1 2], isoTable, colors);
ionAdd(spec2, 'C+', [1 2], isoTable, colors);
ionAdd(spec2, 'Mn+', [1 2], isoTable, colors);
ionAdd(spec2, 'Si+', [1 2], isoTable, colors);

% Hydrogen (important for this study)
ionAdd(spec2, 'H+', 1, isoTable, colors);
ionAdd(spec2, 'H2+', 1, isoTable, colors);
ionAdd(spec2, 'H3+', 1, isoTable, colors);

% Oxygen and nitrogen (possible contaminants)
ionAdd(spec2, 'O+', [1 2], isoTable, colors);
ionAdd(spec2, 'N+', [1 2], isoTable, colors);

% Molecular ions
ionAdd(spec2, 'NiO+', [1 2], isoTable, colors);
ionAdd(spec2, 'CrO+', [1 2], isoTable, colors);
ionAdd(spec2, 'FeO+', [1 2], isoTable, colors);
ionAdd(spec2, 'TiO+', [1 2], isoTable, colors);
ionAdd(spec2, 'AlO+', [1 2], isoTable, colors);

% Hydrides (important for hydrogen study)
ionAdd(spec2, 'NiH+', [1 2], isoTable, colors);
ionAdd(spec2, 'CrH+', [1 2], isoTable, colors);
ionAdd(spec2, 'FeH+', [1 2], isoTable, colors);
ionAdd(spec2, 'TiH+', [1 2], isoTable, colors);

% Save mass spectrum with ions
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ions.png'));
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ions.fig'));
fprintf('Saved mass spectrum with ions\n');

%% Step 7: Add ranges using rangeAddAll
fprintf('\n=== Step 7: Adding ranges ===\n');
rangeAddAll(spec2);

% Extract ranges from the mass spectrum
ranges = rangesExtractFromMassSpec(spec2);
fprintf('Added %d ranges from mass spectrum\n', height(ranges));

% Save ranges table
writetable(ranges, fullfile(outputDir, 'ranges.csv'));
save(fullfile(outputDir, 'ranges.mat'), 'ranges');
fprintf('Saved ranges to CSV and MAT files\n');

% Save mass spectrum with ranges
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ranges.png'));
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ranges.fig'));
fprintf('Saved mass spectrum with ranges\n');

%% Step 8: Apply ranges to pos table
fprintf('\n=== Step 8: Applying ranges to pos table ===\n');
posRanged = rangesFromPos(pos, ranges);
fprintf('Ranged atoms: %d out of %d (%.1f%%)\n', ...
    sum(~ismissing(posRanged.ion)), height(posRanged), ...
    100 * sum(~ismissing(posRanged.ion)) / height(posRanged));

% Summary of ranged species
if ismember('atom', posRanged.Properties.VariableNames)
    atomCounts = groupcounts(posRanged, 'atom');
    atomCounts = sortrows(atomCounts, 'GroupCount', 'descend');
    fprintf('\nAtom counts:\n');
    disp(atomCounts(1:min(15, height(atomCounts)), :));
end

%% Step 9: Create scatter plot visualization
fprintf('\n=== Step 9: Creating scatter plot visualization ===\n');

% Create scatter plot of all ions
fig3 = figure('Position', [100, 100, 1200, 900], 'Visible', 'on');
scatterPlotPosData(posRanged);
title('3D Atom Map - A718 (Inconel 718)');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
view(3);
axis equal;

% Save scatter plot
saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.png'));
saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.fig'));
fprintf('Saved scatter plot\n');

%% Step 10: Create individual element scatter plots
fprintf('\n=== Step 10: Creating individual element scatter plots ===\n');

% Get unique atoms
if ismember('atom', posRanged.Properties.VariableNames)
    uniqueAtoms = categories(posRanged.atom);
    uniqueAtoms = uniqueAtoms(~cellfun(@isempty, uniqueAtoms));

    % Create scatter plots for major elements
    majorElements = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al'};

    for i = 1:length(majorElements)
        elem = majorElements{i};
        if ismember(elem, uniqueAtoms)
            fig_elem = figure('Position', [100, 100, 800, 600], 'Visible', 'on');
            posElem = posRanged(posRanged.atom == elem, :);
            if height(posElem) > 0
                scatter3(posElem.x, posElem.y, posElem.z, 1, '.');
                title(sprintf('%s atoms (%d)', elem, height(posElem)));
                xlabel('x (nm)');
                ylabel('y (nm)');
                zlabel('z (nm)');
                axis equal;
                view(3);
                saveas(fig_elem, fullfile(outputDir, sprintf('scatter_%s.png', elem)));
                close(fig_elem);
                fprintf('Saved scatter plot for %s\n', elem);
            end
        end
    end
end

%% Step 11: Save all data
fprintf('\n=== Step 11: Saving all data ===\n');

% Save pos and posRanged
save(fullfile(outputDir, 'pos.mat'), 'pos');
save(fullfile(outputDir, 'posRanged.mat'), 'posRanged');

% Save peak analysis results
if exist('peakTable', 'var')
    writetable(peakTable, fullfile(outputDir, 'peaks_detected.csv'));
    save(fullfile(outputDir, 'peakAnalysis.mat'), 'peakTable', 'peakInfo');
end

% Save ion assignment results
if exist('ionTable', 'var') && height(ionTable) > 0
    save(fullfile(outputDir, 'ionAssignment.mat'), 'ionTable', 'assignInfo');
end

% Save auto-generated ranges
if exist('rangeTable', 'var') && height(rangeTable) > 0
    writetable(rangeTable, fullfile(outputDir, 'ranges_auto.csv'));
    save(fullfile(outputDir, 'rangesAuto.mat'), 'rangeTable', 'rangeInfo');
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('All outputs saved to: %s\n', outputDir);

% List saved files
fprintf('\nSaved files:\n');
files = dir(outputDir);
for i = 1:length(files)
    if ~files(i).isdir
        fprintf('  - %s\n', files(i).name);
    end
end
