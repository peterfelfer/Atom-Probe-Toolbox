%% Step 6-11: Add ions, ranges, and create visualizations
% Continuing from step 5 (peaks and auto-assignments already done)

clearvars; close all;

% Add toolbox to path
addpath(genpath(pwd));

%% File paths
eposFile = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/recons/recon-v01/default/R56_09162-v01.epos';

[parentDir, ~, ~] = fileparts(eposFile);
outputDir = fullfile(parentDir, 'analysis_output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fprintf('Output directory: %s\n', outputDir);

%% Load pos
fprintf('\n=== Loading data ===\n');
pos = posLoad(eposFile);
fprintf('Loaded %d atoms\n', height(pos));

%% Load isotope table (correct file)
data = load('isotopeTable_naturalAbundances.mat');
isoTable = data.isotopeTable;
fprintf('Loaded isotope table with %d entries\n', height(isoTable));

%% Create color scheme
% Need to create an ion table first for colorSchemeCreate
ionNames = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'Co', 'C', 'Mn', 'Si', ...
            'H', 'O', 'N', 'NiO', 'CrO', 'FeO', 'TiO', 'AlO', 'NiH', 'CrH', 'FeH', 'TiH'};
ionTableForColors = table(categorical(ionNames'), 'VariableNames', {'ionName'});
colors = colorSchemeCreate(ionTableForColors);

%% ionAdd parameters
sumMargin = 0.1;    % margin for peak summing
minAbundance = 0.01; % minimum abundance threshold (1%)

%% Step 6: Create mass spectrum and add ions
fprintf('\n=== Step 6: Adding ions to mass spectrum ===\n');
fig2 = figure('Position', [100, 100, 1600, 800], 'Visible', 'off');
spec = massSpecPlot(pos, 0.01);
xlim([0 100]);
title('Mass Spectrum with Ion Assignments - A718 (Inconel 718)');

% Add ions based on A718 composition
% ionAdd(spec, ion, chargeState, isotopeTable, colorScheme, sumMargin, minAbundance)

% Major elements: Ni, Cr, Fe
ionAdd(spec, 'Ni+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Ni+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Ni+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Cr+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Cr+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Cr+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Fe+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Fe+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Fe+', 3, isoTable, colors, sumMargin, minAbundance);

% Minor elements: Nb, Mo, Ti, Al
ionAdd(spec, 'Nb+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Nb+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Nb+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Mo+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Mo+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Mo+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Ti+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Ti+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Ti+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Al+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Al+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Al+', 3, isoTable, colors, sumMargin, minAbundance);

% Trace elements
ionAdd(spec, 'Co+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Co+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'C+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'C+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Mn+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Mn+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Si+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'Si+', 2, isoTable, colors, sumMargin, minAbundance);

% Hydrogen (important for this hydrogen study)
ionAdd(spec, 'H+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'H2+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'H3+', 1, isoTable, colors, sumMargin, minAbundance);

% Oxygen and nitrogen (possible contaminants)
ionAdd(spec, 'O+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'O+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'N+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'N+', 2, isoTable, colors, sumMargin, minAbundance);

% Molecular ions - oxides
ionAdd(spec, 'NiO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'NiO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'CrO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'CrO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'FeO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'FeO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'TiO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'TiO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'AlO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'AlO+', 2, isoTable, colors, sumMargin, minAbundance);

% Hydrides (important for hydrogen study)
ionAdd(spec, 'NiH+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'NiH+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'CrH+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'CrH+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'FeH+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'FeH+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'TiH+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec, 'TiH+', 2, isoTable, colors, sumMargin, minAbundance);

fprintf('Added ions to mass spectrum\n');

% Save mass spectrum with ions
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ions.png'));
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ions.fig'));
fprintf('Saved mass spectrum with ions\n');

%% Step 7: Add ranges using rangeAddAll
fprintf('\n=== Step 7: Adding ranges ===\n');
rangeAddAll(spec);

% Extract ranges from the mass spectrum
ranges = rangesExtractFromMassSpec(spec);
fprintf('Added %d ranges from mass spectrum\n', height(ranges));

% Save ranges
writetable(ranges, fullfile(outputDir, 'ranges.csv'));
save(fullfile(outputDir, 'ranges.mat'), 'ranges');
fprintf('Saved ranges to CSV and MAT files\n');

% Save mass spectrum with ranges
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ranges.png'));
saveas(fig2, fullfile(outputDir, 'mass_spectrum_with_ranges.fig'));
fprintf('Saved mass spectrum with ranges\n');

%% Create zoomed mass spectrum views
% Low mass region (0-20 Da)
fig_low = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_low = massSpecPlot(pos, 0.01);
xlim([0 20]);
title('Mass Spectrum - Low Mass Region (0-20 Da)');
ionAdd(spec_low, 'H+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'H2+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'H3+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'C+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'C+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'N+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'N+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'O+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'O+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Al+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Al+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Si+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Ti+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Cr+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Fe+', 3, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_low, 'Ni+', 3, isoTable, colors, sumMargin, minAbundance);
saveas(fig_low, fullfile(outputDir, 'mass_spectrum_0-20Da.png'));
close(fig_low);

% Mid mass region (20-40 Da)
fig_mid = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_mid = massSpecPlot(pos, 0.01);
xlim([20 40]);
title('Mass Spectrum - Mid Mass Region (20-40 Da)');
ionAdd(spec_mid, 'Al+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Cr+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Mn+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Fe+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Co+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Ni+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Ti+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_mid, 'Nb+', 3, isoTable, colors, sumMargin, minAbundance);
saveas(fig_mid, fullfile(outputDir, 'mass_spectrum_20-40Da.png'));
close(fig_mid);

% High mass region (40-100 Da)
fig_high = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_high = massSpecPlot(pos, 0.01);
xlim([40 100]);
title('Mass Spectrum - High Mass Region (40-100 Da)');
ionAdd(spec_high, 'Ti+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Cr+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Mn+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Fe+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Co+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Ni+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Nb+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Nb+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Mo+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'Mo+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'NiO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'NiO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'CrO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'CrO+', 2, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'FeO+', 1, isoTable, colors, sumMargin, minAbundance);
ionAdd(spec_high, 'FeO+', 2, isoTable, colors, sumMargin, minAbundance);
saveas(fig_high, fullfile(outputDir, 'mass_spectrum_40-100Da.png'));
close(fig_high);

fprintf('Saved zoomed mass spectrum views\n');

%% Step 8: Apply ranges to pos table
fprintf('\n=== Step 8: Applying ranges to pos table ===\n');
posRanged = rangesFromPos(pos, ranges);
numRanged = sum(~ismissing(posRanged.ion));
fprintf('Ranged atoms: %d out of %d (%.1f%%)\n', ...
    numRanged, height(posRanged), 100 * numRanged / height(posRanged));

% Summary of ranged species
if ismember('atom', posRanged.Properties.VariableNames)
    atomCounts = groupcounts(posRanged, 'atom');
    atomCounts = sortrows(atomCounts, 'GroupCount', 'descend');
    fprintf('\nAtom counts (top 15):\n');
    disp(atomCounts(1:min(15, height(atomCounts)), :));
end

% Summary of ion counts
if ismember('ion', posRanged.Properties.VariableNames)
    ionCounts = groupcounts(posRanged, 'ion');
    ionCounts = sortrows(ionCounts, 'GroupCount', 'descend');
    fprintf('\nIon counts (top 20):\n');
    disp(ionCounts(1:min(20, height(ionCounts)), :));
end

%% Step 9: Create scatter plot visualization
fprintf('\n=== Step 9: Creating scatter plot visualization ===\n');

% Create scatter plot of all ranged ions
fig3 = figure('Position', [100, 100, 1200, 900], 'Visible', 'off');
scatterPlotPosData(posRanged);
title('3D Atom Map - A718 (Inconel 718)');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
view(3);
axis equal;

saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.png'));
saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.fig'));
fprintf('Saved scatter plot of all ions\n');

%% Step 10: Create individual element scatter plots
fprintf('\n=== Step 10: Creating individual element scatter plots ===\n');

if ismember('atom', posRanged.Properties.VariableNames)
    uniqueAtoms = categories(posRanged.atom);
    uniqueAtoms = uniqueAtoms(~cellfun(@isempty, uniqueAtoms));

    % Major and minor elements to plot
    elementsToPlot = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'H', 'O', 'C'};

    for i = 1:length(elementsToPlot)
        elem = elementsToPlot{i};
        if ismember(elem, uniqueAtoms)
            posElem = posRanged(posRanged.atom == elem, :);
            if height(posElem) > 0
                fig_elem = figure('Position', [100, 100, 800, 600], 'Visible', 'off');

                % Subsample if too many points
                if height(posElem) > 100000
                    idx = randperm(height(posElem), 100000);
                    scatter3(posElem.x(idx), posElem.y(idx), posElem.z(idx), 1, '.');
                else
                    scatter3(posElem.x, posElem.y, posElem.z, 1, '.');
                end

                title(sprintf('%s atoms (%d total)', elem, height(posElem)));
                xlabel('x (nm)');
                ylabel('y (nm)');
                zlabel('z (nm)');
                axis equal;
                view(3);

                saveas(fig_elem, fullfile(outputDir, sprintf('scatter_%s.png', elem)));
                close(fig_elem);
                fprintf('Saved scatter plot for %s (%d atoms)\n', elem, height(posElem));
            end
        end
    end
end

%% Step 11: Save all data
fprintf('\n=== Step 11: Saving all data ===\n');

% Save pos and posRanged
save(fullfile(outputDir, 'pos.mat'), 'pos', '-v7.3');
save(fullfile(outputDir, 'posRanged.mat'), 'posRanged', '-v7.3');
fprintf('Saved pos.mat and posRanged.mat\n');

% Calculate and save concentration summary
fprintf('\n=== Concentration Summary ===\n');
if ismember('atom', posRanged.Properties.VariableNames)
    atomCounts = groupcounts(posRanged, 'atom');
    totalRanged = sum(atomCounts.GroupCount);
    atomCounts.Concentration_pct = 100 * atomCounts.GroupCount / totalRanged;
    atomCounts = sortrows(atomCounts, 'GroupCount', 'descend');

    fprintf('\nAtomic concentrations:\n');
    for i = 1:min(15, height(atomCounts))
        fprintf('  %s: %.2f%% (%d atoms)\n', ...
            string(atomCounts.atom(i)), ...
            atomCounts.Concentration_pct(i), ...
            atomCounts.GroupCount(i));
    end

    writetable(atomCounts, fullfile(outputDir, 'concentration_summary.csv'));
    fprintf('\nSaved concentration_summary.csv\n');
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('All outputs saved to: %s\n', outputDir);

% List saved files
fprintf('\nSaved files:\n');
files = dir(outputDir);
for i = 1:length(files)
    if ~files(i).isdir
        fprintf('  - %s (%.1f KB)\n', files(i).name, files(i).bytes/1024);
    end
end

close all;
