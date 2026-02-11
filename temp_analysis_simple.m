%% Simplified Analysis Script for A718 EPOS data
% Uses only elemental ions to avoid potential issues with complex ions

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

%% Load isotope table
data = load('isotopeTable_naturalAbundances.mat');
isoTable = data.isotopeTable;
fprintf('Loaded isotope table with %d entries\n', height(isoTable));

%% Create color scheme
ionNames = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'Co', 'C', 'Mn', 'Si', 'H', 'O', 'N'};
ionTableForColors = table(categorical(ionNames'), 'VariableNames', {'ionName'});
colors = colorSchemeCreate(ionTableForColors);

%% ionAdd parameters
sumMargin = 0.1;
minAbundance = 0.01;

%% Step 1: Create mass spectrum and add ions
fprintf('\n=== Creating mass spectrum and adding ions ===\n');
fig1 = figure('Position', [100, 100, 1600, 800], 'Visible', 'off');
spec = massSpecPlot(pos, 0.01);
xlim([0 100]);
title('Mass Spectrum with Ion Assignments - A718 (Inconel 718)');

% Add elemental ions only (avoiding complex ions for now)
elementsAndCharges = {
    'Ni', [1 2 3];
    'Cr', [1 2 3];
    'Fe', [1 2 3];
    'Nb', [1 2 3];
    'Mo', [1 2 3];
    'Ti', [1 2 3];
    'Al', [1 2 3];
    'Co', [1 2];
    'C', [1 2];
    'Mn', [1 2];
    'Si', [1 2];
    'H', 1;
    'O', [1 2];
    'N', [1 2];
};

for i = 1:size(elementsAndCharges, 1)
    elem = elementsAndCharges{i, 1};
    charges = elementsAndCharges{i, 2};
    for cs = charges
        try
            ionAdd(spec, [elem '+'], cs, isoTable, colors, sumMargin, minAbundance);
            fprintf('Added %s%d+\n', elem, cs);
        catch ME
            fprintf('Warning: Could not add %s%d+: %s\n', elem, cs, ME.message);
        end
    end
end

% Try to add H2 and H3
try
    ionAdd(spec, 'H2+', 1, isoTable, colors, sumMargin, minAbundance);
    fprintf('Added H2+\n');
catch ME
    fprintf('Warning: Could not add H2+: %s\n', ME.message);
end

try
    ionAdd(spec, 'H3+', 1, isoTable, colors, sumMargin, minAbundance);
    fprintf('Added H3+\n');
catch ME
    fprintf('Warning: Could not add H3+: %s\n', ME.message);
end

fprintf('Finished adding ions to mass spectrum\n');

% Save mass spectrum with ions
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ions.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ions.fig'));
fprintf('Saved mass spectrum with ions\n');

%% Step 2: Add ranges
fprintf('\n=== Adding ranges ===\n');
rangeAddAll(spec);

% Extract ranges from mass spectrum
ranges = rangesExtractFromMassSpec(spec);
fprintf('Added %d ranges from mass spectrum\n', height(ranges));

% Save ranges
writetable(ranges, fullfile(outputDir, 'ranges.csv'));
save(fullfile(outputDir, 'ranges.mat'), 'ranges');
fprintf('Saved ranges\n');

% Save mass spectrum with ranges
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ranges.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ranges.fig'));
fprintf('Saved mass spectrum with ranges\n');

%% Step 3: Create zoomed views
fprintf('\n=== Creating zoomed mass spectrum views ===\n');

% Low mass (0-20 Da)
fig_low = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_low = massSpecPlot(pos, 0.01);
xlim([0 20]);
title('Mass Spectrum - Low Mass Region (0-20 Da)');
for elem = {'H', 'C', 'N', 'O', 'Al', 'Si', 'Ti', 'Cr', 'Fe', 'Ni'}
    for cs = 1:3
        try
            ionAdd(spec_low, [elem{1} '+'], cs, isoTable, colors, sumMargin, minAbundance);
        catch
        end
    end
end
saveas(fig_low, fullfile(outputDir, 'mass_spectrum_0-20Da.png'));
close(fig_low);

% Mid mass (20-40 Da)
fig_mid = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_mid = massSpecPlot(pos, 0.01);
xlim([20 40]);
title('Mass Spectrum - Mid Mass Region (20-40 Da)');
for elem = {'Al', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Nb'}
    for cs = 1:3
        try
            ionAdd(spec_mid, [elem{1} '+'], cs, isoTable, colors, sumMargin, minAbundance);
        catch
        end
    end
end
saveas(fig_mid, fullfile(outputDir, 'mass_spectrum_20-40Da.png'));
close(fig_mid);

% High mass (40-100 Da)
fig_high = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
spec_high = massSpecPlot(pos, 0.01);
xlim([40 100]);
title('Mass Spectrum - High Mass Region (40-100 Da)');
for elem = {'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Nb', 'Mo'}
    for cs = 1:2
        try
            ionAdd(spec_high, [elem{1} '+'], cs, isoTable, colors, sumMargin, minAbundance);
        catch
        end
    end
end
saveas(fig_high, fullfile(outputDir, 'mass_spectrum_40-100Da.png'));
close(fig_high);

fprintf('Saved zoomed mass spectrum views\n');

%% Step 4: Apply ranges to pos
fprintf('\n=== Applying ranges to pos table ===\n');
posRanged = rangesFromPos(pos, ranges);
numRanged = sum(~ismissing(posRanged.ion));
fprintf('Ranged atoms: %d out of %d (%.1f%%)\n', ...
    numRanged, height(posRanged), 100 * numRanged / height(posRanged));

% Summary
if ismember('atom', posRanged.Properties.VariableNames)
    atomCounts = groupcounts(posRanged, 'atom');
    atomCounts = sortrows(atomCounts, 'GroupCount', 'descend');
    fprintf('\nAtom counts (top 15):\n');
    disp(atomCounts(1:min(15, height(atomCounts)), :));
end

%% Step 5: Create scatter plot
fprintf('\n=== Creating scatter plot visualization ===\n');

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

%% Step 6: Individual element plots
fprintf('\n=== Creating individual element scatter plots ===\n');

if ismember('atom', posRanged.Properties.VariableNames)
    uniqueAtoms = categories(posRanged.atom);
    uniqueAtoms = uniqueAtoms(~cellfun(@isempty, uniqueAtoms));
    elementsToPlot = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'H', 'O', 'C'};

    for i = 1:length(elementsToPlot)
        elem = elementsToPlot{i};
        if ismember(elem, uniqueAtoms)
            posElem = posRanged(posRanged.atom == elem, :);
            if height(posElem) > 0
                fig_elem = figure('Position', [100, 100, 800, 600], 'Visible', 'off');
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

%% Step 7: Save all data
fprintf('\n=== Saving all data ===\n');

save(fullfile(outputDir, 'pos.mat'), 'pos', '-v7.3');
save(fullfile(outputDir, 'posRanged.mat'), 'posRanged', '-v7.3');
fprintf('Saved pos.mat and posRanged.mat\n');

% Concentration summary
if ismember('atom', posRanged.Properties.VariableNames)
    atomCounts = groupcounts(posRanged, 'atom');
    totalRanged = sum(atomCounts.GroupCount);
    atomCounts.Concentration_pct = 100 * atomCounts.GroupCount / totalRanged;
    atomCounts = sortrows(atomCounts, 'GroupCount', 'descend');

    fprintf('\n=== Concentration Summary ===\n');
    for i = 1:min(15, height(atomCounts))
        fprintf('  %s: %.2f%% (%d atoms)\n', ...
            string(atomCounts.atom(i)), ...
            atomCounts.Concentration_pct(i), ...
            atomCounts.GroupCount(i));
    end

    writetable(atomCounts, fullfile(outputDir, 'concentration_summary.csv'));
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('All outputs saved to: %s\n', outputDir);

% List files
fprintf('\nSaved files:\n');
files = dir(outputDir);
for i = 1:length(files)
    if ~files(i).isdir
        fprintf('  - %s (%.1f KB)\n', files(i).name, files(i).bytes/1024);
    end
end

close all;
