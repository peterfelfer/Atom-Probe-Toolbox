%% Automated Analysis Script for A718 EPOS data
% Uses automatic ranging based on detected peaks

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

%% Step 1: Find peaks using new massSpecFindPeaks
fprintf('\n=== Step 1: Finding peaks ===\n');
[peakTable, peakInfo] = massSpecFindPeaks(pos, 'binWidth', 0.02, 'mcRange', [0 100]);
fprintf('Found %d peaks\n', height(peakTable));

%% Step 2: Auto-assign ions using massSpecAutoAssignIons
fprintf('\n=== Step 2: Auto-assigning ions ===\n');
elements = {'Ni', 'Cr', 'Fe', 'Nb', 'Mo', 'Ti', 'Al', 'Co', 'C', 'Mn', 'Si', 'H', 'O', 'N'};

[ionTable, assignInfo] = massSpecAutoAssignIons(peakTable, elements, ...
    'maxAtoms', 2, ...
    'maxCharge', 3, ...
    'conservative', true);

fprintf('Assigned %d ions\n', height(ionTable));
if height(ionTable) > 0
    disp(ionTable);
end

%% Step 3: Generate ranges using massSpecAutoRanges
fprintf('\n=== Step 3: Generating ranges ===\n');
if ~isempty(assignInfo.assigned)
    [rangeTableAuto, rangeInfo] = massSpecAutoRanges(pos, peakTable, assignInfo.assigned, ...
        'binWidth', 0.02);
    fprintf('Generated %d ranges\n', height(rangeTableAuto));
else
    fprintf('No ions assigned.\n');
    rangeTableAuto = table();
end

%% Step 4: Convert auto-ranges to standard range format for posAllocateRange
fprintf('\n=== Step 4: Converting ranges to standard format ===\n');

% Load isotope table for ionConvertName
data = load('isotopeTable_naturalAbundances.mat');
isoTable = data.isotopeTable;

if height(rangeTableAuto) > 0
    numRanges = height(rangeTableAuto);
    mcbegin = rangeTableAuto.mcMin;
    mcend = rangeTableAuto.mcMax;
    chargeState = rangeTableAuto.chargeState;

    % Create rangeName from ion names
    rangeName = categorical(string(rangeTableAuto.ionName));

    % Create ion cell array (element/isotope tables) using ionConvertName
    ionCell = cell(numRanges, 1);
    for i = 1:numRanges
        ionStr = rangeTableAuto.ion{i};
        try
            ionCell{i} = ionConvertName(ionStr);
        catch
            % Fallback: create simple table
            ionCell{i} = table(categorical({ionStr}), 0, 'VariableNames', {'element', 'isotope'});
        end
    end

    % Create color (dummy values)
    color = repmat([0.5 0.5 0.5], numRanges, 1);

    % Build ranges table in format expected by posAllocateRange
    ranges = table(mcbegin, mcend, rangeName, chargeState, ionCell, color, ...
        'VariableNames', {'mcbegin', 'mcend', 'rangeName', 'chargeState', 'ion', 'color'});

    fprintf('Created %d standard ranges\n', height(ranges));
    disp(ranges(1:min(20, height(ranges)), {'mcbegin', 'mcend', 'rangeName', 'chargeState'}));
else
    ranges = table();
    fprintf('No ranges created.\n');
end

%% Step 5: Save ranges
writetable(ranges, fullfile(outputDir, 'ranges.csv'));
save(fullfile(outputDir, 'ranges.mat'), 'ranges');
fprintf('Saved ranges\n');

%% Step 6: Create mass spectrum figure
fprintf('\n=== Step 5: Creating mass spectrum figure ===\n');
fig1 = figure('Position', [100, 100, 1600, 800], 'Visible', 'off');

% Create histogram
binWidth = 0.01;
edges = 0:binWidth:100;
counts = histcounts(pos.mc, edges);
centers = edges(1:end-1) + binWidth/2;

area(centers, counts, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
hold on;
set(gca, 'YScale', 'log');
xlabel('mass-to-charge-state ratio [Da]');
ylabel('counts');
title('Mass Spectrum - A718 (Inconel 718)');

% Add range markers
if height(ranges) > 0
    yLim = ylim;
    uniqueIons = unique(ranges.rangeName);
    colorMap = lines(numel(uniqueIons));

    for i = 1:height(ranges)
        ionIdx = find(uniqueIons == ranges.rangeName(i));
        rangeColor = colorMap(ionIdx, :);

        % Draw range box
        x = [ranges.mcbegin(i), ranges.mcend(i), ranges.mcend(i), ranges.mcbegin(i)];
        y = [yLim(1), yLim(1), yLim(2), yLim(2)];
        patch(x, y, rangeColor, 'FaceAlpha', 0.2, 'EdgeColor', rangeColor, 'LineWidth', 0.5);
    end
end

hold off;
xlim([0 100]);

saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ranges.png'));
saveas(fig1, fullfile(outputDir, 'mass_spectrum_with_ranges.fig'));
fprintf('Saved mass spectrum with ranges\n');

%% Create zoomed views
% Low mass
fig_low = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
area(centers, counts, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
xlim([0 20]);
xlabel('mass-to-charge-state ratio [Da]');
ylabel('counts');
title('Mass Spectrum - Low Mass Region (0-20 Da)');
saveas(fig_low, fullfile(outputDir, 'mass_spectrum_0-20Da.png'));
close(fig_low);

% Mid mass
fig_mid = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
area(centers, counts, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
xlim([20 40]);
xlabel('mass-to-charge-state ratio [Da]');
ylabel('counts');
title('Mass Spectrum - Mid Mass Region (20-40 Da)');
saveas(fig_mid, fullfile(outputDir, 'mass_spectrum_20-40Da.png'));
close(fig_mid);

% High mass
fig_high = figure('Position', [100, 100, 1200, 500], 'Visible', 'off');
area(centers, counts, 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
xlim([40 100]);
xlabel('mass-to-charge-state ratio [Da]');
ylabel('counts');
title('Mass Spectrum - High Mass Region (40-100 Da)');
saveas(fig_high, fullfile(outputDir, 'mass_spectrum_40-100Da.png'));
close(fig_high);

fprintf('Saved zoomed mass spectrum views\n');

%% Step 7: Apply ranges to pos using posAllocateRange
fprintf('\n=== Step 6: Applying ranges to pos table ===\n');

if height(ranges) > 0
    try
        posRanged = posAllocateRange(pos, ranges, 'decompose');
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
    catch ME
        fprintf('Warning: Could not apply ranges: %s\n', ME.message);
        posRanged = pos;
    end
else
    posRanged = pos;
    fprintf('No ranges to apply.\n');
end

%% Step 8: Create scatter plot
fprintf('\n=== Step 7: Creating scatter plot visualization ===\n');

fig3 = figure('Position', [100, 100, 1200, 900], 'Visible', 'off');

if ismember('atom', posRanged.Properties.VariableNames) && sum(~ismissing(posRanged.atom)) > 0
    % Use scatterPlotPosData if atoms are ranged
    scatterPlotPosData(posRanged);
else
    % Basic scatter plot
    if height(pos) > 500000
        idx = randperm(height(pos), 500000);
        scatter3(pos.x(idx), pos.y(idx), pos.z(idx), 1, pos.mc(idx), '.');
    else
        scatter3(pos.x, pos.y, pos.z, 1, pos.mc, '.');
    end
    colorbar;
    title('Colored by mass-to-charge');
end

title('3D Atom Map - A718 (Inconel 718)');
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
view(3);
axis equal;

saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.png'));
saveas(fig3, fullfile(outputDir, 'scatter_plot_all_ions.fig'));
fprintf('Saved scatter plot of all ions\n');

%% Step 9: Individual element plots
fprintf('\n=== Step 8: Creating individual element scatter plots ===\n');

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

%% Step 10: Save all data
fprintf('\n=== Step 9: Saving all data ===\n');

save(fullfile(outputDir, 'pos.mat'), 'pos', '-v7.3');
save(fullfile(outputDir, 'posRanged.mat'), 'posRanged', '-v7.3');
save(fullfile(outputDir, 'peakAnalysis.mat'), 'peakTable', 'peakInfo');
save(fullfile(outputDir, 'ionAssignment.mat'), 'ionTable', 'assignInfo');
if exist('rangeTableAuto', 'var')
    save(fullfile(outputDir, 'rangesAuto.mat'), 'rangeTableAuto', 'rangeInfo');
end
fprintf('Saved all MAT files\n');

% Concentration summary
if ismember('atom', posRanged.Properties.VariableNames) && sum(~ismissing(posRanged.atom)) > 0
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
