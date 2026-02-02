function [excess, info] = interfacialExcessGUI(info)
% interfacialExcessGUI interactive GUI for interfacial excess calculation
%
% [excess, info] = interfacialExcessGUI(info)
%
% This function creates a programmatic uifigure GUI for interactively
% adjusting the fit limits for interfacial excess calculation.
%
% INPUT
% info:     struct from posCalculateInterfacialExcess containing:
%           - cumulative: cumulative curves (numAtoms x numSpecies)
%           - dist_sorted: sorted distance array
%           - sortIdx: sort indices
%           - interfaceLoc: interface location index
%           - area: interface area (or NaN)
%           - species: cell array of species names
%           - speciesFound: logical array
%           - mode: 'atomic', 'ionic', or 'isotopic'
%           - side: 'both', 'positive', or 'negative'
%           - limits: [lim1 lim2 lim3 lim4] initial limits
%
% OUTPUT
% excess:   table with calculated excess values
% info:     updated info struct with GUI handle
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Extract data from info
cumulative = info.cumulative;
numAtoms = size(cumulative, 1);
numSpecies = length(info.species);
limits = info.limits;
interfaceLoc = info.interfaceLoc;
area = info.area;
isSurface = strcmp(info.side, 'positive') || strcmp(info.side, 'negative');

%% Create figure
fig = uifigure('Name', 'Interfacial Excess Calculator', ...
    'Position', [100 100 800 600], ...
    'CloseRequestFcn', @closeFigure);

% Store data in figure UserData
fig.UserData.info = info;
fig.UserData.limits = limits;
fig.UserData.adjustMode = false;

%% Create layout
mainGrid = uigridlayout(fig, [3 1]);
mainGrid.RowHeight = {'1x', 180, 40};

%% Create axes panel
axPanel = uipanel(mainGrid);
axPanel.Layout.Row = 1;
ax = uiaxes(axPanel, 'Position', [50 50 700 350]);
ax.XLabel.String = 'Cumulative number of atoms';
ax.YLabel.String = 'Cumulative number of species atoms';
ax.Title.String = 'Cumulative Diagram';
hold(ax, 'on');
fig.UserData.ax = ax;

%% Define colors for species
colors = lines(numSpecies);

%% Plot cumulative curves
cumulativePlots = gobjects(numSpecies, 1);
fitPlotsLower = gobjects(numSpecies, 1);
fitPlotsUpper = gobjects(numSpecies, 1);

x = 1:numAtoms;

for s = 1:numSpecies
    if info.speciesFound(s)
        cumulativePlots(s) = plot(ax, x, cumulative(:, s), ...
            'LineWidth', 2, 'Color', colors(s, :), ...
            'DisplayName', info.species{s});

        % Placeholder for fit lines
        fitPlotsLower(s) = plot(ax, [1 numAtoms], [0 0], ...
            ':', 'LineWidth', 1, 'Color', colors(s, :), ...
            'HandleVisibility', 'off');
        fitPlotsUpper(s) = plot(ax, [1 numAtoms], [0 0], ...
            ':', 'LineWidth', 1, 'Color', colors(s, :), ...
            'HandleVisibility', 'off');
    end
end

fig.UserData.cumulativePlots = cumulativePlots;
fig.UserData.fitPlotsLower = fitPlotsLower;
fig.UserData.fitPlotsUpper = fitPlotsUpper;
fig.UserData.colors = colors;

%% Plot limit lines
yMax = max(cumulative(end, info.speciesFound));
limitLines = gobjects(4, 1);
limitColors = [0.3 0.3 0.3; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.3 0.3 0.3];

for i = 1:4
    limitLines(i) = xline(ax, limits(i), '-', 'Color', limitColors(i, :), ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end

fig.UserData.limitLines = limitLines;

%% Plot interface location
interfaceLine = xline(ax, interfaceLoc, '--r', 'LineWidth', 2, ...
    'Label', 'Interface', 'HandleVisibility', 'off');
fig.UserData.interfaceLine = interfaceLine;

%% Add legend
legend(ax, 'Location', 'northwest');

%% Create results panel
resultsPanel = uipanel(mainGrid, 'Title', 'Results');
resultsPanel.Layout.Row = 2;

% Create table
columnNames = {'Species', 'Excess (at)', 'Excess (at/nm^2)', 'Partial 1', 'Partial 2'};
columnWidth = {100, 100, 120, 100, 100};

resultsTable = uitable(resultsPanel, ...
    'Position', [10 10 760 140], ...
    'ColumnName', columnNames, ...
    'ColumnWidth', columnWidth, ...
    'RowName', {});

fig.UserData.resultsTable = resultsTable;

%% Create button panel
buttonPanel = uipanel(mainGrid);
buttonPanel.Layout.Row = 3;

buttonGrid = uigridlayout(buttonPanel, [1 5]);
buttonGrid.ColumnWidth = {'1x', '1x', '1x', '1x', '1x'};

adjustBtn = uibutton(buttonGrid, 'Text', 'Adjust Limits', ...
    'ButtonPushedFcn', @adjustLimitsCallback);
adjustBtn.Layout.Column = 1;

resetBtn = uibutton(buttonGrid, 'Text', 'Reset Limits', ...
    'ButtonPushedFcn', @resetLimitsCallback);
resetBtn.Layout.Column = 2;

saveBtn = uibutton(buttonGrid, 'Text', 'Save Figure', ...
    'ButtonPushedFcn', @saveFigureCallback);
saveBtn.Layout.Column = 3;

exportBtn = uibutton(buttonGrid, 'Text', 'Export Table', ...
    'ButtonPushedFcn', @exportTableCallback);
exportBtn.Layout.Column = 4;

doneBtn = uibutton(buttonGrid, 'Text', 'Done', ...
    'ButtonPushedFcn', @doneCallback);
doneBtn.Layout.Column = 5;

fig.UserData.adjustBtn = adjustBtn;

%% Initial calculation
recalculateExcess(fig);

%% Store handle in info
info.guiHandle = fig;

%% Wait for user to close or click Done
uiwait(fig);

%% Return results
if isvalid(fig)
    excess = fig.UserData.excess;
    info = fig.UserData.info;
    info.limits = fig.UserData.limits;
    delete(fig);
else
    % Figure was closed
    excess = calculateExcessFromLimits(info);
    info.limits = limits;
end

end

%% ========== Callback Functions ==========

function adjustLimitsCallback(~, ~)
    fig = gcbf;
    fig.UserData.adjustMode = true;
    fig.UserData.adjustBtn.Text = 'Click twice on plot...';
    fig.UserData.adjustBtn.Enable = 'off';

    % Get two clicks
    ax = fig.UserData.ax;
    [xClicks, ~] = ginput(2);

    if length(xClicks) >= 2
        % Find nearest limit to first click
        limits = fig.UserData.limits;
        dist = abs(limits - xClicks(1));
        [~, targetIdx] = min(dist);

        % Move that limit to second click position
        newLimit = round(xClicks(2));
        newLimit = max(1, min(size(fig.UserData.info.cumulative, 1), newLimit));

        limits(targetIdx) = newLimit;

        % Ensure limits are in order
        limits = sort(limits);

        fig.UserData.limits = limits;

        % Update plot
        recalculateExcess(fig);
    end

    fig.UserData.adjustMode = false;
    fig.UserData.adjustBtn.Text = 'Adjust Limits';
    fig.UserData.adjustBtn.Enable = 'on';
end

function resetLimitsCallback(~, ~)
    fig = gcbf;
    info = fig.UserData.info;

    % Recalculate initial limits
    numAtoms = size(info.cumulative, 1);
    interfaceLoc = info.interfaceLoc;
    isSurface = strcmp(info.side, 'positive') || strcmp(info.side, 'negative');

    delta = round(numAtoms / 10);

    if isSurface
        if strcmp(info.side, 'positive')
            limits = [1, 1, interfaceLoc + delta, numAtoms];
        else
            limits = [1, interfaceLoc - delta, numAtoms, numAtoms];
        end
    else
        limL = max(1, interfaceLoc - delta);
        limU = min(numAtoms, interfaceLoc + delta);
        limits = [1, limL, limU, numAtoms];
    end

    limits = round(limits);
    limits = max(1, min(numAtoms, limits));

    fig.UserData.limits = limits;
    recalculateExcess(fig);
end

function saveFigureCallback(~, ~)
    fig = gcbf;
    ax = fig.UserData.ax;

    % Create new figure with plot
    newFig = figure('Name', 'Interfacial Excess - Cumulative Diagram', ...
        'Color', 'white');
    newAx = copyobj(ax, newFig);
    newAx.Units = 'normalized';
    newAx.Position = [0.1 0.15 0.85 0.75];

    % Add text with results
    excess = fig.UserData.excess;
    textStr = '';
    for s = 1:height(excess)
        if ~isnan(excess.excess_atoms(s))
            textStr = [textStr sprintf('%s: %.1f at', excess.species{s}, excess.excess_atoms(s))];
            if ~isnan(excess.excess_per_nm2(s))
                textStr = [textStr sprintf(' (%.3f at/nm^2)', excess.excess_per_nm2(s))];
            end
            textStr = [textStr newline];
        end
    end

    annotation(newFig, 'textbox', [0.15 0.02 0.7 0.1], ...
        'String', textStr, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', 'FontSize', 10);

    % Prompt for save
    [file, path] = uiputfile({'*.fig'; '*.png'; '*.pdf'}, 'Save Figure');
    if file
        saveas(newFig, fullfile(path, file));
        delete(newFig);
    end
end

function exportTableCallback(~, ~)
    fig = gcbf;
    excess = fig.UserData.excess;

    [file, path] = uiputfile({'*.csv'; '*.xlsx'}, 'Export Table');
    if file
        writetable(excess, fullfile(path, file));
    end
end

function doneCallback(~, ~)
    fig = gcbf;
    uiresume(fig);
end

function closeFigure(src, ~)
    uiresume(src);
    delete(src);
end

%% ========== Helper Functions ==========

function recalculateExcess(fig)
% Recalculate excess values and update display

    info = fig.UserData.info;
    limits = fig.UserData.limits;
    cumulative = info.cumulative;
    numAtoms = size(cumulative, 1);
    numSpecies = length(info.species);
    interfaceLoc = info.interfaceLoc;
    area = info.area;
    isSurface = strcmp(info.side, 'positive') || strcmp(info.side, 'negative');

    x = (1:numAtoms)';

    % Update limit lines
    limitLines = fig.UserData.limitLines;
    for i = 1:4
        limitLines(i).Value = limits(i);
    end

    % Calculate excess for each species
    speciesNames = info.species(:);
    excess_atoms = zeros(numSpecies, 1);
    excess_per_nm2 = zeros(numSpecies, 1);
    partial1 = zeros(numSpecies, 1);
    partial2 = zeros(numSpecies, 1);

    fitPlotsLower = fig.UserData.fitPlotsLower;
    fitPlotsUpper = fig.UserData.fitPlotsUpper;

    for s = 1:numSpecies
        if ~info.speciesFound(s)
            excess_atoms(s) = NaN;
            excess_per_nm2(s) = NaN;
            partial1(s) = NaN;
            partial2(s) = NaN;
            continue;
        end

        cum = cumulative(:, s);

        % Linear regression
        if isSurface
            if strcmp(info.side, 'positive')
                a1 = 0; b1 = 0;
                fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
                a2 = fit2(1); b2 = fit2(2);
            else
                fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
                a1 = fit1(1); b1 = fit1(2);
                a2 = 0; b2 = cum(end);
            end
        else
            fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
            a1 = fit1(1); b1 = fit1(2);
            fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
            a2 = fit2(1); b2 = fit2(2);
        end

        % Update fit line plots
        fitPlotsLower(s).XData = [1 numAtoms];
        fitPlotsLower(s).YData = [a1*1+b1, a1*numAtoms+b1];
        fitPlotsUpper(s).XData = [1 numAtoms];
        fitPlotsUpper(s).YData = [a2*1+b2, a2*numAtoms+b2];

        % Calculate intersections at interface
        intersect1 = a1 * interfaceLoc + b1;
        intersect2 = a2 * interfaceLoc + b2;

        % Excess
        excess_atoms(s) = intersect2 - intersect1;
        if ~isnan(area) && area > 0
            excess_per_nm2(s) = excess_atoms(s) / area;
        else
            excess_per_nm2(s) = NaN;
        end

        % Partial excesses
        cumAtInterface = cum(round(interfaceLoc));
        partial1(s) = cumAtInterface - intersect1;
        partial2(s) = intersect2 - cumAtInterface;

        if ~isnan(area) && area > 0
            partial1(s) = partial1(s) / area;
            partial2(s) = partial2(s) / area;
        end
    end

    % Build output table
    fit_limits = repmat(limits, numSpecies, 1);
    excess = table(speciesNames, excess_atoms, excess_per_nm2, partial1, partial2, fit_limits, ...
        'VariableNames', {'species', 'excess_atoms', 'excess_per_nm2', 'partial1', 'partial2', 'fit_limits'});

    fig.UserData.excess = excess;

    % Update results table display
    tableData = cell(numSpecies, 5);
    for s = 1:numSpecies
        tableData{s, 1} = info.species{s};
        if isnan(excess_atoms(s))
            tableData{s, 2} = 'N/A';
            tableData{s, 3} = 'N/A';
            tableData{s, 4} = 'N/A';
            tableData{s, 5} = 'N/A';
        else
            tableData{s, 2} = sprintf('%.1f', excess_atoms(s));
            if isnan(excess_per_nm2(s))
                tableData{s, 3} = 'N/A';
                tableData{s, 4} = sprintf('%.1f at', partial1(s) * area);
                tableData{s, 5} = sprintf('%.1f at', partial2(s) * area);
            else
                tableData{s, 3} = sprintf('%.4f', excess_per_nm2(s));
                tableData{s, 4} = sprintf('%.4f', partial1(s));
                tableData{s, 5} = sprintf('%.4f', partial2(s));
            end
        end
    end

    fig.UserData.resultsTable.Data = tableData;
end

function excess = calculateExcessFromLimits(info)
% Calculate excess when GUI is closed without clicking Done

    limits = info.limits;
    cumulative = info.cumulative;
    numAtoms = size(cumulative, 1);
    numSpecies = length(info.species);
    interfaceLoc = info.interfaceLoc;
    area = info.area;
    isSurface = strcmp(info.side, 'positive') || strcmp(info.side, 'negative');

    x = (1:numAtoms)';

    speciesNames = info.species(:);
    excess_atoms = zeros(numSpecies, 1);
    excess_per_nm2 = zeros(numSpecies, 1);
    partial1 = zeros(numSpecies, 1);
    partial2 = zeros(numSpecies, 1);

    for s = 1:numSpecies
        if ~info.speciesFound(s)
            excess_atoms(s) = NaN;
            excess_per_nm2(s) = NaN;
            partial1(s) = NaN;
            partial2(s) = NaN;
            continue;
        end

        cum = cumulative(:, s);

        if isSurface
            if strcmp(info.side, 'positive')
                a1 = 0; b1 = 0;
                fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
                a2 = fit2(1); b2 = fit2(2);
            else
                fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
                a1 = fit1(1); b1 = fit1(2);
                a2 = 0; b2 = cum(end);
            end
        else
            fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
            a1 = fit1(1); b1 = fit1(2);
            fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
            a2 = fit2(1); b2 = fit2(2);
        end

        intersect1 = a1 * interfaceLoc + b1;
        intersect2 = a2 * interfaceLoc + b2;

        excess_atoms(s) = intersect2 - intersect1;
        if ~isnan(area) && area > 0
            excess_per_nm2(s) = excess_atoms(s) / area;
        else
            excess_per_nm2(s) = NaN;
        end

        cumAtInterface = cum(round(interfaceLoc));
        partial1(s) = cumAtInterface - intersect1;
        partial2(s) = intersect2 - cumAtInterface;

        if ~isnan(area) && area > 0
            partial1(s) = partial1(s) / area;
            partial2(s) = partial2(s) / area;
        end
    end

    fit_limits = repmat(limits, numSpecies, 1);
    excess = table(speciesNames, excess_atoms, excess_per_nm2, partial1, partial2, fit_limits, ...
        'VariableNames', {'species', 'excess_atoms', 'excess_per_nm2', 'partial1', 'partial2', 'fit_limits'});
end
