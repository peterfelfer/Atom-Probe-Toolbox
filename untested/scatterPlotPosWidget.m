function [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme, options)
% SCATTERPLOTPOSWIDGET Interactive scatter plot with visibility/sample controls.
%
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme)
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme, 'groupBy', 'ion')
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, [], 'axes', ax)
%
% Creates a 3D scatter plot of APT data and a separate control window that
% toggles species visibility and adjusts the displayed sampling fraction.
%
% INPUT:
%   pos         table with x, y, z and a species column (ion/atom/isotope)
%   colorScheme table with columns: ion (categorical/string) and color (Nx3)
%
% OPTIONS:
%   'axes'         - Target axes for scatter plot (default: new figure)
%   'groupBy'      - 'ion' | 'atom' | 'auto' (default: 'auto')
%   'sample'       - Initial sample fraction (0..1) or count (>1) (default: 0.1)
%   'markerSize'   - Marker size for scatter (default: 15)
%   'markerSizeGlobal' - Global marker size (default: markerSize)
%   'randomSeed'   - Random seed for reproducible sampling (default: 1)
%   'showUnranged' - Include unranged/missing category (default: true)
%   'controlTitle' - Title for control window (default: 'APT Scatter Controls')
%   'splitIsotope' - Split by isotope if available (default: false)
%   'splitCharge'  - Split by charge state if available (default: false)
%   'fixBoundingBox' - Fix axis limits to dataset bounds (default: true)
%   'clipToBoundingBox' - Clip points to bounding box (default: true)
%   'bboxUseSpan'  - Use center/span edits for bounding box (default: false)
%   'bboxPadding'  - Fractional padding for bounding box (default: 0)
%
% OUTPUT:
%   scatterHandles - array of scatter handles (one per species)
%   ax             - axes handle
%   controlFig     - control window figure handle
%   info           - struct with fields: speciesNames, speciesCounts, colors
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    colorScheme = []
    options.axes = []
    options.groupBy (1,1) string = "auto"
    options.sample (1,1) double = NaN
    options.markerSize (1,1) double {mustBePositive} = 15
    options.markerSizeGlobal (1,1) double = NaN
    options.randomSeed (1,1) double = 1
    options.showUnranged (1,1) logical = true
    options.controlTitle (1,1) string = "APT Scatter Controls"
    options.splitIsotope (1,1) logical = false
    options.splitCharge (1,1) logical = false
    options.fixBoundingBox (1,1) logical = true
    options.clipToBoundingBox (1,1) logical = true
    options.bboxUseSpan (1,1) logical = false
    options.bboxPadding (1,1) double {mustBeNonnegative} = 0
end

mode = resolveMode(pos, options.groupBy);
coords = [pos.x, pos.y, pos.z];

[speciesNames, baseNames, displayNames, speciesIndices, speciesCounts] = computeGroups( ...
    pos, mode, options.splitIsotope, options.splitCharge, options.showUnranged);
colors = mapColors(baseNames, colorScheme);

if isempty(options.axes)
    fig = figure('Name', 'APT Scatter', 'Color', 'w');
    ax = axes(fig);
else
    ax = options.axes;
    fig = ancestor(ax, 'figure');
end

markerSizeGlobal = resolveMarkerSizeGlobal(options.markerSizeGlobal, options.markerSize);
markerSizes = initMarkerSizes(numel(speciesNames), markerSizeGlobal);

axes(ax);
holdState = ishold(ax);
hold(ax, 'on');
scatterHandles = createScatterHandles(ax, displayNames, colors, markerSizes);
if ~holdState
    hold(ax, 'off');
end

axisSpatialAptify;

controlFig = createControlWindow(options.controlTitle);

data = struct();
data.pos = pos;
data.colorScheme = colorScheme;
data.coordsOriginal = coords;
data.coords = coords;
data.speciesNames = speciesNames;
data.baseNames = baseNames;
data.displayNames = displayNames;
data.speciesIndices = speciesIndices;
data.speciesCounts = speciesCounts;
data.colors = colors;
data.visible = true(numel(speciesNames), 1);
data.sampleValue = resolveInitialSample(height(pos), options.sample);
data.sampleFractions = initSampleFractions(speciesCounts, data.sampleValue);
data.scatterHandles = scatterHandles;
data.ax = ax;
data.controlFig = controlFig;
data.markerSize = markerSizeGlobal;
data.markerSizeGlobal = markerSizeGlobal;
data.markerSizes = markerSizes;
data.randomSeed = options.randomSeed;
data.rotationAngle = 0;
data.showUnranged = options.showUnranged;
data.mode = mode;
data.splitIsotope = options.splitIsotope;
data.splitCharge = options.splitCharge;
data.hasIsotope = ismember('isotope', pos.Properties.VariableNames);
data.hasCharge = ismember('chargeState', pos.Properties.VariableNames);
data.bboxOriginal = computeBoundingBox(coords, options.bboxPadding);
data.bbox = data.bboxOriginal;
data.fixBoundingBox = options.fixBoundingBox;
data.clipToBoundingBox = options.clipToBoundingBox;
data.bboxUseSpan = options.bboxUseSpan;
data.bboxPadding = options.bboxPadding;
data.bboxUserEdited = false;
data.fn = struct( ...
    'computeGroups', @computeGroups, ...
    'mapColors', @mapColors, ...
    'createScatterHandles', @createScatterHandles, ...
    'updateScatter', @updateScatter, ...
    'updateTable', @updateTable, ...
    'buildTableData', @buildTableData, ...
    'initSampleFractions', @initSampleFractions, ...
    'initMarkerSizes', @initMarkerSizes, ...
    'rotateCoordsZ', @rotateCoordsZ);

tableData = buildTableData(data);

setupControls(controlFig, tableData, data);
updateScatter(controlFig);
applyBoundingBox(ax, data.bbox, data.fixBoundingBox);

info = struct();
info.speciesNames = speciesNames;
info.speciesCounts = speciesCounts;
info.colors = colors;
info.groupBy = mode;
end

function mode = resolveMode(pos, groupBy)
    if groupBy == "auto"
        if ismember('ion', pos.Properties.VariableNames)
            mode = "ionic";
        elseif ismember('atom', pos.Properties.VariableNames)
            mode = "atomic";
        else
            error('scatterPlotPosWidget:missingGroupField', ...
                'pos must contain ion or atom column.');
        end
    else
        if groupBy == "ion"
            mode = "ionic";
        elseif groupBy == "atom"
            mode = "atomic";
        else
            mode = "ionic";
        end
    end
end

function [groupNames, baseNames, displayNames, groupIndices, groupCounts] = computeGroups( ...
    pos, mode, splitIsotope, splitCharge, showUnranged)
    if mode == "ionic"
        if ismember('ion', pos.Properties.VariableNames)
            base = string(pos.ion);
        elseif ismember('atom', pos.Properties.VariableNames)
            base = string(pos.atom);
        else
            error('scatterPlotPosWidget:missingGroupField', ...
                'pos must contain ion or atom column.');
        end
    else
        if ismember('atom', pos.Properties.VariableNames)
            base = string(pos.atom);
        elseif ismember('ion', pos.Properties.VariableNames)
            base = string(pos.ion);
        else
            error('scatterPlotPosWidget:missingGroupField', ...
                'pos must contain ion or atom column.');
        end
    end

    if showUnranged
        base(ismissing(base)) = "unranged";
    else
        keep = ~ismissing(base);
        base = base(keep);
    end

    groupKey = base;
    if splitIsotope && ismember('isotope', pos.Properties.VariableNames)
        iso = string(pos.isotope);
        iso(ismissing(iso)) = "NaN";
        groupKey = groupKey + "-" + iso;
    end
    if splitCharge && ismember('chargeState', pos.Properties.VariableNames)
        cs = pos.chargeState;
        csString = strings(size(cs));
        csString(:) = "";
        valid = ~isnan(cs);
        for k = 1:numel(cs)
            if valid(k)
                n = round(cs(k));
                if n <= 0
                    csString(k) = "";
                else
                    csString(k) = repmat("+", 1, n);
                end
            end
        end
        groupKey = groupKey + csString;
    end

    [groupNames, ~, groupIdx] = unique(groupKey, 'stable');
    baseNames = strings(size(groupNames));
    displayNames = strings(size(groupNames));
    groupIndices = cell(numel(groupNames), 1);
    groupCounts = zeros(numel(groupNames), 1);

    for i = 1:numel(groupNames)
        idx = find(groupIdx == i);
        groupIndices{i} = idx;
        groupCounts(i) = numel(idx);
        baseNames(i) = base(idx(1));
        displayNames(i) = formatDisplayName(base(idx(1)), pos, idx(1), splitIsotope, splitCharge, showUnranged);
    end
end

function colors = mapColors(baseNames, colorScheme)
    n = numel(baseNames);
    defaultColors = lines(max(n, 7));
    colors = zeros(n, 3);

    hasScheme = istable(colorScheme) && ...
        all(ismember({'ion', 'color'}, colorScheme.Properties.VariableNames));
    if hasScheme
        ionNames = string(colorScheme.ion);
    else
        ionNames = string.empty(0, 1);
    end

    for i = 1:n
        name = baseNames(i);
        if hasScheme
            idx = find(ionNames == name, 1, 'first');
            if isempty(idx) && name == "unranged"
                idx = find(ionNames == "background", 1, 'first');
            end
            if ~isempty(idx)
                colors(i, :) = colorScheme.color(idx, :);
                continue;
            end
        end
        colors(i, :) = defaultColors(mod(i-1, size(defaultColors, 1)) + 1, :);
    end
end

function controlFig = createControlWindow(titleText)
    controlFig = figure('Name', char(titleText), ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'none', ...
        'Color', 'w', ...
        'Units', 'pixels', ...
        'Position', [100 100 360 680]);
end

function setupControls(controlFig, tableData, data)
    setappdata(controlFig, 'scatterPlotPosWidget', data);

    modeGroup = uibuttongroup(controlFig, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.90 0.9 0.08], ...
        'SelectionChangedFcn', @(src, evd) onGroupingChanged(controlFig, evd, src));

    rbIonic = uicontrol(modeGroup, 'Style', 'radiobutton', ...
        'String', 'Ionic', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.1 0.4 0.8]);

    rbAtomic = uicontrol(modeGroup, 'Style', 'radiobutton', ...
        'String', 'Atomic', ...
        'Units', 'normalized', ...
        'Position', [0.5 0.1 0.4 0.8]);

    if data.mode == "atomic"
        modeGroup.SelectedObject = rbAtomic;
    else
        modeGroup.SelectedObject = rbIonic;
    end

    splitIsotope = uicontrol(controlFig, 'Style', 'checkbox', ...
        'String', 'Split isotopes', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.85 0.45 0.04], ...
        'Value', data.splitIsotope, ...
        'Enable', boolToOnOff(data.hasIsotope), ...
        'Callback', @(~, ~) onGroupingChanged(controlFig));

    splitCharge = uicontrol(controlFig, 'Style', 'checkbox', ...
        'String', 'Split charge state', ...
        'Units', 'normalized', ...
        'Position', [0.52 0.85 0.43 0.04], ...
        'Value', data.splitCharge, ...
        'Enable', boolToOnOff(data.hasCharge), ...
        'Callback', @(~, ~) onGroupingChanged(controlFig));

    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Sample fraction', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.79 0.4 0.05], ...
        'HorizontalAlignment', 'left');

    slider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', min(data.sampleValue, 1), ...
        'Units', 'normalized', ...
        'Position', [0.05 0.75 0.65 0.04], ...
        'Callback', @(src, ~) onSlider(src, controlFig));

    editBox = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.sampleValue, '%.3f'), ...
        'Units', 'normalized', ...
        'Position', [0.72 0.75 0.2 0.05], ...
        'Callback', @(src, ~) onEdit(src, controlFig));

    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Marker size', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.69 0.4 0.05], ...
        'HorizontalAlignment', 'left');

    markerMax = max(50, data.markerSizeGlobal);
    markerSlider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 1, 'Max', markerMax, 'Value', min(data.markerSizeGlobal, markerMax), ...
        'Units', 'normalized', ...
        'Position', [0.05 0.65 0.65 0.04], ...
        'Callback', @(src, ~) onMarkerSizeSlider(src, controlFig));

    markerEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.markerSizeGlobal, '%.2f'), ...
        'Units', 'normalized', ...
        'Position', [0.72 0.65 0.2 0.05], ...
        'Callback', @(src, ~) onMarkerSizeEdit(src, controlFig));

    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Random seed', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.60 0.4 0.05], ...
        'HorizontalAlignment', 'left');

    seedEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.randomSeed, '%.0f'), ...
        'Units', 'normalized', ...
        'Position', [0.72 0.60 0.2 0.05], ...
        'Callback', @(src, ~) onRandomSeedEdit(src, controlFig));

    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Rotation (deg)', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.55 0.4 0.05], ...
        'HorizontalAlignment', 'left');

    rotationSlider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 0, 'Max', 360, 'Value', data.rotationAngle, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.51 0.65 0.04], ...
        'Callback', @(src, ~) onRotationSlider(src, controlFig));

    rotationEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.rotationAngle, '%.1f'), ...
        'Units', 'normalized', ...
        'Position', [0.72 0.51 0.2 0.05], ...
        'Callback', @(src, ~) onRotationEdit(src, controlFig));

    showAllBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'Show All', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.47 0.42 0.05], ...
        'Callback', @(~, ~) onShowAll(controlFig));

    hideAllBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'Hide All', ...
        'Units', 'normalized', ...
        'Position', [0.53 0.47 0.42 0.05], ...
        'Callback', @(~, ~) onHideAll(controlFig));

    bboxPanel = uipanel(controlFig, 'Title', 'Bounding Box', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.30 0.9 0.16]);

    bboxFix = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Fix', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.72 0.12 0.22], ...
        'Value', data.fixBoundingBox, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'toggle'));

    bboxClip = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Clip', ...
        'Units', 'normalized', ...
        'Position', [0.14 0.72 0.12 0.22], ...
        'Value', data.clipToBoundingBox, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'clip'));

    bboxSpan = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Span', ...
        'Units', 'normalized', ...
        'Position', [0.26 0.72 0.16 0.22], ...
        'Value', data.bboxUseSpan, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'span'));

    bboxUseCurrent = uicontrol(bboxPanel, 'Style', 'pushbutton', ...
        'String', 'Use Current', ...
        'Units', 'normalized', ...
        'Position', [0.42 0.72 0.24 0.22], ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'current'));

    bboxReset = uicontrol(bboxPanel, 'Style', 'pushbutton', ...
        'String', 'Reset', ...
        'Units', 'normalized', ...
        'Position', [0.68 0.72 0.30 0.22], ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'reset'));

    bboxEdits = struct();
    [bboxEdits.xmin, bboxEdits.xmax] = makeBBoxRow(bboxPanel, 'X', 0.44, data.bbox.xmin, data.bbox.xmax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));
    [bboxEdits.ymin, bboxEdits.ymax] = makeBBoxRow(bboxPanel, 'Y', 0.22, data.bbox.ymin, data.bbox.ymax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));
    [bboxEdits.zmin, bboxEdits.zmax] = makeBBoxRow(bboxPanel, 'Z', 0.00, data.bbox.zmin, data.bbox.zmax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));

    exportPanel = uipanel(controlFig, 'Title', 'Export', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.22 0.9 0.07]);

    exportImagesBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Export Images', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.15 0.42 0.7], ...
        'Callback', @(~, ~) onExportImages(controlFig));

    exportTurntableBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Turntable', ...
        'Units', 'normalized', ...
        'Position', [0.53 0.15 0.42 0.7], ...
        'Callback', @(~, ~) onExportTurntable(controlFig));

    tbl = uitable(controlFig, ...
        'Data', tableData, ...
        'ColumnName', {'Show', 'Species', 'Count', 'Number Shown', 'Fraction Shown', 'Marker Size'}, ...
        'ColumnEditable', [true false false true true true], ...
        'Units', 'normalized', ...
        'Position', [0.05 0.05 0.9 0.16], ...
        'CellEditCallback', @(src, evd) onTableEdit(src, evd, controlFig));

    setappdata(controlFig, 'scatterPlotPosWidgetControls', struct( ...
        'slider', slider, 'editBox', editBox, 'table', tbl, ...
        'markerSlider', markerSlider, 'markerEdit', markerEdit, ...
        'seedEdit', seedEdit, ...
        'rotationSlider', rotationSlider, 'rotationEdit', rotationEdit, ...
        'exportPanel', exportPanel, ...
        'exportImagesBtn', exportImagesBtn, 'exportTurntableBtn', exportTurntableBtn, ...
        'modeGroup', modeGroup, 'splitIsotope', splitIsotope, 'splitCharge', splitCharge, ...
        'bboxPanel', bboxPanel, 'bboxFix', bboxFix, 'bboxClip', bboxClip, 'bboxSpan', bboxSpan, 'bboxUseCurrent', bboxUseCurrent, ...
        'bboxReset', bboxReset, 'bboxEdits', bboxEdits));

    drawnow;
end

function updateScatter(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    if ~isfield(data, 'markerSizes') || numel(data.markerSizes) ~= numel(data.speciesNames)
        data.markerSizes = initMarkerSizes(numel(data.speciesNames), data.markerSizeGlobal);
        setappdata(controlFig, 'scatterPlotPosWidget', data);
    end
    if isfield(data, 'randomSeed') && ~isnan(data.randomSeed)
        prevRng = rng;
        rng(data.randomSeed, 'twister');
        rngCleanup = onCleanup(@() rng(prevRng)); %#ok<NASGU>
    else
        rngCleanup = [];
    end

    for i = 1:numel(data.speciesNames)
        h = data.scatterHandles(i);
        if ~isvalid(h)
            continue;
        end
        if data.visible(i)
            idx = data.speciesIndices{i};
            if data.clipToBoundingBox
                inside = isInsideBBox(data.coords(idx, :), data.bbox);
                idx = idx(inside);
            end
            n = computeSampleCount(numel(idx), data.sampleFractions(i));
            if n == 0
                set(h, 'XData', [], 'YData', [], 'ZData', [], 'Visible', 'off');
                continue;
            end
            if n < numel(idx)
                sampleIdx = idx(randperm(numel(idx), n));
            else
                sampleIdx = idx;
            end
            pts = data.coords(sampleIdx, :);
            set(h, 'XData', pts(:, 1), 'YData', pts(:, 2), 'ZData', pts(:, 3), ...
                'Visible', 'on', 'SizeData', data.markerSizes(i));
        else
            set(h, 'Visible', 'off', 'SizeData', data.markerSizes(i));
        end
    end
    applyBoundingBox(data.ax, data.bbox, data.fixBoundingBox);
end

function n = computeSampleCount(total, sampleValue)
    if sampleValue <= 1
        n = round(total * sampleValue);
    else
        n = round(sampleValue);
    end
    n = max(0, min(total, n));
end

function onSlider(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    value = max(0, min(1, src.Value));
    data.sampleValue = value;
    data.sampleFractions = initSampleFractions(data.speciesCounts, value);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'editBox')
        ctrls.editBox.String = num2str(value, '%.3f');
    end
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onEdit(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    value = str2double(src.String);
    if isnan(value) || value < 0
        value = data.sampleValue;
    end
    data.sampleValue = value;
    data.sampleFractions = initSampleFractions(data.speciesCounts, value);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'slider')
        ctrls.slider.Value = min(value, 1);
    end
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onMarkerSizeSlider(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    value = max(1, src.Value);
    data.markerSizeGlobal = value;
    data.markerSize = value;
    data.markerSizes = initMarkerSizes(numel(data.speciesNames), value);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'markerEdit')
        ctrls.markerEdit.String = num2str(value, '%.2f');
    end
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onMarkerSizeEdit(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    value = str2double(src.String);
    if isnan(value) || value <= 0
        value = data.markerSizeGlobal;
    end
    data.markerSizeGlobal = value;
    data.markerSize = value;
    data.markerSizes = initMarkerSizes(numel(data.speciesNames), value);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'markerSlider') && isgraphics(ctrls.markerSlider)
        if value > ctrls.markerSlider.Max
            ctrls.markerSlider.Max = value;
        end
        if value < ctrls.markerSlider.Min
            ctrls.markerSlider.Min = value;
        end
        ctrls.markerSlider.Value = min(max(value, ctrls.markerSlider.Min), ctrls.markerSlider.Max);
    end
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onRandomSeedEdit(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    value = str2double(src.String);
    if isnan(value)
        value = data.randomSeed;
    end
    value = round(value);
    data.randomSeed = value;
    src.String = num2str(value, '%.0f');
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateScatter(controlFig);
end

function onRotationSlider(src, controlFig)
    angleDeg = max(0, min(360, src.Value));
    applyRotation(controlFig, angleDeg);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'rotationEdit') && isgraphics(ctrls.rotationEdit)
        ctrls.rotationEdit.String = num2str(angleDeg, '%.1f');
    end
end

function onRotationEdit(src, controlFig)
    angleDeg = str2double(src.String);
    if isnan(angleDeg)
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        angleDeg = data.rotationAngle;
    end
    angleDeg = max(0, min(360, angleDeg));
    applyRotation(controlFig, angleDeg);
    src.String = num2str(angleDeg, '%.1f');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'rotationSlider') && isgraphics(ctrls.rotationSlider)
        ctrls.rotationSlider.Value = angleDeg;
    end
end

function applyRotation(controlFig, angleDeg)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data.rotationAngle = angleDeg;
    if isfield(data, 'coordsOriginal')
        data.coords = rotateCoordsZ(data.coordsOriginal, angleDeg);
    end
    if ~isfield(data, 'bboxUserEdited')
        data.bboxUserEdited = false;
    end
    if ~data.bboxUserEdited
        data.bboxOriginal = computeBoundingBox(data.coords, data.bboxPadding);
        data.bbox = data.bboxOriginal;
    end
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateBoundingBoxControls(controlFig);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onExportImages(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end

    filterSpec = {'*.png;*.tif;*.jpg;*.pdf;*.eps;*.svg', 'All Supported (*.png, *.tif, *.jpg, *.pdf, *.eps, *.svg)'; ...
        '*.png', 'PNG (*.png)'; '*.tif', 'TIFF (*.tif)'; '*.jpg', 'JPEG (*.jpg)'; ...
        '*.pdf', 'PDF (*.pdf)'; '*.eps', 'EPS (*.eps)'; '*.svg', 'SVG (*.svg)'};

    [fileName, pathName, filterIndex] = uiputfile(filterSpec, 'Export Images');
    if isequal(fileName, 0)
        return;
    end

    [baseName, ext] = fileparts(fileName);
    ext = lower(string(ext));
    if ~ismember(ext, [".png", ".tif", ".tiff", ".jpg", ".jpeg", ".pdf", ".eps", ".svg"])
        ext = resolveExportExtension(filterIndex);
    end

    visibleIdx = find(data.visible);
    if isempty(visibleIdx)
        return;
    end

    originalVisible = data.visible;
    fig = ancestor(data.ax, 'figure');
    if isgraphics(fig)
        figure(fig);
    end

    for k = 1:numel(visibleIdx)
        i = visibleIdx(k);
        data.visible(:) = false;
        data.visible(i) = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateScatter(controlFig);
        drawnow;

        tag = sanitizeFileTag(data.displayNames(i));
        if isempty(baseName)
            outName = tag;
        else
            outName = baseName + "_" + tag;
        end
        outFile = fullfile(pathName, outName + ext);

        if ismember(ext, [".pdf", ".eps", ".svg"])
            exportgraphics(data.ax, outFile, 'ContentType', 'vector');
        else
            exportgraphics(data.ax, outFile, 'Resolution', 300);
        end
    end

    data.visible = originalVisible;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onExportTurntable(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end

    filterSpec = {'*.avi;*.mp4', 'Video Files (*.avi, *.mp4)'; ...
        '*.avi', 'AVI (*.avi)'; '*.mp4', 'MP4 (*.mp4)'};
    [fileName, pathName] = uiputfile(filterSpec, 'Export Turntable');
    if isequal(fileName, 0)
        return;
    end

    [baseName, ext] = fileparts(fileName);
    if isempty(ext)
        ext = '.avi';
    end

    visibleIdx = find(data.visible);
    if isempty(visibleIdx)
        return;
    end

    originalVisible = data.visible;
    fig = ancestor(data.ax, 'figure');
    if isgraphics(fig)
        figure(fig);
    end

    for k = 1:numel(visibleIdx)
        i = visibleIdx(k);
        data.visible(:) = false;
        data.visible(i) = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateScatter(controlFig);
        drawnow;

        tag = sanitizeFileTag(data.displayNames(i));
        if isempty(baseName)
            outName = tag;
        else
            outName = baseName + "_" + tag;
        end
        outFile = fullfile(pathName, outName + ext);

        movieCreateTurntableAnimation(0.5, 30, outFile);
    end

    data.visible = originalVisible;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onTableEdit(src, evd, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    row = evd.Indices(1);
    col = evd.Indices(2);
    if col == 1
        data.visible(row) = logical(evd.NewData);
    elseif col == 4
        newCount = max(0, round(evd.NewData));
        total = data.speciesCounts(row);
        data.sampleFractions(row) = min(1, newCount / max(total, 1));
    elseif col == 5
        newFrac = max(0, min(1, evd.NewData));
        data.sampleFractions(row) = newFrac;
    elseif col == 6
        newSize = evd.NewData;
        if ~isnumeric(newSize) || isnan(newSize)
            newSize = data.markerSizes(row);
        end
        newSize = max(0.1, newSize);
        data.markerSizes(row) = newSize;
    end
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onShowAll(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data.visible(:) = true;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onHideAll(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data.visible(:) = false;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function onGroupingChanged(controlFig, ~, modeGroup)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(ctrls)
        return;
    end
    if nargin >= 3 && ~isempty(modeGroup) && isgraphics(modeGroup)
        selected = modeGroup.SelectedObject;
        if isgraphics(selected) && strcmpi(selected.String, 'Atomic')
            data.mode = "atomic";
        else
            data.mode = "ionic";
        end
    elseif isfield(ctrls, 'modeGroup') && isgraphics(ctrls.modeGroup)
        selected = ctrls.modeGroup.SelectedObject;
        if isgraphics(selected) && strcmpi(selected.String, 'Atomic')
            data.mode = "atomic";
        else
            data.mode = "ionic";
        end
    end
    if isfield(ctrls, 'splitIsotope') && isgraphics(ctrls.splitIsotope)
        data.splitIsotope = logical(ctrls.splitIsotope.Value) && data.hasIsotope;
    end
    if isfield(ctrls, 'splitCharge') && isgraphics(ctrls.splitCharge)
        data.splitCharge = logical(ctrls.splitCharge.Value) && data.hasCharge;
    end

    [speciesNames, baseNames, displayNames, speciesIndices, speciesCounts] = computeGroups( ...
        data.pos, data.mode, data.splitIsotope, data.splitCharge, data.showUnranged);
    colors = mapColors(baseNames, data.colorScheme);
    data.markerSizes = initMarkerSizes(numel(speciesNames), data.markerSizeGlobal);

    if ~isempty(data.scatterHandles)
        delete(data.scatterHandles(ishghandle(data.scatterHandles)));
    end
    axes(data.ax);
    holdState = ishold(data.ax);
    hold(data.ax, 'on');
    data.scatterHandles = createScatterHandles(data.ax, displayNames, colors, data.markerSizes);
    if ~holdState
        hold(data.ax, 'off');
    end

    data.speciesNames = speciesNames;
    data.baseNames = baseNames;
    data.displayNames = displayNames;
    data.speciesIndices = speciesIndices;
    data.speciesCounts = speciesCounts;
    data.colors = colors;
    data.visible = true(numel(speciesNames), 1);
    data.sampleFractions = initSampleFractions(data.speciesCounts, data.sampleValue);

    updateTable(controlFig);

    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateScatter(controlFig);
end

function handles = createScatterHandles(ax, displayNames, colors, markerSizes)
    handles = gobjects(numel(displayNames), 1);
    for i = 1:numel(displayNames)
        if numel(markerSizes) >= i
            sizeValue = markerSizes(i);
        else
            sizeValue = markerSizes(1);
        end
        handles(i) = scatter3(ax, nan, nan, nan, sizeValue, ...
            'filled', 'MarkerFaceColor', colors(i, :), ...
            'MarkerEdgeColor', colors(i, :), ...
            'DisplayName', char(displayNames(i)));
    end
end

function bbox = computeBoundingBox(coords, paddingFraction)
    minVals = min(coords, [], 1);
    maxVals = max(coords, [], 1);
    delta = maxVals - minVals;
    if paddingFraction > 0
        minVals = minVals - delta * paddingFraction;
        maxVals = maxVals + delta * paddingFraction;
    end
    bbox = struct('xmin', minVals(1), 'xmax', maxVals(1), ...
        'ymin', minVals(2), 'ymax', maxVals(2), ...
        'zmin', minVals(3), 'zmax', maxVals(3));
end

function coordsRot = rotateCoordsZ(coords, angleDeg)
    c = cosd(angleDeg);
    s = sind(angleDeg);
    rot = [c -s 0; s c 0; 0 0 1];
    coordsRot = coords * rot;
end

function applyBoundingBox(ax, bbox, fixBoundingBox)
    if isempty(ax) || ~isgraphics(ax)
        return;
    end
    if fixBoundingBox
        xlim(ax, [bbox.xmin bbox.xmax]);
        ylim(ax, [bbox.ymin bbox.ymax]);
        zlim(ax, [bbox.zmin bbox.zmax]);
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual';
    else
        ax.XLimMode = 'auto';
        ax.YLimMode = 'auto';
        ax.ZLimMode = 'auto';
    end
end

function inside = isInsideBBox(coords, bbox)
    inside = coords(:, 1) >= bbox.xmin & coords(:, 1) <= bbox.xmax & ...
        coords(:, 2) >= bbox.ymin & coords(:, 2) <= bbox.ymax & ...
        coords(:, 3) >= bbox.zmin & coords(:, 3) <= bbox.zmax;
end

function [center, span] = bboxToCenterSpan(bbox)
    center = [(bbox.xmin + bbox.xmax) / 2, ...
        (bbox.ymin + bbox.ymax) / 2, ...
        (bbox.zmin + bbox.zmax) / 2];
    span = [bbox.xmax - bbox.xmin, ...
        bbox.ymax - bbox.ymin, ...
        bbox.zmax - bbox.zmin];
end

function bbox = centerSpanToBbox(center, span)
    halfSpan = span / 2;
    bbox = struct('xmin', center(1) - halfSpan(1), 'xmax', center(1) + halfSpan(1), ...
        'ymin', center(2) - halfSpan(2), 'ymax', center(2) + halfSpan(2), ...
        'zmin', center(3) - halfSpan(3), 'zmax', center(3) + halfSpan(3));
end

function [editMin, editMax] = makeBBoxRow(parent, labelText, y, valueMin, valueMax, callbackFn)
    uicontrol(parent, 'Style', 'text', ...
        'String', labelText, ...
        'Units', 'normalized', ...
        'Position', [0.02 y 0.05 0.2], ...
        'HorizontalAlignment', 'left');
    editMin = uicontrol(parent, 'Style', 'edit', ...
        'String', num2str(valueMin, '%.4f'), ...
        'Units', 'normalized', ...
        'Position', [0.08 y 0.4 0.22], ...
        'Callback', callbackFn);
    editMax = uicontrol(parent, 'Style', 'edit', ...
        'String', num2str(valueMax, '%.4f'), ...
        'Units', 'normalized', ...
        'Position', [0.55 y 0.4 0.22], ...
        'Callback', callbackFn);
end

function onBoundingBoxChanged(controlFig, action)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(data) || isempty(ctrls)
        return;
    end

    switch action
        case 'toggle'
            if isfield(ctrls, 'bboxFix') && isgraphics(ctrls.bboxFix)
                data.fixBoundingBox = logical(ctrls.bboxFix.Value);
            end
        case 'clip'
            if isfield(ctrls, 'bboxClip') && isgraphics(ctrls.bboxClip)
                data.clipToBoundingBox = logical(ctrls.bboxClip.Value);
            end
        case 'span'
            if isfield(ctrls, 'bboxSpan') && isgraphics(ctrls.bboxSpan)
                data.bboxUseSpan = logical(ctrls.bboxSpan.Value);
            end
        case 'current'
            if isgraphics(data.ax)
                data.bbox = struct('xmin', data.ax.XLim(1), 'xmax', data.ax.XLim(2), ...
                    'ymin', data.ax.YLim(1), 'ymax', data.ax.YLim(2), ...
                    'zmin', data.ax.ZLim(1), 'zmax', data.ax.ZLim(2));
                data.bboxUserEdited = true;
            end
        case 'reset'
            data.bbox = data.bboxOriginal;
            data.bboxUserEdited = false;
        case 'edit'
            data.bbox = readBoundingBoxEdits(ctrls, data.bbox, data.bboxUseSpan);
            data.bboxUserEdited = true;
    end

    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateBoundingBoxControls(controlFig);
    applyBoundingBox(data.ax, data.bbox, data.fixBoundingBox);
    updateTable(controlFig);
    updateScatter(controlFig);
end

function bbox = readBoundingBoxEdits(ctrls, bbox, useSpan)
    if ~isfield(ctrls, 'bboxEdits')
        return;
    end
    edits = ctrls.bboxEdits;
    if useSpan
        center = [readEditValue(edits.xmin, mean([bbox.xmin bbox.xmax])), ...
            readEditValue(edits.ymin, mean([bbox.ymin bbox.ymax])), ...
            readEditValue(edits.zmin, mean([bbox.zmin bbox.zmax]))];
        span = [readEditValue(edits.xmax, bbox.xmax - bbox.xmin), ...
            readEditValue(edits.ymax, bbox.ymax - bbox.ymin), ...
            readEditValue(edits.zmax, bbox.zmax - bbox.zmin)];
        span = max(span, 0);
        bbox = centerSpanToBbox(center, span);
    else
        bbox.xmin = readEditValue(edits.xmin, bbox.xmin);
        bbox.xmax = readEditValue(edits.xmax, bbox.xmax);
        bbox.ymin = readEditValue(edits.ymin, bbox.ymin);
        bbox.ymax = readEditValue(edits.ymax, bbox.ymax);
        bbox.zmin = readEditValue(edits.zmin, bbox.zmin);
        bbox.zmax = readEditValue(edits.zmax, bbox.zmax);
        if bbox.xmin > bbox.xmax
            [bbox.xmin, bbox.xmax] = deal(bbox.xmax, bbox.xmin);
        end
        if bbox.ymin > bbox.ymax
            [bbox.ymin, bbox.ymax] = deal(bbox.ymax, bbox.ymin);
        end
        if bbox.zmin > bbox.zmax
            [bbox.zmin, bbox.zmax] = deal(bbox.zmax, bbox.zmin);
        end
    end
end

function val = readEditValue(h, fallback)
    if ~isgraphics(h)
        val = fallback;
        return;
    end
    val = str2double(h.String);
    if isnan(val)
        val = fallback;
    end
end

function updateBoundingBoxControls(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(data) || isempty(ctrls) || ~isfield(ctrls, 'bboxEdits')
        return;
    end
    edits = ctrls.bboxEdits;
    if data.bboxUseSpan
        [center, span] = bboxToCenterSpan(data.bbox);
        if isgraphics(edits.xmin), edits.xmin.String = num2str(center(1), '%.4f'); end
        if isgraphics(edits.xmax), edits.xmax.String = num2str(span(1), '%.4f'); end
        if isgraphics(edits.ymin), edits.ymin.String = num2str(center(2), '%.4f'); end
        if isgraphics(edits.ymax), edits.ymax.String = num2str(span(2), '%.4f'); end
        if isgraphics(edits.zmin), edits.zmin.String = num2str(center(3), '%.4f'); end
        if isgraphics(edits.zmax), edits.zmax.String = num2str(span(3), '%.4f'); end
    else
        if isgraphics(edits.xmin), edits.xmin.String = num2str(data.bbox.xmin, '%.4f'); end
        if isgraphics(edits.xmax), edits.xmax.String = num2str(data.bbox.xmax, '%.4f'); end
        if isgraphics(edits.ymin), edits.ymin.String = num2str(data.bbox.ymin, '%.4f'); end
        if isgraphics(edits.ymax), edits.ymax.String = num2str(data.bbox.ymax, '%.4f'); end
        if isgraphics(edits.zmin), edits.zmin.String = num2str(data.bbox.zmin, '%.4f'); end
        if isgraphics(edits.zmax), edits.zmax.String = num2str(data.bbox.zmax, '%.4f'); end
    end
    if isfield(ctrls, 'bboxFix') && isgraphics(ctrls.bboxFix)
        ctrls.bboxFix.Value = data.fixBoundingBox;
    end
    if isfield(ctrls, 'bboxClip') && isgraphics(ctrls.bboxClip)
        ctrls.bboxClip.Value = data.clipToBoundingBox;
    end
    if isfield(ctrls, 'bboxSpan') && isgraphics(ctrls.bboxSpan)
        ctrls.bboxSpan.Value = data.bboxUseSpan;
    end
end

function sampleValue = resolveInitialSample(totalCount, requestedSample)
    if ~isnan(requestedSample)
        if requestedSample < 0
            error('scatterPlotPosWidget:invalidSample', ...
                'Sample must be nonnegative.');
        end
        sampleValue = requestedSample;
        return;
    end
    if totalCount <= 1e6
        sampleValue = 1;
    else
        sampleValue = 1e6 / totalCount;
    end
end

function displayName = formatDisplayName(baseName, pos, rowIdx, splitIsotope, splitCharge, showUnranged)
    name = string(baseName);
    if ~showUnranged && (name == "" || name == "missing")
        displayName = "";
        return;
    end
    if name == "" || name == "missing"
        name = "unranged";
    end
    if splitCharge
        name = stripChargeFromName(name);
    end

    isoPart = "";
    if splitIsotope && ismember('isotope', pos.Properties.VariableNames)
        iso = pos.isotope(rowIdx);
        if ~ismissing(iso)
            isoPart = string(iso);
        end
    end

    chargePart = "";
    if splitCharge && ismember('chargeState', pos.Properties.VariableNames)
        cs = pos.chargeState(rowIdx);
        if ~isnan(cs)
            n = round(cs);
            if n > 0
                chargePart = repmat("+", 1, n);
            end
        end
    end

    displayName = isoPart + name + chargePart;
end

function nameOut = stripChargeFromName(nameIn)
    nameOut = string(nameIn);
    nameOut = regexprep(nameOut, '(\\d*\\+)+$', '');
    nameOut = regexprep(nameOut, '\\++$', '');
    nameOut = strtrim(nameOut);
end

function state = boolToOnOff(value)
    if value
        state = 'on';
    else
        state = 'off';
    end
end

function fractions = initSampleFractions(counts, sampleValue)
    if sampleValue <= 1
        fractions = repmat(sampleValue, size(counts));
    else
        fractions = min(1, sampleValue ./ max(counts, 1));
    end
end

function markerSizeGlobal = resolveMarkerSizeGlobal(requested, fallback)
    if ~isnan(requested)
        if requested <= 0
            error('scatterPlotPosWidget:invalidMarkerSize', ...
                'Marker size must be positive.');
        end
        markerSizeGlobal = requested;
    else
        markerSizeGlobal = fallback;
    end
end

function sizes = initMarkerSizes(count, markerSizeGlobal)
    if isempty(markerSizeGlobal) || isnan(markerSizeGlobal) || markerSizeGlobal <= 0
        markerSizeGlobal = 1;
    end
    sizes = repmat(markerSizeGlobal, count, 1);
end

function tag = sanitizeFileTag(nameIn)
    tag = string(nameIn);
    if strlength(tag) == 0
        tag = "species";
    end
    tag = regexprep(tag, '[^a-zA-Z0-9]+', '_');
    tag = regexprep(tag, '^_+|_+$', '');
    if strlength(tag) == 0
        tag = "species";
    end
end

function ext = resolveExportExtension(filterIndex)
    switch filterIndex
        case 2
            ext = '.png';
        case 3
            ext = '.tif';
        case 4
            ext = '.jpg';
        case 5
            ext = '.pdf';
        case 6
            ext = '.eps';
        case 7
            ext = '.svg';
        otherwise
            ext = '.png';
    end
end

function tableData = buildTableData(data)
    tableData = cell(numel(data.speciesNames), 6);
    for i = 1:numel(data.speciesNames)
        total = data.speciesCounts(i);
        frac = data.sampleFractions(i);
        if data.clipToBoundingBox
            idx = data.speciesIndices{i};
            inside = isInsideBBox(data.coords(idx, :), data.bbox);
            totalShown = sum(inside);
        else
            totalShown = total;
        end
        shown = computeSampleCount(totalShown, frac);
        tableData{i, 1} = data.visible(i);
    if isfield(data, 'displayNames') && numel(data.displayNames) >= i
        tableData{i, 2} = char(data.displayNames(i));
    else
        tableData{i, 2} = char(data.speciesNames(i));
    end
        tableData{i, 3} = total;
        tableData{i, 4} = shown;
        tableData{i, 5} = frac;
        if isfield(data, 'markerSizes') && numel(data.markerSizes) >= i
            tableData{i, 6} = data.markerSizes(i);
        else
            tableData{i, 6} = data.markerSizeGlobal;
        end
    end
end

function updateTable(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(ctrls) || ~isfield(ctrls, 'table')
        return;
    end
    ctrls.table.Data = buildTableData(data);
    updateBoundingBoxControls(controlFig);
    drawnow;
end
