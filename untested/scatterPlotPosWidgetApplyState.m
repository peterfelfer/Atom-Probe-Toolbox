function controlFig = scatterPlotPosWidgetApplyState(controlFig, state)
% SCATTERPLOTPOSWIDGETAPPLYSTATE Apply a saved widget state to the UI.
%
% controlFig = scatterPlotPosWidgetApplyState(controlFig, state)
%
% Requires an active scatterPlotPosWidget control figure. To continue from
% a saved visualization axis, reopen the widget on that axis first:
%   [~, ~, controlFig] = scatterPlotPosWidget(pos, colorScheme, 'axes', ax);
%   scatterPlotPosWidgetApplyState(controlFig, state);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 2
    error('scatterPlotPosWidgetApplyState:missingState', ...
        'State struct is required.');
end

if nargin < 1 || isempty(controlFig)
    controlFig = gcf;
end

if isgraphics(controlFig, 'axes') || isgraphics(controlFig, 'patch')
    controlFig = ancestor(controlFig, 'figure');
end

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data) || ~isfield(data, 'fn')
    error('scatterPlotPosWidgetApplyState:invalidHandle', ...
        'No widget state found on the specified figure.');
end

ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
fn = data.fn;

if isfield(state, 'mode')
    data.mode = string(state.mode);
end
if isfield(state, 'splitIsotope')
    data.splitIsotope = logical(state.splitIsotope) && data.hasIsotope;
end
if isfield(state, 'splitCharge')
    data.splitCharge = logical(state.splitCharge) && data.hasCharge;
end
if isfield(state, 'showUnranged')
    data.showUnranged = logical(state.showUnranged);
end
if isfield(state, 'sampleValue')
    data.sampleValue = double(state.sampleValue);
end
if isfield(state, 'markerSize')
    data.markerSize = double(state.markerSize);
end
if isfield(state, 'markerSizeGlobal')
    data.markerSizeGlobal = double(state.markerSizeGlobal);
elseif isfield(state, 'markerSize')
    data.markerSizeGlobal = double(state.markerSize);
end
if isfield(data, 'markerSizeGlobal')
    data.markerSize = data.markerSizeGlobal;
end
if isfield(state, 'randomSeed')
    data.randomSeed = double(state.randomSeed);
end
if isfield(state, 'rotationAngle')
    data.rotationAngle = double(state.rotationAngle);
end
if isfield(state, 'liveUpdate')
    data.liveUpdate = logical(state.liveUpdate);
end
if isfield(state, 'showAxisLabels')
    data.showAxisLabels = logical(state.showAxisLabels);
end
if isfield(state, 'showAxisTicks')
    data.showAxisTicks = logical(state.showAxisTicks);
end
if isfield(state, 'showScaleCube')
    data.showScaleCube = logical(state.showScaleCube);
end
if isfield(state, 'scaleCubeSize')
    data.scaleCubeSize = double(state.scaleCubeSize);
end
if isfield(state, 'searchFilter')
    data.searchFilter = string(state.searchFilter);
end
if isfield(state, 'sortMode')
    data.sortMode = string(state.sortMode);
end
if isfield(state, 'boundingBox')
    data.bbox = state.boundingBox;
end
if isfield(state, 'fixBoundingBox')
    data.fixBoundingBox = logical(state.fixBoundingBox);
end
if isfield(state, 'clipToBoundingBox')
    data.clipToBoundingBox = logical(state.clipToBoundingBox);
end
if isfield(state, 'bboxUseSpan')
    data.bboxUseSpan = logical(state.bboxUseSpan);
end
if isfield(state, 'boundingBox')
    data.bboxUserEdited = true;
else
    data.bboxUserEdited = false;
end

% Update UI controls
if ~isempty(ctrls)
    if isfield(ctrls, 'modeGroup') && isgraphics(ctrls.modeGroup)
        children = ctrls.modeGroup.Children;
        for i = 1:numel(children)
            if strcmpi(children(i).String, 'Atomic') && data.mode == "atomic"
                ctrls.modeGroup.SelectedObject = children(i);
            elseif strcmpi(children(i).String, 'Ionic') && data.mode == "ionic"
                ctrls.modeGroup.SelectedObject = children(i);
            end
        end
    end
    if isfield(ctrls, 'splitIsotope') && isgraphics(ctrls.splitIsotope)
        ctrls.splitIsotope.Value = data.splitIsotope;
    end
    if isfield(ctrls, 'splitCharge') && isgraphics(ctrls.splitCharge)
        ctrls.splitCharge.Value = data.splitCharge;
    end
    if isfield(ctrls, 'slider') && isgraphics(ctrls.slider)
        ctrls.slider.Value = min(data.sampleValue, 1);
    end
    if isfield(ctrls, 'editBox') && isgraphics(ctrls.editBox)
        ctrls.editBox.String = num2str(data.sampleValue, '%.3f');
    end
    if isfield(ctrls, 'markerSlider') && isgraphics(ctrls.markerSlider)
        if data.markerSizeGlobal > ctrls.markerSlider.Max
            ctrls.markerSlider.Max = data.markerSizeGlobal;
        end
        ctrls.markerSlider.Value = min(max(data.markerSizeGlobal, ctrls.markerSlider.Min), ctrls.markerSlider.Max);
    end
    if isfield(ctrls, 'markerEdit') && isgraphics(ctrls.markerEdit)
        ctrls.markerEdit.String = num2str(data.markerSizeGlobal, '%.2f');
    end
    if isfield(ctrls, 'seedEdit') && isgraphics(ctrls.seedEdit)
        ctrls.seedEdit.String = num2str(data.randomSeed, '%.0f');
    end
    if isfield(ctrls, 'liveUpdateCb') && isgraphics(ctrls.liveUpdateCb)
        ctrls.liveUpdateCb.Value = data.liveUpdate;
    end
    if isfield(ctrls, 'rotationSlider') && isgraphics(ctrls.rotationSlider)
        ctrls.rotationSlider.Value = max(0, min(360, data.rotationAngle));
    end
    if isfield(ctrls, 'rotationEdit') && isgraphics(ctrls.rotationEdit)
        ctrls.rotationEdit.String = num2str(max(0, min(360, data.rotationAngle)), '%.1f');
    end
    if isfield(ctrls, 'bboxClip') && isgraphics(ctrls.bboxClip)
        ctrls.bboxClip.Value = data.clipToBoundingBox;
    end
    if isfield(ctrls, 'bboxSpan') && isgraphics(ctrls.bboxSpan)
        ctrls.bboxSpan.Value = data.bboxUseSpan;
    end
    if isfield(ctrls, 'showAxisLabels') && isgraphics(ctrls.showAxisLabels)
        ctrls.showAxisLabels.Value = data.showAxisLabels;
    end
    if isfield(ctrls, 'showAxisTicks') && isgraphics(ctrls.showAxisTicks)
        ctrls.showAxisTicks.Value = data.showAxisTicks;
    end
    if isfield(ctrls, 'showScaleCube') && isgraphics(ctrls.showScaleCube)
        ctrls.showScaleCube.Value = data.showScaleCube;
    end
    if isfield(ctrls, 'scaleCubeSizeEdit') && isgraphics(ctrls.scaleCubeSizeEdit)
        ctrls.scaleCubeSizeEdit.String = num2str(data.scaleCubeSize);
    end
    if isfield(ctrls, 'searchBox') && isgraphics(ctrls.searchBox)
        ctrls.searchBox.String = char(data.searchFilter);
    end
    if isfield(ctrls, 'sortDropdown') && isgraphics(ctrls.sortDropdown)
        sortModes = {'none', 'name_asc', 'name_desc', 'count_desc', 'count_asc', 'mass_desc', 'mass_asc'};
        idxSort = find(strcmp(sortModes, char(data.sortMode)), 1, 'first');
        if isempty(idxSort)
            idxSort = 1;
        end
        ctrls.sortDropdown.Value = idxSort;
    end
end

% Recompute grouping
[speciesNames, baseNames, displayNames, speciesIndices, speciesCounts] = fn.computeGroups( ...
    data.pos, data.mode, data.splitIsotope, data.splitCharge, data.showUnranged);
colors = fn.mapColors(baseNames, data.colorScheme);
data.markerSizes = fn.initMarkerSizes(numel(speciesNames), data.markerSizeGlobal);
data.markerSizeLinked = true(numel(speciesNames), 1);
if isfield(data, 'coordsOriginal')
    data.coords = fn.rotateCoordsZ(data.coordsOriginal, data.rotationAngle);
end

% Rebuild scatter handles
if ~isempty(data.scatterHandles)
    delete(data.scatterHandles(ishghandle(data.scatterHandles)));
end
axes(data.ax);
holdState = ishold(data.ax);
hold(data.ax, 'on');
data.scatterHandles = fn.createScatterHandles(data.ax, displayNames, colors, data.markerSizes);
if ~holdState
    hold(data.ax, 'off');
end

data.speciesNames = speciesNames;
data.baseNames = baseNames;
data.displayNames = displayNames;
data.speciesIndices = speciesIndices;
data.speciesCounts = speciesCounts;
data.colors = colors;

% Apply per-species visibility/fractions by name
data.visible = true(numel(speciesNames), 1);
data.sampleFractions = fn.initSampleFractions(speciesCounts, data.sampleValue);
savedMarkerSize = [];
savedMarkerLinked = [];

if isfield(state, 'species') && isfield(state.species, 'name')
    savedNames = string(state.species.name);
    savedVisible = logical(state.species.visible);
    savedFraction = double(state.species.fraction);
    if isfield(state.species, 'markerSize')
        savedMarkerSize = double(state.species.markerSize);
    else
        savedMarkerSize = [];
    end
    if isfield(state.species, 'markerSizeLinked')
        savedMarkerLinked = logical(state.species.markerSizeLinked);
    else
        savedMarkerLinked = [];
    end
    if isfield(state.species, 'color')
        savedColors = double(state.species.color);
    else
        savedColors = [];
    end
    for i = 1:numel(speciesNames)
        idx = find(savedNames == speciesNames(i), 1, 'first');
        if ~isempty(idx)
            if idx <= numel(savedVisible)
                data.visible(i) = savedVisible(idx);
            end
            if idx <= numel(savedFraction) && ~isnan(savedFraction(idx))
                data.sampleFractions(i) = savedFraction(idx);
            end
            if idx <= numel(savedMarkerSize) && ~isnan(savedMarkerSize(idx))
                data.markerSizes(i) = savedMarkerSize(idx);
            end
            if idx <= numel(savedMarkerLinked)
                data.markerSizeLinked(i) = savedMarkerLinked(idx);
            end
            if ~isempty(savedColors) && size(savedColors, 1) >= idx && size(savedColors, 2) == 3
                data.colors(i, :) = savedColors(idx, :);
            end
        end
    end
end

% Backward compatibility: infer marker-size links from saved per-species sizes
% when an older state did not store explicit link flags.
if isempty(savedMarkerLinked) && ~isempty(savedMarkerSize)
    tol = max(1e-9, 1e-6 * max(1, data.markerSizeGlobal));
    data.markerSizeLinked = abs(data.markerSizes(:) - data.markerSizeGlobal) <= tol;
end

% Linked species always follow the current global marker size.
data.markerSizeLinked = logical(data.markerSizeLinked(:));
linkedIdx = data.markerSizeLinked;
if any(linkedIdx)
    data.markerSizes(linkedIdx) = data.markerSizeGlobal;
end

for i = 1:numel(data.scatterHandles)
    if isgraphics(data.scatterHandles(i))
        set(data.scatterHandles(i), ...
            'MarkerFaceColor', data.colors(i, :), ...
            'MarkerEdgeColor', data.colors(i, :), ...
            'SizeData', data.markerSizes(i));
    end
end

setappdata(controlFig, 'scatterPlotPosWidget', data);
fn.updateTable(controlFig);
fn.updateScatter(controlFig);
if isfield(data, 'bbox') && isgraphics(data.ax)
    if data.fixBoundingBox
        xlim(data.ax, [data.bbox.xmin data.bbox.xmax]);
        ylim(data.ax, [data.bbox.ymin data.bbox.ymax]);
        zlim(data.ax, [data.bbox.zmin data.bbox.zmax]);
        data.ax.XLimMode = 'manual';
        data.ax.YLimMode = 'manual';
        data.ax.ZLimMode = 'manual';
    else
        data.ax.XLimMode = 'auto';
        data.ax.YLimMode = 'auto';
        data.ax.ZLimMode = 'auto';
    end
end

% Restore camera
if isfield(state, 'camera') && isgraphics(data.ax)
    cam = state.camera;
    if isfield(cam, 'CameraPosition'), data.ax.CameraPosition = cam.CameraPosition; end
    if isfield(cam, 'CameraTarget'), data.ax.CameraTarget = cam.CameraTarget; end
    if isfield(cam, 'CameraUpVector'), data.ax.CameraUpVector = cam.CameraUpVector; end
    if isfield(cam, 'CameraViewAngle'), data.ax.CameraViewAngle = cam.CameraViewAngle; end
    if isfield(cam, 'ViewAngle'), data.ax.ViewAngle = cam.ViewAngle; end
    if isfield(cam, 'Projection'), data.ax.Projection = cam.Projection; end
    if isfield(cam, 'View') && isnumeric(cam.View) && numel(cam.View) == 2
        view(data.ax, cam.View(1), cam.View(2));
    end
end

% Apply axis labels/ticks via callbacks where available.
if ~isempty(ctrls)
    if isfield(ctrls, 'showAxisLabels') && isgraphics(ctrls.showAxisLabels)
        try
            feval(ctrls.showAxisLabels.Callback, ctrls.showAxisLabels, []);
        catch
        end
    end
    if isfield(ctrls, 'showAxisTicks') && isgraphics(ctrls.showAxisTicks)
        try
            feval(ctrls.showAxisTicks.Callback, ctrls.showAxisTicks, []);
        catch
        end
    end
    if isfield(ctrls, 'scaleCubeSizeEdit') && isgraphics(ctrls.scaleCubeSizeEdit)
        try
            feval(ctrls.scaleCubeSizeEdit.Callback, ctrls.scaleCubeSizeEdit, []);
        catch
        end
    end
    if isfield(ctrls, 'showScaleCube') && isgraphics(ctrls.showScaleCube)
        try
            feval(ctrls.showScaleCube.Callback, ctrls.showScaleCube, []);
        catch
        end
    end
    if isfield(ctrls, 'liveUpdateCb') && isgraphics(ctrls.liveUpdateCb)
        try
            feval(ctrls.liveUpdateCb.Callback, ctrls.liveUpdateCb, []);
        catch
        end
    end
    if isfield(ctrls, 'sortDropdown') && isgraphics(ctrls.sortDropdown)
        try
            feval(ctrls.sortDropdown.Callback, ctrls.sortDropdown, []);
        catch
        end
    end
    if isfield(ctrls, 'searchBox') && isgraphics(ctrls.searchBox)
        try
            feval(ctrls.searchBox.Callback, ctrls.searchBox, []);
        catch
        end
    end
end

% Persist applied state in the associated axes.
if isgraphics(data.ax, 'axes')
    try
        currentState = scatterPlotPosWidgetGetState(controlFig);
        writeStateToAxis(data.ax, currentState);
    catch
    end
end
end

function writeStateToAxis(ax, state)
if ~isgraphics(ax, 'axes')
    return;
end
ud = ax.UserData;
if isempty(ud) || isstruct(ud)
    if isempty(ud)
        ud = struct();
    end
    ud.scatterPlotPosWidgetState = state;
    ax.UserData = ud;
else
    setappdata(ax, 'scatterPlotPosWidgetState', state);
end
end
