function controlFig = scatterPlotPosWidgetApplyState(controlFig, state)
% SCATTERPLOTPOSWIDGETAPPLYSTATE Apply a saved widget state to the UI.
%
% controlFig = scatterPlotPosWidgetApplyState(controlFig, state)
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
end

% Recompute grouping
[speciesNames, baseNames, displayNames, speciesIndices, speciesCounts] = fn.computeGroups( ...
    data.pos, data.mode, data.splitIsotope, data.splitCharge, data.showUnranged);
colors = fn.mapColors(baseNames, data.colorScheme);
data.markerSizes = fn.initMarkerSizes(numel(speciesNames), data.markerSizeGlobal);
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

if isfield(state, 'species') && isfield(state.species, 'name')
    savedNames = string(state.species.name);
    savedVisible = logical(state.species.visible);
    savedFraction = double(state.species.fraction);
    if isfield(state.species, 'markerSize')
        savedMarkerSize = double(state.species.markerSize);
    else
        savedMarkerSize = [];
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
        end
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
    if isfield(cam, 'ViewAngle'), data.ax.ViewAngle = cam.ViewAngle; end
    if isfield(cam, 'Projection'), data.ax.Projection = cam.Projection; end
end
end
