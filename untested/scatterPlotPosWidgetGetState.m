function state = scatterPlotPosWidgetGetState(controlFig)
% SCATTERPLOTPOSWIDGETGETSTATE Return widget state as a struct.
%
% state = scatterPlotPosWidgetGetState(controlFig)
% state = scatterPlotPosWidgetGetState()
%
% Returns a struct suitable for storing in the workspace and later
% restoring with scatterPlotPosWidgetApplyState.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 1 || isempty(controlFig)
    controlFig = gcf;
end

if isgraphics(controlFig, 'axes') || isgraphics(controlFig, 'patch')
    controlFig = ancestor(controlFig, 'figure');
end

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data)
    error('scatterPlotPosWidgetGetState:invalidHandle', ...
        'No widget state found on the specified figure.');
end

state = struct();
state.version = 1;
state.mode = data.mode;
state.splitIsotope = data.splitIsotope;
state.splitCharge = data.splitCharge;
state.showUnranged = data.showUnranged;
state.sampleValue = data.sampleValue;
state.markerSize = data.markerSize;
if isfield(data, 'markerSizeGlobal')
    state.markerSizeGlobal = data.markerSizeGlobal;
end
if isfield(data, 'randomSeed')
    state.randomSeed = data.randomSeed;
end
if isfield(data, 'rotationAngle')
    state.rotationAngle = data.rotationAngle;
end
if isfield(data, 'bbox')
    state.boundingBox = data.bbox;
    state.fixBoundingBox = data.fixBoundingBox;
    if isfield(data, 'clipToBoundingBox')
        state.clipToBoundingBox = data.clipToBoundingBox;
    end
    if isfield(data, 'bboxUseSpan')
        state.bboxUseSpan = data.bboxUseSpan;
    end
end

species = struct();
species.name = string(data.speciesNames);
species.visible = logical(data.visible);
species.fraction = data.sampleFractions;
if isfield(data, 'markerSizes')
    species.markerSize = data.markerSizes;
end
state.species = species;

if isgraphics(data.ax)
    cam = struct();
    cam.CameraPosition = data.ax.CameraPosition;
    cam.CameraTarget = data.ax.CameraTarget;
    cam.CameraUpVector = data.ax.CameraUpVector;
    cam.ViewAngle = data.ax.ViewAngle;
    cam.Projection = data.ax.Projection;
    state.camera = cam;
end
end
