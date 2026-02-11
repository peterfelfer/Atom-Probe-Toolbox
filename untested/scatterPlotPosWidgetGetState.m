function state = scatterPlotPosWidgetGetState(controlFig)
% SCATTERPLOTPOSWIDGETGETSTATE Return widget state as a struct.
%
% state = scatterPlotPosWidgetGetState(controlFig)
% state = scatterPlotPosWidgetGetState()
% state = scatterPlotPosWidgetGetState(ax)
%
% Returns a struct suitable for storing in the workspace and later
% restoring with scatterPlotPosWidgetApplyState. If an axes handle is
% provided, returns the state persisted on that axes (if present).
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 1 || isempty(controlFig)
    controlFig = gcf;
end

ax = [];
if isgraphics(controlFig, 'axes') || isgraphics(controlFig, 'patch')
    ax = ancestor(controlFig, 'axes');
    controlFig = ancestor(controlFig, 'figure');
end

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data) && ~isempty(ax)
    state = readStateFromAxis(ax);
    if ~isempty(state)
        return;
    end
end
if isempty(data)
    error('scatterPlotPosWidgetGetState:invalidHandle', ...
        'No widget state found on the specified handle.');
end

state = struct();
state.version = 3;
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
if isfield(data, 'liveUpdate')
    state.liveUpdate = data.liveUpdate;
end
if isfield(data, 'showAxisLabels')
    state.showAxisLabels = data.showAxisLabels;
end
if isfield(data, 'showAxisTicks')
    state.showAxisTicks = data.showAxisTicks;
end
if isfield(data, 'showScaleCube')
    state.showScaleCube = data.showScaleCube;
end
if isfield(data, 'scaleCubeSize')
    state.scaleCubeSize = data.scaleCubeSize;
end
if isfield(data, 'searchFilter')
    state.searchFilter = data.searchFilter;
end
if isfield(data, 'sortMode')
    state.sortMode = data.sortMode;
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
if isfield(data, 'markerSizeLinked')
    species.markerSizeLinked = logical(data.markerSizeLinked);
end
if isfield(data, 'colors')
    species.color = data.colors;
end
state.species = species;
if isfield(data, 'pos')
    state.posSignature = localComputePosSignature(data.pos);
end

if isgraphics(data.ax)
    cam = struct();
    cam.CameraPosition = data.ax.CameraPosition;
    cam.CameraTarget = data.ax.CameraTarget;
    cam.CameraUpVector = data.ax.CameraUpVector;
    cam.CameraViewAngle = data.ax.CameraViewAngle;
    cam.Projection = data.ax.Projection;
    [az, el] = view(data.ax);
    cam.View = [az, el];
    state.camera = cam;
end
end

function state = readStateFromAxis(ax)
state = [];
if ~isgraphics(ax, 'axes')
    return;
end
ud = ax.UserData;
if isstruct(ud) && isfield(ud, 'scatterPlotPosWidgetState') && isstruct(ud.scatterPlotPosWidgetState)
    state = ud.scatterPlotPosWidgetState;
    return;
end
if isappdata(ax, 'scatterPlotPosWidgetState')
    candidate = getappdata(ax, 'scatterPlotPosWidgetState');
    if isstruct(candidate)
        state = candidate;
    end
end
end

function sig = localComputePosSignature(pos)
sig = struct();
sig.height = height(pos);
sig.width = width(pos);
sig.names = string(pos.Properties.VariableNames);
sampleCols = intersect({'x','y','z','mc','detx','dety','ionIdx'}, pos.Properties.VariableNames, 'stable');
sampleIdx = unique([1:min(5,height(pos)), max(1,height(pos)-4):height(pos)]);
acc = 0;
for i = 1:numel(sampleCols)
    v = pos.(sampleCols{i});
    if isnumeric(v)
        acc = acc + sum(double(v(sampleIdx)), 'omitnan');
    end
end
sig.sample = acc;
end
