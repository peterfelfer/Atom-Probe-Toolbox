function [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme, options)
% SCATTERPLOTPOSWIDGET Interactive scatter plot with visibility/sample controls.
%
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme)
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, colorScheme, 'groupBy', 'ion')
% [scatterHandles, ax, controlFig, info] = scatterPlotPosWidget(pos, [], 'axes', ax)
%
% Creates a 3D scatter plot of APT data and a separate control window that
% toggles species visibility and adjusts the displayed sampling fraction.
% The widget persists visualization state on the target axes so reopening
% the widget on the same axes resumes the previous configuration.
%
% INPUT:
%   pos         table with x, y, z and a species column (ion/atom/isotope)
%   colorScheme table with columns: ion (categorical/string) and color (Nx3)
%
% OPTIONS:
%   'axes'         - Target axes for scatter plot (default: new figure)
%   'groupBy'      - 'ion' | 'atom' | 'auto' (default: 'auto')
%   'sample'       - Initial sample fraction (0..1) or count (>1)
%                    (default: 1 up to 1e6 ions, else 1e6/total ions)
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
%   'restoreStateFromAxis' - Restore saved state from ax.UserData (default: true)
%   'persistStateToAxis' - Persist state to ax.UserData during interaction (default: true)
%
% MARKER-SIZE LINKING:
%   Editing a species marker size in the table unlinks that species from
%   the global marker size slider. Use the table "Link" column to relink
%   species back to the global value.
%
% STATE PERSISTENCE:
%   Reopen on existing axes to continue a saved visualization:
%     [~, ax] = scatterPlotPosData(pos, colorScheme);
%     scatterPlotPosWidget(pos, colorScheme, 'axes', ax);
%
%   Save/restore explicit state structs:
%     state = scatterPlotPosWidgetGetState(controlFig);
%     save('state.mat', 'state');
%     load('state.mat', 'state');
%     scatterPlotPosWidgetApplyState(controlFig, state);
%
%   Load a saved profile/state from workspace into a running widget:
%     scatterPlotPosWidgetLoadState(controlFig, 'visualisationProfile');
%
%   Reattach widget to a previously saved/opened visualization axes:
%     fig = openfig('savedScatter.fig');
%     ax = findobj(fig, 'Type', 'axes');
%     scatterPlotPosWidget(pos, colorScheme, 'axes', ax(1));
%
% KEYBOARD SHORTCUTS (when control window is focused):
%   Space  - Toggle visibility of selected species
%   A      - Show all species
%   H      - Hide all species
%   R      - Reset view to default
%   1-9    - Toggle visibility of species 1-9
%   Ctrl+Z - Undo last change
%   Ctrl+Y - Redo last undone change
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
    options.linkToBase (1,1) logical = false
    options.linkedVar (1,1) string = ""
    options.restoreStateFromAxis (1,1) logical = true
    options.persistStateToAxis (1,1) logical = true
end

baseVarName = inputname(1);

% Input validation
validatePosTable(pos);

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
cleanupWidgetScatterHandles(ax);

markerSizeGlobal = resolveMarkerSizeGlobal(options.markerSizeGlobal, options.markerSize);
markerSizes = initMarkerSizes(numel(speciesNames), markerSizeGlobal);
markerSizeLinked = true(numel(speciesNames), 1);

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
data.scatterFig = fig;
data.controlFig = controlFig;
data.markerSize = markerSizeGlobal;
data.markerSizeGlobal = markerSizeGlobal;
data.markerSizes = markerSizes;
data.markerSizeLinked = markerSizeLinked;
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
data.linkToBase = options.linkToBase;
data.linkedVar = options.linkedVar;
data.linkKey = "";
data.linkSignature = [];
data.linkUpdating = false;
data.linkRefreshTimer = [];
data.restoreStateFromAxis = options.restoreStateFromAxis;
data.persistStateToAxis = options.persistStateToAxis;
data.suspendStatePersistence = false;

% Performance: deferred updates and sampling cache
data.liveUpdate = true;
data.samplingCache = cell(numel(speciesNames), 1);
data.samplingCacheValid = false(numel(speciesNames), 1);

% Undo/Redo history
data.undoStack = {};
data.redoStack = {};
data.maxUndoLevels = 50;

% Species ordering for z-order control
data.speciesOrder = (1:numel(speciesNames))';

% Search filter
data.searchFilter = "";
data.filteredIndices = (1:numel(speciesNames))';

% Sorting state: 'none', 'count_asc', 'count_desc', 'mass_asc', 'mass_desc', 'name_asc', 'name_desc'
data.sortMode = "none";

% Selected row in table
data.selectedRow = [];

% Axis labels visibility
data.showAxisLabels = true;
data.showAxisTicks = true;

% Scale cube settings
data.showScaleCube = false;
data.scaleCubeSize = 10;  % nm default
data.scaleCubeFig = [];
data.scaleCubeAx = [];
data.scaleCubePatches = [];
data.scaleCubeLabel = [];
data.viewSyncListener = [];
data.xlimSyncListener = [];
data.ylimSyncListener = [];
data.zlimSyncListener = [];
data.camViewAngleSyncListener = [];
data.projectionSyncListener = [];
data.positionSyncListener = [];
data.dataAspectSyncListener = [];

data.fn = struct( ...
    'computeGroups', @computeGroups, ...
    'mapColors', @mapColors, ...
    'createScatterHandles', @createScatterHandles, ...
    'updateScatter', @updateScatter, ...
    'updateTable', @updateTable, ...
    'buildTableData', @buildTableData, ...
    'initSampleFractions', @initSampleFractions, ...
    'initMarkerSizes', @initMarkerSizes, ...
    'rotateCoordsZ', @rotateCoordsZ, ...
    'pushUndo', @pushUndo, ...
    'updateStatusBar', @updateStatusBar);

tableData = buildTableData(data);

setupControls(controlFig, tableData, data);
configureLinkData(controlFig, baseVarName);

didRestore = tryRestoreStateFromAxis(controlFig);
if ~didRestore
    updateScatter(controlFig);
    applyBoundingBox(ax, data.bbox, data.fixBoundingBox);
end
persistStateToAxis(controlFig);

% Setup cleanup listeners for graceful handle deletion
setupCleanupListeners(controlFig, ax, fig);

info = struct();
info.speciesNames = speciesNames;
info.speciesCounts = speciesCounts;
info.colors = colors;
info.groupBy = mode;
end

%% Input Validation
function validatePosTable(pos)
    % Check for required columns
    requiredCols = {'x', 'y', 'z'};
    missingCols = setdiff(requiredCols, pos.Properties.VariableNames);
    if ~isempty(missingCols)
        error('scatterPlotPosWidget:missingColumns', ...
            'Missing required columns: %s', strjoin(missingCols, ', '));
    end

    % Check that x, y, z are numeric
    if ~isnumeric(pos.x) || ~isnumeric(pos.y) || ~isnumeric(pos.z)
        error('scatterPlotPosWidget:nonNumericCoords', ...
            'Columns x, y, z must be numeric.');
    end

    % Check for at least one grouping column
    hasIon = ismember('ion', pos.Properties.VariableNames);
    hasAtom = ismember('atom', pos.Properties.VariableNames);
    if ~hasIon && ~hasAtom
        error('scatterPlotPosWidget:missingGroupField', ...
            'pos must contain at least one of: ion, atom column.');
    end

    % Check for empty table
    if height(pos) == 0
        error('scatterPlotPosWidget:emptyTable', ...
            'pos table is empty.');
    end

    % Check for NaN/Inf in coordinates
    coordData = [pos.x, pos.y, pos.z];
    if any(isnan(coordData(:)))
        warning('scatterPlotPosWidget:nanCoords', ...
            'pos table contains NaN values in coordinates. These points will be excluded.');
    end
    if any(isinf(coordData(:)))
        warning('scatterPlotPosWidget:infCoords', ...
            'pos table contains Inf values in coordinates. These points will be excluded.');
    end
end

%% Cleanup Listeners
function setupCleanupListeners(controlFig, ax, scatterFig)
    % When scatter figure is closed, close control window and scale cube
    if isgraphics(scatterFig)
        addlistener(scatterFig, 'ObjectBeingDestroyed', @(~,~) cleanupAll(controlFig));
    end

    % When control window is closed, clean up scatter handles and scale cube
    if isgraphics(controlFig)
        addlistener(controlFig, 'ObjectBeingDestroyed', @(~,~) cleanupOnControlClose(controlFig, ax));
    end
end

function cleanupAll(controlFig)
    % Clean up scale cube (axes and listeners) first
    if isscalar(controlFig) && isgraphics(controlFig)
        persistStateToAxis(controlFig);
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if ~isempty(data)
            data = stopLinkData(data);
            data = destroyScaleCube(data);
            setappdata(controlFig, 'scatterPlotPosWidget', data);
        end
        % Then close control window
        delete(controlFig);
    end
end

function cleanupOnControlClose(controlFig, ~)
    persistStateToAxis(controlFig);
    % Clean up scale cube (axes and listeners)
    if isscalar(controlFig) && isgraphics(controlFig)
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if ~isempty(data)
            data = stopLinkData(data);
            data = destroyScaleCube(data);
            setappdata(controlFig, 'scatterPlotPosWidget', data);
        end
    end
    % Keep scatter handles in axes to allow reopening widget on same plot.
end

function safeClose(fig)
    if isValidGraphics(fig)
        delete(fig);
    end
end

function valid = isValidGraphics(h)
    % Safely check if h is a valid scalar graphics handle
    valid = ~isempty(h) && isscalar(h) && isgraphics(h);
end

function cleanupScatterHandles(ax)
    if isgraphics(ax)
        % Find and delete scatter handles created by this widget
        children = ax.Children;
        for i = 1:numel(children)
            if isa(children(i), 'matlab.graphics.chart.primitive.Scatter')
                delete(children(i));
            end
        end
    end
end

function cleanupWidgetScatterHandles(ax)
    if ~isgraphics(ax, 'axes')
        return;
    end
    widgetScatters = findobj(ax, 'Type', 'Scatter', 'Tag', 'scatterPlotPosWidgetScatter');
    if ~isempty(widgetScatters)
        delete(widgetScatters(ishghandle(widgetScatters)));
    end
end

function configureLinkData(controlFig, baseVarName)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    if ~data.linkToBase
        return;
    end
    if strlength(data.linkedVar) == 0
        if ~isempty(baseVarName)
            data.linkedVar = string(baseVarName);
        else
            warning('scatterPlotPosWidget:linkVarMissing', ...
                'linkToBase is true but no linkedVar name is available. Linking disabled.');
            data.linkToBase = false;
            setappdata(controlFig, 'scatterPlotPosWidget', data);
            return;
        end
    end
    if ~isvarname(char(data.linkedVar))
        warning('scatterPlotPosWidget:invalidLinkVar', ...
            'linkedVar "%s" is not a valid variable name. Linking disabled.', data.linkedVar);
        data.linkToBase = false;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        return;
    end

    data.linkSignature = computePosSignature(data.pos);
    if strlength(data.linkKey) == 0
        rawId = char(java.util.UUID.randomUUID);
        data.linkKey = "scatterPlotPosWidget_link_" + string(strrep(rawId, '-', '_'));
    end

    linkState = struct();
    linkState.controlFig = controlFig;
    linkState.linkedVar = char(data.linkedVar);
    linkState.lastSig = data.linkSignature;
    linkState.refreshFcn = @(newPos) refreshFromPos(controlFig, newPos);
    linkState.sigFcn = @(newPos) computePosSignature(newPos);
    linkState.updating = false;
    setappdata(0, char(data.linkKey), linkState);

    setappdata(controlFig, 'scatterPlotPosWidget', data);
    applyLinkDataSources(controlFig);
    startLinkRefreshTimer(controlFig);
    if isgraphics(data.scatterFig)
        linkdata(data.scatterFig, 'on');
        refreshdata(data.scatterFig);
    end
end

function data = stopLinkData(data)
    if isfield(data, 'scatterFig') && isgraphics(data.scatterFig)
        try
            linkdata(data.scatterFig, 'off');
        catch
        end
    end
    if isfield(data, 'linkKey') && strlength(data.linkKey) > 0
        try
            rmappdata(0, char(data.linkKey));
        catch
        end
        data.linkKey = "";
    end
    if isfield(data, 'linkRefreshTimer') && ~isempty(data.linkRefreshTimer) && isvalid(data.linkRefreshTimer)
        try
            stop(data.linkRefreshTimer);
        catch
        end
        try
            delete(data.linkRefreshTimer);
        catch
        end
    end
    data.linkRefreshTimer = [];
    data.linkUpdating = false;
end

function applyLinkDataSources(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~data.linkToBase
        return;
    end
    if strlength(data.linkKey) == 0 || strlength(data.linkedVar) == 0
        return;
    end
    varName = char(data.linkedVar);
    if ~isvarname(varName)
        return;
    end
    for i = 1:numel(data.scatterHandles)
        h = data.scatterHandles(i);
        if ~isgraphics(h)
            continue;
        end
        exprX = sprintf('scatterPlotPosWidget_linkedData(''%s'', %d, ''x'', %s)', ...
            char(data.linkKey), i, varName);
        exprY = sprintf('scatterPlotPosWidget_linkedData(''%s'', %d, ''y'', %s)', ...
            char(data.linkKey), i, varName);
        exprZ = sprintf('scatterPlotPosWidget_linkedData(''%s'', %d, ''z'', %s)', ...
            char(data.linkKey), i, varName);
        h.XDataSource = exprX;
        h.YDataSource = exprY;
        h.ZDataSource = exprZ;
    end
end

function startLinkRefreshTimer(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~data.linkToBase
        return;
    end
    if ~isfield(data, 'scatterFig') || ~isgraphics(data.scatterFig)
        return;
    end
    if isfield(data, 'linkRefreshTimer') && ~isempty(data.linkRefreshTimer) && isvalid(data.linkRefreshTimer)
        return;
    end
    data.linkRefreshTimer = timer( ...
        'ExecutionMode', 'fixedSpacing', ...
        'Period', 0.5, ...
        'TimerFcn', @(t,~) safeRefreshLinkData(controlFig, t));
    start(data.linkRefreshTimer);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
end

function safeRefreshLinkData(controlFig, timerObj)
    if isempty(controlFig) || ~isgraphics(controlFig)
        stopAndDeleteTimer(timerObj);
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~data.linkToBase
        stopAndDeleteTimer(timerObj);
        return;
    end
    if ~isfield(data, 'scatterFig') || ~isgraphics(data.scatterFig)
        stopAndDeleteTimer(timerObj);
        return;
    end
    try
        refreshdata(data.scatterFig);
    catch
    end
end

function stopAndDeleteTimer(timerObj)
    if isempty(timerObj) || ~isvalid(timerObj)
        return;
    end
    try
        stop(timerObj);
    catch
    end
    try
        delete(timerObj);
    catch
    end
end

function sig = computePosSignature(pos)
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

function didRestore = tryRestoreStateFromAxis(controlFig)
    didRestore = false;
    if ~isgraphics(controlFig)
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~isfield(data, 'restoreStateFromAxis') || ~data.restoreStateFromAxis
        return;
    end
    if ~isfield(data, 'ax') || ~isgraphics(data.ax, 'axes')
        return;
    end

    state = readStateFromAxis(data.ax);
    if isempty(state)
        return;
    end

    if isfield(state, 'posSignature')
        sigCurrent = computePosSignature(data.pos);
        if ~isequaln(state.posSignature, sigCurrent)
            return;
        end
    end

    try
        scatterPlotPosWidgetApplyState(controlFig, state);
        didRestore = true;
    catch ME
        warning('scatterPlotPosWidget:restoreStateFailed', ...
            'Could not restore saved axis state: %s', ME.message);
        didRestore = false;
    end
end

function persistStateToAxis(controlFig)
    if ~isgraphics(controlFig)
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    if ~isfield(data, 'persistStateToAxis') || ~data.persistStateToAxis
        return;
    end
    if isfield(data, 'suspendStatePersistence') && data.suspendStatePersistence
        return;
    end
    if ~isfield(data, 'ax') || ~isgraphics(data.ax, 'axes')
        return;
    end

    try
        state = buildWidgetState(data);
        writeStateToAxis(data.ax, state);
    catch
    end
end

function state = buildWidgetState(data)
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
    if isfield(data, 'bbox')
        state.boundingBox = data.bbox;
        state.fixBoundingBox = data.fixBoundingBox;
        state.clipToBoundingBox = data.clipToBoundingBox;
        state.bboxUseSpan = data.bboxUseSpan;
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

    state.posSignature = computePosSignature(data.pos);

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

function restoreStatePersistence(controlFig, originalSuspend)
    if ~isgraphics(controlFig)
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data.suspendStatePersistence = originalSuspend;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
end

function refreshFromPos(controlFig, newPos)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end

    validatePosTable(newPos);
    data.pos = newPos;
    data.coordsOriginal = [newPos.x, newPos.y, newPos.z];
    data.coords = rotateCoordsZ(data.coordsOriginal, data.rotationAngle);
    data.hasIsotope = ismember('isotope', newPos.Properties.VariableNames);
    data.hasCharge = ismember('chargeState', newPos.Properties.VariableNames);

    if ~data.bboxUserEdited
        data.bboxOriginal = computeBoundingBox(data.coords, data.bboxPadding);
        data.bbox = data.bboxOriginal;
    end

    [speciesNames, baseNames, displayNames, speciesIndices, speciesCounts] = computeGroups( ...
        data.pos, data.mode, data.splitIsotope, data.splitCharge, data.showUnranged);
    colors = mapColors(baseNames, data.colorScheme);
    markerSizes = initMarkerSizes(numel(speciesNames), data.markerSizeGlobal);

    allowRebuild = ~(isfield(data, 'linkUpdating') && data.linkUpdating);
    reuseHandles = numel(speciesNames) == numel(data.speciesNames) && ...
        all(string(displayNames) == string(data.displayNames));

    if ~allowRebuild && ~reuseHandles
        oldSpeciesNames = data.speciesNames;
        oldDisplayNames = data.displayNames;
        oldBaseNames = data.baseNames;
        oldColors = data.colors;
        oldMarkerSizes = data.markerSizes;

        [isMatch, idx] = ismember(oldSpeciesNames, speciesNames);
        alignedIndices = cell(numel(oldSpeciesNames), 1);
        alignedCounts = zeros(numel(oldSpeciesNames), 1);
        for i = 1:numel(oldSpeciesNames)
            if isMatch(i)
                alignedIndices{i} = speciesIndices{idx(i)};
                alignedCounts(i) = speciesCounts(idx(i));
            else
                alignedIndices{i} = [];
                alignedCounts(i) = 0;
            end
        end

        speciesNames = oldSpeciesNames;
        baseNames = oldBaseNames;
        displayNames = oldDisplayNames;
        colors = oldColors;
        markerSizes = oldMarkerSizes;
        speciesIndices = alignedIndices;
        speciesCounts = alignedCounts;

        visible = data.visible;
        sampleFractions = data.sampleFractions;
    else
        visible = true(numel(speciesNames), 1);
        sampleFractions = initSampleFractions(speciesCounts, data.sampleValue);

        [isMatch, idx] = ismember(speciesNames, data.speciesNames);
        if any(isMatch)
            visible(isMatch) = data.visible(idx(isMatch));
            sampleFractions(isMatch) = data.sampleFractions(idx(isMatch));
            if numel(data.markerSizes) == numel(data.speciesNames)
                markerSizes(isMatch) = data.markerSizes(idx(isMatch));
            end
        end
    end

    if allowRebuild && ~reuseHandles
        if ~isempty(data.scatterHandles)
            delete(data.scatterHandles(ishghandle(data.scatterHandles)));
        end
        axes(data.ax);
        holdState = ishold(data.ax);
        hold(data.ax, 'on');
        data.scatterHandles = createScatterHandles(data.ax, displayNames, colors, markerSizes);
        if ~holdState
            hold(data.ax, 'off');
        end
    elseif reuseHandles
        for i = 1:numel(data.scatterHandles)
            h = data.scatterHandles(i);
            if ~isgraphics(h)
                continue;
            end
            set(h, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), ...
                'DisplayName', char(displayNames(i)));
        end
    end

    data.speciesNames = speciesNames;
    data.baseNames = baseNames;
    data.displayNames = displayNames;
    data.speciesIndices = speciesIndices;
    data.speciesCounts = speciesCounts;
    data.colors = colors;
    data.visible = visible;
    data.sampleFractions = sampleFractions;
    data.markerSizes = markerSizes;
    data.linkSignature = computePosSignature(newPos);

    data.samplingCache = cell(numel(speciesNames), 1);
    data.samplingCacheValid = false(numel(speciesNames), 1);

    setappdata(controlFig, 'scatterPlotPosWidget', data);
    if data.linkToBase && strlength(data.linkKey) > 0 && isappdata(0, char(data.linkKey))
        linkState = getappdata(0, char(data.linkKey));
        linkState.lastSig = data.linkSignature;
        linkState.controlFig = controlFig;
        setappdata(0, char(data.linkKey), linkState);
    end
    applyLinkDataSources(controlFig);
    updateTable(controlFig);
    updateScatter(controlFig);
    updateBoundingBoxControls(controlFig);
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
    posIdx = (1:height(pos))';

    if mode == "ionic"
        if ismember('ion', pos.Properties.VariableNames)
            baseRaw = string(pos.ion);
        elseif ismember('atom', pos.Properties.VariableNames)
            baseRaw = string(pos.atom);
        else
            error('scatterPlotPosWidget:missingGroupField', ...
                'pos must contain ion or atom column.');
        end
    else
        if ismember('atom', pos.Properties.VariableNames)
            baseRaw = string(pos.atom);
        elseif ismember('ion', pos.Properties.VariableNames)
            baseRaw = string(pos.ion);
        else
            error('scatterPlotPosWidget:missingGroupField', ...
                'pos must contain ion or atom column.');
        end
    end

    if showUnranged
        baseRaw(ismissing(baseRaw)) = "unranged";
    else
        keep = ~ismissing(baseRaw);
        baseRaw = baseRaw(keep);
        posIdx = posIdx(keep);
    end

    % In ionic mode, use a charge-neutral base name and let splitCharge
    % explicitly control whether charge states are merged or separated.
    if mode == "ionic"
        baseCore = stripChargeFromName(baseRaw);
    else
        baseCore = baseRaw;
    end

    groupKey = baseCore;
    if splitIsotope && ismember('isotope', pos.Properties.VariableNames)
        iso = string(pos.isotope(posIdx));
        iso(ismissing(iso)) = "NaN";
        groupKey = groupKey + "-" + iso;
    end
    if splitCharge && ismember('chargeState', pos.Properties.VariableNames)
        cs = pos.chargeState(posIdx);
        csString = strings(size(cs));
        csString(:) = "";
        valid = ~isnan(cs);
        for k = 1:numel(cs)
            if valid(k)
                n = round(cs(k));
                if n <= 0
                    csString(k) = "";
                else
                    csString(k) = string(repmat('+', 1, n));
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
        groupIndices{i} = posIdx(idx);
        groupCounts(i) = numel(idx);
        baseNames(i) = baseCore(idx(1));
        displayNames(i) = formatDisplayName(baseCore(idx(1)), pos, posIdx(idx(1)), ...
            splitIsotope, splitCharge, showUnranged);
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
        'Color', [0.94 0.95 0.96], ...
        'Units', 'pixels', ...
        'Position', [100 80 520 900], ...
        'Resize', 'on', ...
        'DefaultUicontrolFontName', 'Helvetica', ...
        'DefaultUicontrolFontSize', 10, ...
        'KeyPressFcn', @onKeyPress);
end

function setupControls(controlFig, tableData, data)
    setappdata(controlFig, 'scatterPlotPosWidget', data);

    % UI style tokens for cleaner layout and visual consistency
    figBg = get(controlFig, 'Color');
    textBg = figBg;
    panelBg = [0.98 0.985 0.99];
    btnBg = [1 1 1];
    accentColor = [0.16 0.33 0.53];

    % Standardized geometry
    btnH = 0.030;
    lblH = 0.024;
    gap = 0.005;
    margin = 0.035;

    % Title stripe for orientation
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Grouping and Sampling Controls', ...
        'Units', 'normalized', ...
        'Position', [margin 0.975 0.92 0.02], ...
        'HorizontalAlignment', 'left', ...
        'FontWeight', 'bold', ...
        'ForegroundColor', accentColor, ...
        'BackgroundColor', textBg);

    % Mode selection panel
    modeGroup = uibuttongroup(controlFig, ...
        'Units', 'normalized', ...
        'Position', [margin 0.925 0.92 0.052], ...
        'Title', 'Display Mode', ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold', ...
        'SelectionChangedFcn', @(src, evd) onGroupingChanged(controlFig, evd, src));

    rbIonic = uicontrol(modeGroup, 'Style', 'radiobutton', ...
        'String', 'Ionic', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.15 0.35 0.7], ...
        'BackgroundColor', panelBg, ...
        'TooltipString', 'Group atoms by ionic species');

    rbAtomic = uicontrol(modeGroup, 'Style', 'radiobutton', ...
        'String', 'Atomic', ...
        'Units', 'normalized', ...
        'Position', [0.38 0.15 0.35 0.7], ...
        'BackgroundColor', panelBg, ...
        'TooltipString', 'Group atoms by element');

    if data.mode == "atomic"
        modeGroup.SelectedObject = rbAtomic;
    else
        modeGroup.SelectedObject = rbIonic;
    end

    % Split options row
    y = 0.905;
    splitIsotope = uicontrol(controlFig, 'Style', 'checkbox', ...
        'String', 'Split isotopes', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.28 lblH], ...
        'Value', data.splitIsotope, ...
        'BackgroundColor', textBg, ...
        'Enable', boolToOnOff(data.hasIsotope), ...
        'TooltipString', 'Separate species by isotope mass number', ...
        'Callback', @(~, ~) onGroupingChanged(controlFig));

    splitCharge = uicontrol(controlFig, 'Style', 'checkbox', ...
        'String', 'Split charge state', ...
        'Units', 'normalized', ...
        'Position', [0.33 y 0.28 lblH], ...
        'Value', data.splitCharge, ...
        'BackgroundColor', textBg, ...
        'Enable', boolToOnOff(data.hasCharge), ...
        'TooltipString', 'Separate species by charge state', ...
        'Callback', @(~, ~) onGroupingChanged(controlFig));

    % Live update and undo/redo row
    y = y - lblH - gap;
    liveUpdateCb = uicontrol(controlFig, 'Style', 'checkbox', ...
        'String', 'Live Update', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.22 lblH], ...
        'Value', data.liveUpdate, ...
        'BackgroundColor', textBg, ...
        'TooltipString', 'Update scatter plot immediately on changes. Uncheck to defer updates.', ...
        'Callback', @(src, ~) onLiveUpdateToggle(src, controlFig));

    applyBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'Apply', ...
        'Units', 'normalized', ...
        'Position', [0.27 y 0.18 btnH], ...
        'Enable', 'off', ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Apply pending changes to the scatter plot', ...
        'Callback', @(~, ~) onApplyChanges(controlFig));

    undoBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'Undo', ...
        'Units', 'normalized', ...
        'Position', [0.46 y 0.14 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Undo last change (Ctrl+Z)', ...
        'Callback', @(~, ~) onUndo(controlFig));

    redoBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'Redo', ...
        'Units', 'normalized', ...
        'Position', [0.61 y 0.14 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Redo last undone change (Ctrl+Y)', ...
        'Callback', @(~, ~) onRedo(controlFig));

    helpBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', '?', ...
        'Units', 'normalized', ...
        'Position', [0.76 y 0.08 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Show keyboard shortcuts and help', ...
        'Callback', @(~, ~) showKeyboardHelp());

    % Sample fraction
    y = y - lblH - gap;
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Sample fraction:', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.25 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    y = y - btnH - gap*0.7;
    slider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 0, 'Max', 1, 'Value', min(data.sampleValue, 1), ...
        'Units', 'normalized', ...
        'Position', [margin y 0.68 btnH], ...
        'TooltipString', 'Fraction of points to display (0-1)', ...
        'Callback', @(src, ~) onSlider(src, controlFig));

    editBox = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.sampleValue, '%.3f'), ...
        'Units', 'normalized', ...
        'Position', [0.73 y 0.23 btnH], ...
        'TooltipString', 'Enter exact sample fraction or absolute count', ...
        'Callback', @(src, ~) onEdit(src, controlFig));

    % Marker size
    y = y - lblH - gap*0.8;
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Marker size:', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.25 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    y = y - btnH - gap*0.7;
    markerMax = max(50, data.markerSizeGlobal);
    markerSlider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 1, 'Max', markerMax, 'Value', min(data.markerSizeGlobal, markerMax), ...
        'Units', 'normalized', ...
        'Position', [margin y 0.68 btnH], ...
        'TooltipString', 'Size of scatter plot markers', ...
        'Callback', @(src, ~) onMarkerSizeSlider(src, controlFig));

    markerEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.markerSizeGlobal, '%.1f'), ...
        'Units', 'normalized', ...
        'Position', [0.73 y 0.23 btnH], ...
        'TooltipString', 'Enter exact marker size', ...
        'Callback', @(src, ~) onMarkerSizeEdit(src, controlFig));

    % Random seed and Rotation row
    y = y - lblH - gap*0.8;
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Random seed:', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.22 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg);

    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Rotation (deg):', ...
        'Units', 'normalized', ...
        'Position', [0.50 y 0.25 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg);

    y = y - btnH - gap*0.7;
    seedEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.randomSeed, '%.0f'), ...
        'Units', 'normalized', ...
        'Position', [margin y 0.20 btnH], ...
        'TooltipString', 'Random seed for reproducible sampling', ...
        'Callback', @(src, ~) onRandomSeedEdit(src, controlFig));

    rotationSlider = uicontrol(controlFig, 'Style', 'slider', ...
        'Min', 0, 'Max', 360, 'Value', data.rotationAngle, ...
        'Units', 'normalized', ...
        'Position', [0.50 y 0.35 btnH], ...
        'TooltipString', 'Rotate data around Z-axis', ...
        'Callback', @(src, ~) onRotationSlider(src, controlFig));

    rotationEdit = uicontrol(controlFig, 'Style', 'edit', ...
        'String', num2str(data.rotationAngle, '%.1f'), ...
        'Units', 'normalized', ...
        'Position', [0.86 y 0.10 btnH], ...
        'TooltipString', 'Enter exact rotation angle in degrees', ...
        'Callback', @(src, ~) onRotationEdit(src, controlFig));

    % Bounding box panel
    y = y - 0.12 - gap;
    bboxPanel = uipanel(controlFig, 'Title', 'Bounding Box', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.92 0.11], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    bboxFix = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Fix', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.70 0.14 0.25], ...
        'BackgroundColor', panelBg, ...
        'TooltipString', 'Lock axis limits to bounding box', ...
        'Value', data.fixBoundingBox, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'toggle'));

    bboxClip = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Clip', ...
        'Units', 'normalized', ...
        'Position', [0.16 0.70 0.14 0.25], ...
        'BackgroundColor', panelBg, ...
        'TooltipString', 'Only show points inside bounding box', ...
        'Value', data.clipToBoundingBox, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'clip'));

    bboxSpan = uicontrol(bboxPanel, 'Style', 'checkbox', ...
        'String', 'Span', ...
        'Units', 'normalized', ...
        'Position', [0.30 0.70 0.16 0.25], ...
        'BackgroundColor', panelBg, ...
        'TooltipString', 'Edit as center/span instead of min/max', ...
        'Value', data.bboxUseSpan, ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'span'));

    bboxUseCurrent = uicontrol(bboxPanel, 'Style', 'pushbutton', ...
        'String', 'Use Current', ...
        'Units', 'normalized', ...
        'Position', [0.47 0.70 0.25 0.25], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Set bounding box from current view limits', ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'current'));

    bboxReset = uicontrol(bboxPanel, 'Style', 'pushbutton', ...
        'String', 'Reset', ...
        'Units', 'normalized', ...
        'Position', [0.73 0.70 0.25 0.25], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Reset bounding box to original data extent', ...
        'Callback', @(~, ~) onBoundingBoxChanged(controlFig, 'reset'));

    bboxEdits = struct();
    [bboxEdits.xmin, bboxEdits.xmax] = makeBBoxRow(bboxPanel, 'X', 0.46, data.bbox.xmin, data.bbox.xmax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));
    [bboxEdits.ymin, bboxEdits.ymax] = makeBBoxRow(bboxPanel, 'Y', 0.24, data.bbox.ymin, data.bbox.ymax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));
    [bboxEdits.zmin, bboxEdits.zmax] = makeBBoxRow(bboxPanel, 'Z', 0.02, data.bbox.zmin, data.bbox.zmax, ...
        @(~, ~) onBoundingBoxChanged(controlFig, 'edit'));

    % Display Options panel (axis labels, scale cube)
    y = y - 0.055 - gap;
    displayPanel = uipanel(controlFig, 'Title', 'Display Options', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.92 0.05], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    showAxisLabels = uicontrol(displayPanel, 'Style', 'checkbox', ...
        'String', 'Axis Labels', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.15 0.22 0.7], ...
        'BackgroundColor', panelBg, ...
        'Value', data.showAxisLabels, ...
        'TooltipString', 'Show/hide X, Y, Z axis labels', ...
        'Callback', @(src, ~) onAxisLabelsToggle(src, controlFig));

    showAxisTicks = uicontrol(displayPanel, 'Style', 'checkbox', ...
        'String', 'Ticks', ...
        'Units', 'normalized', ...
        'Position', [0.25 0.15 0.15 0.7], ...
        'BackgroundColor', panelBg, ...
        'Value', data.showAxisTicks, ...
        'TooltipString', 'Show/hide axis tick marks and numbers', ...
        'Callback', @(src, ~) onAxisTicksToggle(src, controlFig));

    showScaleCube = uicontrol(displayPanel, 'Style', 'checkbox', ...
        'String', 'Scale Cube', ...
        'Units', 'normalized', ...
        'Position', [0.42 0.15 0.22 0.7], ...
        'BackgroundColor', panelBg, ...
        'Value', data.showScaleCube, ...
        'TooltipString', 'Show RGB scale cube (X=Red, Y=Green, Z=Blue)', ...
        'Callback', @(src, ~) onScaleCubeToggle(src, controlFig));

    uicontrol(displayPanel, 'Style', 'text', ...
        'String', 'Size:', ...
        'Units', 'normalized', ...
        'Position', [0.65 0.15 0.10 0.7], ...
        'HorizontalAlignment', 'right', ...
        'BackgroundColor', panelBg);

    scaleCubeSizeEdit = uicontrol(displayPanel, 'Style', 'edit', ...
        'String', num2str(data.scaleCubeSize), ...
        'Units', 'normalized', ...
        'Position', [0.76 0.15 0.12 0.7], ...
        'TooltipString', 'Scale cube edge length (nm)', ...
        'Callback', @(src, ~) onScaleCubeSizeEdit(src, controlFig));

    uicontrol(displayPanel, 'Style', 'text', ...
        'String', 'nm', ...
        'Units', 'normalized', ...
        'Position', [0.89 0.15 0.10 0.7], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', panelBg);

    % Export panel
    y = y - 0.055 - gap;
    exportPanel = uipanel(controlFig, 'Title', 'Export', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.92 0.048], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    exportImagesBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Images', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Export individual species images to files', ...
        'Callback', @(~, ~) onExportImages(controlFig));

    exportTurntableBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Turntable', ...
        'Units', 'normalized', ...
        'Position', [0.26 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Export rotating turntable animation', ...
        'Callback', @(~, ~) onExportTurntable(controlFig));

    exportProfileWsBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Save WS', ...
        'Units', 'normalized', ...
        'Position', [0.50 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Export visualisation profile to MATLAB workspace', ...
        'Callback', @(~, ~) onExportProfileToWorkspace(controlFig));

    importProfileWsBtn = uicontrol(exportPanel, 'Style', 'pushbutton', ...
        'String', 'Load WS', ...
        'Units', 'normalized', ...
        'Position', [0.74 0.10 0.24 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Load a visualisation profile/state from MATLAB workspace', ...
        'Callback', @(~, ~) onImportProfileFromWorkspace(controlFig));

    % Visibility panel
    y = y - 0.055 - gap;
    visibilityPanel = uipanel(controlFig, 'Title', 'Visibility', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.92 0.048], ...
        'BackgroundColor', panelBg, ...
        'ForegroundColor', accentColor, ...
        'FontWeight', 'bold');

    showAllBtn = uicontrol(visibilityPanel, 'Style', 'pushbutton', ...
        'String', 'Show All', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Show all species (keyboard: A)', ...
        'Callback', @(~, ~) onShowAll(controlFig));

    hideAllBtn = uicontrol(visibilityPanel, 'Style', 'pushbutton', ...
        'String', 'Hide All', ...
        'Units', 'normalized', ...
        'Position', [0.26 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Hide all species (keyboard: H)', ...
        'Callback', @(~, ~) onHideAll(controlFig));

    resetViewBtn = uicontrol(visibilityPanel, 'Style', 'pushbutton', ...
        'String', 'Reset View', ...
        'Units', 'normalized', ...
        'Position', [0.50 0.10 0.23 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Reset 3D view to default (keyboard: R)', ...
        'Callback', @(~, ~) onResetView(controlFig));

    colorBtn = uicontrol(visibilityPanel, 'Style', 'pushbutton', ...
        'String', 'Color...', ...
        'Units', 'normalized', ...
        'Position', [0.74 0.10 0.24 0.80], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Change color of selected species', ...
        'Callback', @(~, ~) onChangeColor(controlFig));

    % Search and sort row
    y = y - btnH - gap;
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Search:', ...
        'Units', 'normalized', ...
        'Position', [margin y 0.12 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg);

    searchBox = uicontrol(controlFig, 'Style', 'edit', ...
        'String', '', ...
        'Units', 'normalized', ...
        'Position', [margin + 0.12 y 0.28 btnH], ...
        'TooltipString', 'Filter species by name (case-insensitive)', ...
        'Callback', @(src, ~) onSearchChanged(src, controlFig));

    clearSearchBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', 'X', ...
        'Units', 'normalized', ...
        'Position', [margin + 0.41 y 0.06 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Clear search filter', ...
        'Callback', @(~, ~) onClearSearch(controlFig));

    % Sorting dropdown
    uicontrol(controlFig, 'Style', 'text', ...
        'String', 'Sort:', ...
        'Units', 'normalized', ...
        'Position', [margin + 0.49 y 0.08 lblH], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', textBg);

    sortDropdown = uicontrol(controlFig, 'Style', 'popupmenu', ...
        'String', {'None', 'Name A-Z', 'Name Z-A', 'Count (High)', 'Count (Low)', 'Mass (High)', 'Mass (Low)'}, ...
        'Units', 'normalized', ...
        'Position', [margin + 0.57 y 0.26 btnH], ...
        'TooltipString', 'Sort species table by name, count, or atomic mass', ...
        'Callback', @(src, ~) onSortChanged(src, controlFig));

    % Move up/down buttons for z-order
    moveUpBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', char(9650), ...
        'Units', 'normalized', ...
        'Position', [margin + 0.84 y 0.06 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Move selected species up in draw order (appears on top)', ...
        'Callback', @(~, ~) onMoveSpecies(controlFig, -1));

    moveDownBtn = uicontrol(controlFig, 'Style', 'pushbutton', ...
        'String', char(9660), ...
        'Units', 'normalized', ...
        'Position', [margin + 0.91 y 0.06 btnH], ...
        'BackgroundColor', btnBg, ...
        'TooltipString', 'Move selected species down in draw order (appears below)', ...
        'Callback', @(~, ~) onMoveSpecies(controlFig, 1));

    % Species table
    % Keep the top aligned with the current flow layout and expand downward
    % to use the available space above the status bar.
    y = y - 0.20 - gap;
    tableTop = y + 0.20;
    tableBottom = 0.035; % leave a small gap above the status bar
    tableHeight = max(0.20, tableTop - tableBottom);
    tbl = uitable(controlFig, ...
        'Data', tableData, ...
        'ColumnName', {'Show', 'Color', 'Species', 'Count', '# Shown', 'Frac', 'Size', 'Link'}, ...
        'RowName', {}, ...
        'ColumnEditable', [true false false false true true true true], ...
        'ColumnWidth', {35, 35, 'auto', 55, 55, 50, 45, 40}, ...
        'BackgroundColor', [1 1 1; 0.96 0.97 0.99], ...
        'FontName', 'Helvetica', ...
        'FontSize', 10, ...
        'Units', 'normalized', ...
        'Position', [margin tableBottom 0.92 tableHeight], ...
        'TooltipString', 'Species table - toggle visibility, adjust sample fractions and marker sizes', ...
        'CellEditCallback', @(src, evd) onTableEdit(src, evd, controlFig), ...
        'CellSelectionCallback', @(src, evd) onTableSelection(src, evd, controlFig));

    % Status bar at the bottom
    statusBar = uicontrol(controlFig, 'Style', 'text', ...
        'String', '', ...
        'Units', 'normalized', ...
        'Position', [margin 0.005 0.92 0.025], ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 9, ...
        'ForegroundColor', [0.1 0.1 0.1], ...
        'BackgroundColor', [0.90 0.92 0.95]);

    setappdata(controlFig, 'scatterPlotPosWidgetControls', struct( ...
        'slider', slider, 'editBox', editBox, 'table', tbl, ...
        'markerSlider', markerSlider, 'markerEdit', markerEdit, ...
        'seedEdit', seedEdit, ...
        'rotationSlider', rotationSlider, 'rotationEdit', rotationEdit, ...
        'exportPanel', exportPanel, ...
        'exportImagesBtn', exportImagesBtn, 'exportTurntableBtn', exportTurntableBtn, ...
        'modeGroup', modeGroup, 'splitIsotope', splitIsotope, 'splitCharge', splitCharge, ...
        'bboxPanel', bboxPanel, 'bboxFix', bboxFix, 'bboxClip', bboxClip, 'bboxSpan', bboxSpan, ...
        'bboxUseCurrent', bboxUseCurrent, 'bboxReset', bboxReset, 'bboxEdits', bboxEdits, ...
        'displayPanel', displayPanel, 'showAxisLabels', showAxisLabels, 'showAxisTicks', showAxisTicks, ...
        'showScaleCube', showScaleCube, 'scaleCubeSizeEdit', scaleCubeSizeEdit, ...
        'visibilityPanel', visibilityPanel, ...
        'liveUpdateCb', liveUpdateCb, 'applyBtn', applyBtn, ...
        'undoBtn', undoBtn, 'redoBtn', redoBtn, ...
        'searchBox', searchBox, 'clearSearchBtn', clearSearchBtn, ...
        'sortDropdown', sortDropdown, ...
        'moveUpBtn', moveUpBtn, 'moveDownBtn', moveDownBtn, ...
        'showAllBtn', showAllBtn, 'hideAllBtn', hideAllBtn, ...
        'colorBtn', colorBtn, 'resetViewBtn', resetViewBtn, ...
        'exportProfileWsBtn', exportProfileWsBtn, 'importProfileWsBtn', importProfileWsBtn, ...
        'statusBar', statusBar, 'helpBtn', helpBtn));

    % Keep controls fixed-size on resize and let only the species table expand.
    configureControlResizeBehavior(controlFig);

    updateStatusBar(controlFig);
    updateUndoRedoButtons(controlFig);
    drawnow;
end

function configureControlResizeBehavior(controlFig)
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(ctrls) || ~isfield(ctrls, 'table') || ~isfield(ctrls, 'statusBar')
        return;
    end

    figPos = getpixelposition(controlFig);
    figW = figPos(3);
    figH = figPos(4);

    topLevel = findall(controlFig, '-depth', 1);
    topLevel = topLevel(topLevel ~= controlFig);
    keepDynamic = [ctrls.table, ctrls.statusBar];
    fixedHandles = setdiff(topLevel, keepDynamic);
    fixedHandles = fixedHandles(arrayfun(@isgraphics, fixedHandles));

    nFixed = numel(fixedHandles);
    fixedLeft = zeros(nFixed, 1);
    fixedTop = zeros(nFixed, 1);
    fixedWidth = zeros(nFixed, 1);
    fixedHeight = zeros(nFixed, 1);

    for i = 1:nFixed
        h = fixedHandles(i);
        set(h, 'Units', 'pixels');
        p = getpixelposition(h);
        fixedLeft(i) = p(1);
        fixedTop(i) = figH - (p(2) + p(4));
        fixedWidth(i) = p(3);
        fixedHeight(i) = p(4);
    end

    set(ctrls.table, 'Units', 'pixels');
    set(ctrls.statusBar, 'Units', 'pixels');

    tablePos = getpixelposition(ctrls.table);
    statusPos = getpixelposition(ctrls.statusBar);

    layout = struct();
    layout.fixedHandles = fixedHandles;
    layout.fixedLeft = fixedLeft;
    layout.fixedTop = fixedTop;
    layout.fixedWidth = fixedWidth;
    layout.fixedHeight = fixedHeight;
    layout.tableHandle = ctrls.table;
    layout.tableLeft = tablePos(1);
    layout.tableRightMargin = figW - (tablePos(1) + tablePos(3));
    layout.tableTopOffset = figH - (tablePos(2) + tablePos(4));
    layout.tableBottom = tablePos(2);
    layout.tableMinHeight = max(140, min(220, tablePos(4)));
    layout.statusHandle = ctrls.statusBar;
    layout.statusLeft = statusPos(1);
    layout.statusRightMargin = figW - (statusPos(1) + statusPos(3));
    layout.statusBottom = statusPos(2);
    layout.statusHeight = statusPos(4);

    setappdata(controlFig, 'scatterPlotPosWidgetLayout', layout);
    set(controlFig, 'SizeChangedFcn', @(src, ~) onControlWindowResize(src));
    onControlWindowResize(controlFig);
end

function onControlWindowResize(controlFig)
    layout = getappdata(controlFig, 'scatterPlotPosWidgetLayout');
    if isempty(layout) || ~isgraphics(controlFig)
        return;
    end

    figPos = getpixelposition(controlFig);
    figW = figPos(3);
    figH = figPos(4);

    % Reposition fixed controls with top anchoring; keep size unchanged.
    nFixed = numel(layout.fixedHandles);
    for i = 1:nFixed
        h = layout.fixedHandles(i);
        if ~isgraphics(h)
            continue;
        end
        newX = layout.fixedLeft(i);
        newY = figH - layout.fixedTop(i) - layout.fixedHeight(i);
        setpixelposition(h, [newX, newY, layout.fixedWidth(i), layout.fixedHeight(i)]);
    end

    % Expand/shrink only the table.
    if isgraphics(layout.tableHandle)
        tableW = max(180, figW - layout.tableLeft - layout.tableRightMargin);
        tableTop = figH - layout.tableTopOffset;
        tableH = max(layout.tableMinHeight, tableTop - layout.tableBottom);
        setpixelposition(layout.tableHandle, [layout.tableLeft, layout.tableBottom, tableW, tableH]);
    end

    % Keep status bar width matched to window width.
    if isgraphics(layout.statusHandle)
        statusW = max(120, figW - layout.statusLeft - layout.statusRightMargin);
        setpixelposition(layout.statusHandle, [layout.statusLeft, layout.statusBottom, statusW, layout.statusHeight]);
    end
end

%% Keyboard Shortcuts
function onKeyPress(src, evd)
    try
        controlFig = src;
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        % Handle modifier keys
        isCtrl = any(strcmp(evd.Modifier, 'control')) || any(strcmp(evd.Modifier, 'command'));

        switch evd.Key
            case 'space'
                % Toggle selected species
                if ~isempty(data.selectedRow) && data.selectedRow <= numel(data.visible)
                    pushUndo(controlFig);
                    idx = data.filteredIndices(data.selectedRow);
                    data.visible(idx) = ~data.visible(idx);
                    setappdata(controlFig, 'scatterPlotPosWidget', data);
                    conditionalUpdate(controlFig);
                end
            case 'a'
                if ~isCtrl
                    onShowAll(controlFig);
                end
            case 'h'
                onHideAll(controlFig);
            case 'r'
                onResetView(controlFig);
            case 'z'
                if isCtrl
                    onUndo(controlFig);
                end
            case 'y'
                if isCtrl
                    onRedo(controlFig);
                end
            case {'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'}
                % Toggle species 1-20 based on current table order (after filtering/sorting)
                % Keys 1-9 map to species 1-9, key 0 maps to species 10
                % Shift+1-9 map to species 11-19, Shift+0 maps to species 20
                if evd.Key == '0'
                    tableRow = 10;
                else
                    tableRow = str2double(evd.Key);
                end
                % Add 10 if Shift is held
                if ismember('shift', evd.Modifier)
                    tableRow = tableRow + 10;
                end
                if tableRow <= numel(data.filteredIndices)
                    pushUndo(controlFig);
                    actualIdx = data.filteredIndices(tableRow);
                    data.visible(actualIdx) = ~data.visible(actualIdx);
                    setappdata(controlFig, 'scatterPlotPosWidget', data);
                    conditionalUpdate(controlFig);
                end
        end
    catch ME
        warning('scatterPlotPosWidget:keyPressError', 'Key press error: %s', ME.message);
    end
end

function showKeyboardHelp()
    helpText = sprintf([...
        'Keyboard Shortcuts:\n\n' ...
        'Space  - Toggle visibility of selected species\n' ...
        'A      - Show all species\n' ...
        'H      - Hide all species\n' ...
        'R      - Reset view to default\n' ...
        '1-9    - Toggle visibility of species 1-9\n' ...
        'Ctrl+Z - Undo\n' ...
        'Ctrl+Y - Redo\n']);
    msgbox(helpText, 'Keyboard Shortcuts', 'help');
end

%% Undo/Redo Functions
function pushUndo(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data = ensureMarkerSizeLinkState(data);

    % Create state snapshot
    state = struct();
    state.visible = data.visible;
    state.sampleFractions = data.sampleFractions;
    state.markerSizes = data.markerSizes;
    state.markerSizeLinked = data.markerSizeLinked;
    state.colors = data.colors;
    state.speciesOrder = data.speciesOrder;

    % Push to undo stack
    data.undoStack{end+1} = state;

    % Limit stack size
    if numel(data.undoStack) > data.maxUndoLevels
        data.undoStack = data.undoStack(2:end);
    end

    % Clear redo stack on new action
    data.redoStack = {};

    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateUndoRedoButtons(controlFig);
end

function onUndo(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || isempty(data.undoStack)
            return;
        end
        data = ensureMarkerSizeLinkState(data);

        % Save current state to redo stack
        currentState = struct();
        currentState.visible = data.visible;
        currentState.sampleFractions = data.sampleFractions;
        currentState.markerSizes = data.markerSizes;
        currentState.markerSizeLinked = data.markerSizeLinked;
        currentState.colors = data.colors;
        currentState.speciesOrder = data.speciesOrder;
        data.redoStack{end+1} = currentState;

        % Pop from undo stack
        state = data.undoStack{end};
        data.undoStack = data.undoStack(1:end-1);

        % Restore state
        data.visible = state.visible;
        data.sampleFractions = state.sampleFractions;
        data.markerSizes = state.markerSizes;
        if isfield(state, 'markerSizeLinked')
            data.markerSizeLinked = logical(state.markerSizeLinked);
        else
            data.markerSizeLinked = true(numel(data.speciesNames), 1);
        end
        data.colors = state.colors;
        if isfield(state, 'speciesOrder')
            data.speciesOrder = state.speciesOrder;
        end

        % Update scatter colors
        for i = 1:numel(data.scatterHandles)
            if isvalid(data.scatterHandles(i))
                set(data.scatterHandles(i), 'MarkerFaceColor', data.colors(i,:), ...
                    'MarkerEdgeColor', data.colors(i,:));
            end
        end

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
        updateScatter(controlFig);
        updateUndoRedoButtons(controlFig);
    catch ME
        warning('scatterPlotPosWidget:undoError', 'Undo error: %s', ME.message);
    end
end

function onRedo(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || isempty(data.redoStack)
            return;
        end
        data = ensureMarkerSizeLinkState(data);

        % Save current state to undo stack
        currentState = struct();
        currentState.visible = data.visible;
        currentState.sampleFractions = data.sampleFractions;
        currentState.markerSizes = data.markerSizes;
        currentState.markerSizeLinked = data.markerSizeLinked;
        currentState.colors = data.colors;
        currentState.speciesOrder = data.speciesOrder;
        data.undoStack{end+1} = currentState;

        % Pop from redo stack
        state = data.redoStack{end};
        data.redoStack = data.redoStack(1:end-1);

        % Restore state
        data.visible = state.visible;
        data.sampleFractions = state.sampleFractions;
        data.markerSizes = state.markerSizes;
        if isfield(state, 'markerSizeLinked')
            data.markerSizeLinked = logical(state.markerSizeLinked);
        else
            data.markerSizeLinked = true(numel(data.speciesNames), 1);
        end
        data.colors = state.colors;
        if isfield(state, 'speciesOrder')
            data.speciesOrder = state.speciesOrder;
        end

        % Update scatter colors
        for i = 1:numel(data.scatterHandles)
            if isvalid(data.scatterHandles(i))
                set(data.scatterHandles(i), 'MarkerFaceColor', data.colors(i,:), ...
                    'MarkerEdgeColor', data.colors(i,:));
            end
        end

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
        updateScatter(controlFig);
        updateUndoRedoButtons(controlFig);
    catch ME
        warning('scatterPlotPosWidget:redoError', 'Redo error: %s', ME.message);
    end
end

function updateUndoRedoButtons(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(data) || isempty(ctrls)
        return;
    end

    if isfield(ctrls, 'undoBtn') && isgraphics(ctrls.undoBtn)
        if isempty(data.undoStack)
            ctrls.undoBtn.Enable = 'off';
        else
            ctrls.undoBtn.Enable = 'on';
        end
    end

    if isfield(ctrls, 'redoBtn') && isgraphics(ctrls.redoBtn)
        if isempty(data.redoStack)
            ctrls.redoBtn.Enable = 'off';
        else
            ctrls.redoBtn.Enable = 'on';
        end
    end
end

%% Live Update / Deferred Updates
function onLiveUpdateToggle(src, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(data) || isempty(ctrls)
        return;
    end

    data.liveUpdate = logical(src.Value);
    setappdata(controlFig, 'scatterPlotPosWidget', data);

    % Enable/disable apply button
    if isfield(ctrls, 'applyBtn') && isgraphics(ctrls.applyBtn)
        if data.liveUpdate
            ctrls.applyBtn.Enable = 'off';
        else
            ctrls.applyBtn.Enable = 'on';
        end
    end
end

function onApplyChanges(controlFig)
    updateTable(controlFig);
    updateScatter(controlFig);
    updateStatusBar(controlFig);
end

function conditionalUpdate(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end

    updateTable(controlFig);
    if data.liveUpdate
        updateScatter(controlFig);
    end
    updateStatusBar(controlFig);
end

%% Search/Filter Functions
function onSearchChanged(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        searchStr = lower(string(src.String));
        data.searchFilter = searchStr;

        % Apply both search filtering and sorting
        data.filteredIndices = applySorting(data);

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
    catch ME
        warning('scatterPlotPosWidget:searchError', 'Search error: %s', ME.message);
    end
end

function onClearSearch(controlFig)
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if ~isempty(ctrls) && isfield(ctrls, 'searchBox') && isgraphics(ctrls.searchBox)
        ctrls.searchBox.String = '';
        onSearchChanged(ctrls.searchBox, controlFig);
    end
end

%% Sorting Functions
function onSortChanged(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        sortModes = {'none', 'name_asc', 'name_desc', 'count_desc', 'count_asc', 'mass_desc', 'mass_asc'};
        idx = src.Value;
        if idx > numel(sortModes)
            idx = 1;
        end
        data.sortMode = string(sortModes{idx});

        % Apply sorting to filtered indices
        data.filteredIndices = applySorting(data);

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
    catch ME
        warning('scatterPlotPosWidget:sortError', 'Sort error: %s', ME.message);
    end
end

function sortedIndices = applySorting(data)
    % Start with all indices or search-filtered indices
    if data.searchFilter == "" || strlength(data.searchFilter) == 0
        indices = (1:numel(data.speciesNames))';
    else
        matches = contains(lower(data.displayNames), lower(data.searchFilter));
        indices = find(matches);
    end

    if isempty(indices)
        sortedIndices = indices;
        return;
    end

    switch data.sortMode
        case 'name_asc'
            [~, order] = sort(lower(data.displayNames(indices)));
            sortedIndices = indices(order);
        case 'name_desc'
            [~, order] = sort(lower(data.displayNames(indices)), 'descend');
            sortedIndices = indices(order);
        case 'count_desc'
            [~, order] = sort(data.speciesCounts(indices), 'descend');
            sortedIndices = indices(order);
        case 'count_asc'
            [~, order] = sort(data.speciesCounts(indices), 'ascend');
            sortedIndices = indices(order);
        case 'mass_desc'
            masses = arrayfun(@(i) parseAtomicMass(data.displayNames(i)), indices);
            [~, order] = sort(masses, 'descend');
            sortedIndices = indices(order);
        case 'mass_asc'
            masses = arrayfun(@(i) parseAtomicMass(data.displayNames(i)), indices);
            [~, order] = sort(masses, 'ascend');
            sortedIndices = indices(order);
        otherwise
            sortedIndices = indices;
    end
end

function mass = parseAtomicMass(speciesName)
    % Extract atomic mass from species name
    % Handles formats like: "56Fe", "Fe", "56Fe++", "12C", "C2", etc.
    % For molecules like "C2" or "FeO", sums the masses
    %
    % Common atomic masses (approximate)
    atomicMasses = struct( ...
        'H', 1, 'He', 4, 'Li', 7, 'Be', 9, 'B', 11, 'C', 12, 'N', 14, 'O', 16, ...
        'F', 19, 'Ne', 20, 'Na', 23, 'Mg', 24, 'Al', 27, 'Si', 28, 'P', 31, 'S', 32, ...
        'Cl', 35, 'Ar', 40, 'K', 39, 'Ca', 40, 'Sc', 45, 'Ti', 48, 'V', 51, 'Cr', 52, ...
        'Mn', 55, 'Fe', 56, 'Co', 59, 'Ni', 59, 'Cu', 64, 'Zn', 65, 'Ga', 70, 'Ge', 73, ...
        'As', 75, 'Se', 79, 'Br', 80, 'Kr', 84, 'Rb', 85, 'Sr', 88, 'Y', 89, 'Zr', 91, ...
        'Nb', 93, 'Mo', 96, 'Tc', 98, 'Ru', 101, 'Rh', 103, 'Pd', 106, 'Ag', 108, 'Cd', 112, ...
        'In', 115, 'Sn', 119, 'Sb', 122, 'Te', 128, 'I', 127, 'Xe', 131, 'Cs', 133, 'Ba', 137, ...
        'La', 139, 'Ce', 140, 'Pr', 141, 'Nd', 144, 'Pm', 145, 'Sm', 150, 'Eu', 152, 'Gd', 157, ...
        'Tb', 159, 'Dy', 163, 'Ho', 165, 'Er', 167, 'Tm', 169, 'Yb', 173, 'Lu', 175, 'Hf', 178, ...
        'Ta', 181, 'W', 184, 'Re', 186, 'Os', 190, 'Ir', 192, 'Pt', 195, 'Au', 197, 'Hg', 201, ...
        'Tl', 204, 'Pb', 207, 'Bi', 209, 'Po', 209, 'At', 210, 'Rn', 222, 'Fr', 223, 'Ra', 226, ...
        'Ac', 227, 'Th', 232, 'Pa', 231, 'U', 238, 'Np', 237, 'Pu', 244);

    nameStr = char(speciesName);

    % Handle special cases
    if isempty(nameStr) || strcmpi(nameStr, 'unranged') || strcmpi(nameStr, 'background')
        mass = 0;
        return;
    end

    % Remove charge state indicators (+ symbols at the end)
    nameStr = regexprep(nameStr, '\++$', '');

    % Try to extract leading isotope number (e.g., "56" from "56Fe")
    tokens = regexp(nameStr, '^(\d+)([A-Z][a-z]?)', 'tokens');
    if ~isempty(tokens)
        mass = str2double(tokens{1}{1});
        return;
    end

    % No leading number - try to identify element and use standard mass
    % Match element symbols (capital letter optionally followed by lowercase)
    tokens = regexp(nameStr, '([A-Z][a-z]?)(\d*)', 'tokens');

    if isempty(tokens)
        mass = 0;
        return;
    end

    mass = 0;
    for i = 1:numel(tokens)
        element = tokens{i}{1};
        countStr = tokens{i}{2};

        if isempty(countStr)
            count = 1;
        else
            count = str2double(countStr);
            if isnan(count)
                count = 1;
            end
        end

        if isfield(atomicMasses, element)
            mass = mass + atomicMasses.(element) * count;
        end
    end

    % If we couldn't parse anything, return 0
    if mass == 0
        mass = 0;
    end
end

%% Display Options - Axis Labels and Scale Cube
function onAxisLabelsToggle(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || ~isValidGraphics(data.ax)
            return;
        end

        data.showAxisLabels = logical(src.Value);
        setappdata(controlFig, 'scatterPlotPosWidget', data);

        if data.showAxisLabels
            data.ax.XLabel.Visible = 'on';
            data.ax.YLabel.Visible = 'on';
            data.ax.ZLabel.Visible = 'on';
        else
            data.ax.XLabel.Visible = 'off';
            data.ax.YLabel.Visible = 'off';
            data.ax.ZLabel.Visible = 'off';
        end
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:axisLabelsError', 'Axis labels toggle error: %s', ME.message);
    end
end

function onAxisTicksToggle(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || ~isValidGraphics(data.ax)
            return;
        end

        data.showAxisTicks = logical(src.Value);
        setappdata(controlFig, 'scatterPlotPosWidget', data);

        if data.showAxisTicks
            data.ax.XAxis.Visible = 'on';
            data.ax.YAxis.Visible = 'on';
            data.ax.ZAxis.Visible = 'on';
        else
            data.ax.XAxis.Visible = 'off';
            data.ax.YAxis.Visible = 'off';
            data.ax.ZAxis.Visible = 'off';
        end
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:axisTicksError', 'Axis ticks toggle error: %s', ME.message);
    end
end

function onScaleCubeToggle(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        data.showScaleCube = logical(src.Value);

        if data.showScaleCube
            data = createScaleCube(data, controlFig);
        else
            data = destroyScaleCube(data);
        end

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:scaleCubeToggleError', 'Scale cube toggle error: %s', ME.message);
    end
end

function onScaleCubeSizeEdit(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        newSize = str2double(src.String);
        if isnan(newSize) || newSize <= 0
            src.String = num2str(data.scaleCubeSize);
            return;
        end

        data.scaleCubeSize = newSize;
        setappdata(controlFig, 'scatterPlotPosWidget', data);

        % Update scale cube if visible
        if data.showScaleCube && isValidGraphics(data.scaleCubeAx)
            updateScaleCubeGeometry(controlFig);
        end
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:scaleCubeSizeError', 'Scale cube size error: %s', ME.message);
    end
end

function data = createScaleCube(data, controlFig)
    % Create scale cube axes as an overlay in the same figure (bottom right corner)
    if ~isValidGraphics(data.scatterFig)
        return;
    end

    % Delete existing scale cube axes if present
    if isValidGraphics(data.scaleCubeAx)
        delete(data.scaleCubeAx);
    end
    if isfield(data, 'scaleCubeLabel') && isValidGraphics(data.scaleCubeLabel)
        delete(data.scaleCubeLabel);
        data.scaleCubeLabel = [];
    end

    % Delete existing listeners if present
    if isfield(data, 'viewSyncListener') && ~isempty(data.viewSyncListener)
        try delete(data.viewSyncListener); catch, end
    end
    if isfield(data, 'xlimSyncListener') && ~isempty(data.xlimSyncListener)
        try delete(data.xlimSyncListener); catch, end
    end
    if isfield(data, 'ylimSyncListener') && ~isempty(data.ylimSyncListener)
        try delete(data.ylimSyncListener); catch, end
    end
    if isfield(data, 'zlimSyncListener') && ~isempty(data.zlimSyncListener)
        try delete(data.zlimSyncListener); catch, end
    end
    if isfield(data, 'camViewAngleSyncListener') && ~isempty(data.camViewAngleSyncListener)
        try delete(data.camViewAngleSyncListener); catch, end
    end
    if isfield(data, 'projectionSyncListener') && ~isempty(data.projectionSyncListener)
        try delete(data.projectionSyncListener); catch, end
    end
    if isfield(data, 'positionSyncListener') && ~isempty(data.positionSyncListener)
        try delete(data.positionSyncListener); catch, end
    end
    if isfield(data, 'dataAspectSyncListener') && ~isempty(data.dataAspectSyncListener)
        try delete(data.dataAspectSyncListener); catch, end
    end

    % Create new axes in the bottom right corner of the scatter figure
    % Position: [left bottom width height] in normalized units
    data.scaleCubeAx = axes(data.scatterFig, ...
        'Position', [0.70 0.055 0.28 0.245], ...
        'ActivePositionProperty', 'position', ...
        'DataAspectRatio', [1 1 1], ...
        'XTick', [], 'YTick', [], 'ZTick', [], ...
        'Box', 'off', ...
        'Color', 'none', ...           % Transparent background
        'XColor', 'none', ...          % Hide X axis line
        'YColor', 'none', ...          % Hide Y axis line
        'ZColor', 'none', ...          % Hide Z axis line
        'Projection', data.ax.Projection);
    hold(data.scaleCubeAx, 'on');
    axis(data.scaleCubeAx, 'vis3d');
    axis(data.scaleCubeAx, 'off');

    % Make sure scale cube axes is on top
    uistack(data.scaleCubeAx, 'top');

    % Store reference (no separate figure needed)
    data.scaleCubeFig = [];

    % Create cube geometry
    data = drawScaleCube(data);

    % Sync view and limits with main axes
    syncScaleCubeView(data);

    % Setup view and limits synchronization listeners
    if isValidGraphics(data.ax)
        data.viewSyncListener = tryAddAxesPropertyListener(data.ax, 'View', ...
            @(~,~) onMainViewChanged(controlFig));
        % Also sync when axis limits change (zoom/pan)
        data.xlimSyncListener = tryAddAxesPropertyListener(data.ax, 'XLim', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.ylimSyncListener = tryAddAxesPropertyListener(data.ax, 'YLim', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.zlimSyncListener = tryAddAxesPropertyListener(data.ax, 'ZLim', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.camViewAngleSyncListener = tryAddAxesPropertyListener(data.ax, 'CameraViewAngle', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.projectionSyncListener = tryAddAxesPropertyListener(data.ax, 'Projection', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.positionSyncListener = tryAddAxesPropertyListener(data.ax, 'Position', ...
            @(~,~) onMainLimitsChanged(controlFig));
        data.dataAspectSyncListener = tryAddAxesPropertyListener(data.ax, 'DataAspectRatio', ...
            @(~,~) onMainLimitsChanged(controlFig));
    end
end

function data = drawScaleCube(data)
    % Remove old children from scale cube axes
    if isValidGraphics(data.scaleCubeAx)
        cla(data.scaleCubeAx);
    end
    data.scaleCubePatches = [];

    if ~isValidGraphics(data.scaleCubeAx) || ~isValidGraphics(data.ax)
        return;
    end

    data = applyScaleCubeAxisLimits(data);

    s = data.scaleCubeSize / 2;  % Half-size

    % Define cube vertices centered at origin
    vertices = [-s -s -s;
                 s -s -s;
                 s  s -s;
                -s  s -s;
                -s -s  s;
                 s -s  s;
                 s  s  s;
                -s  s  s];

    % Define faces
    % Face colors: X faces = Red, Y faces = Green, Z faces = Blue
    faces = [2 3 7 6;   % +X (red)
             1 5 8 4;   % -X (dark red)
             3 4 8 7;   % +Y (green)
             1 2 6 5;   % -Y (dark green)
             5 6 7 8;   % +Z (blue)
             1 4 3 2];  % -Z (dark blue)

    faceColors = [1.0 0.3 0.3;   % +X red
                  0.6 0.1 0.1;   % -X dark red
                  0.3 1.0 0.3;   % +Y green
                  0.1 0.6 0.1;   % -Y dark green
                  0.3 0.3 1.0;   % +Z blue
                  0.1 0.1 0.6];  % -Z dark blue

    data.scaleCubePatches = gobjects(6, 1);

    hold(data.scaleCubeAx, 'on');
    for i = 1:6
        data.scaleCubePatches(i) = patch(data.scaleCubeAx, ...
            'Vertices', vertices, ...
            'Faces', faces(i,:), ...
            'FaceColor', faceColors(i,:), ...
            'EdgeColor', 'k', ...
            'LineWidth', 1.5, ...
            'FaceAlpha', 1.0);
    end

    data = updateScaleCubeLabelOverlay(data);
end

function updateScaleCubeGeometry(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~isValidGraphics(data.scaleCubeAx)
        return;
    end

    % Clear and redraw
    cla(data.scaleCubeAx);
    hold(data.scaleCubeAx, 'on');
    data = drawScaleCube(data);
    syncScaleCubeView(data);

    % Save updated data (patches may have changed)
    setappdata(controlFig, 'scatterPlotPosWidget', data);
end

function data = destroyScaleCube(data)
    % Delete all sync listeners
    if isfield(data, 'viewSyncListener') && ~isempty(data.viewSyncListener)
        try delete(data.viewSyncListener); catch, end
        data.viewSyncListener = [];
    end
    if isfield(data, 'xlimSyncListener') && ~isempty(data.xlimSyncListener)
        try delete(data.xlimSyncListener); catch, end
        data.xlimSyncListener = [];
    end
    if isfield(data, 'ylimSyncListener') && ~isempty(data.ylimSyncListener)
        try delete(data.ylimSyncListener); catch, end
        data.ylimSyncListener = [];
    end
    if isfield(data, 'zlimSyncListener') && ~isempty(data.zlimSyncListener)
        try delete(data.zlimSyncListener); catch, end
        data.zlimSyncListener = [];
    end
    if isfield(data, 'camViewAngleSyncListener') && ~isempty(data.camViewAngleSyncListener)
        try delete(data.camViewAngleSyncListener); catch, end
        data.camViewAngleSyncListener = [];
    end
    if isfield(data, 'projectionSyncListener') && ~isempty(data.projectionSyncListener)
        try delete(data.projectionSyncListener); catch, end
        data.projectionSyncListener = [];
    end
    if isfield(data, 'positionSyncListener') && ~isempty(data.positionSyncListener)
        try delete(data.positionSyncListener); catch, end
        data.positionSyncListener = [];
    end
    if isfield(data, 'dataAspectSyncListener') && ~isempty(data.dataAspectSyncListener)
        try delete(data.dataAspectSyncListener); catch, end
        data.dataAspectSyncListener = [];
    end

    % Delete the scale cube axes (this also deletes all children including patches)
    if isValidGraphics(data.scaleCubeAx)
        delete(data.scaleCubeAx);
    end
    if isfield(data, 'scaleCubeLabel') && isValidGraphics(data.scaleCubeLabel)
        delete(data.scaleCubeLabel);
    end

    data.scaleCubePatches = [];
    data.scaleCubeFig = [];
    data.scaleCubeAx = [];
    data.scaleCubeLabel = [];
end

function onMainViewChanged(controlFig)
    if ~isValidGraphics(controlFig)
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~data.showScaleCube
        return;
    end
    syncScaleCubeView(data);
end

function onMainLimitsChanged(controlFig)
    % Called when main axes limits change (zoom/pan)
    % Just sync the view - limits will be recalculated
    if ~isValidGraphics(controlFig)
        return;
    end
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~data.showScaleCube
        return;
    end
    syncScaleCubeView(data);
end

function syncScaleCubeView(data)
    if ~isValidGraphics(data.ax)
        return;
    end
    if ~isValidGraphics(data.scaleCubeAx)
        return;
    end

    applyScaleCubeAxisLimits(data);
    applyScaleCubeCamera(data);
    data = updateScaleCubeLabelOverlay(data);
    if isfield(data, 'controlFig') && isValidGraphics(data.controlFig)
        setappdata(data.controlFig, 'scatterPlotPosWidget', data);
    end
end

function data = updateScaleCubeLabelOverlay(data)
    if ~isValidGraphics(data.scatterFig) || ~isValidGraphics(data.scaleCubeAx)
        return;
    end

    sizeStr = sprintf('Scale: %g nm', data.scaleCubeSize);
    labelPos = getScaleCubeLabelPosition(data.scaleCubeAx);

    if ~isfield(data, 'scaleCubeLabel') || ~isValidGraphics(data.scaleCubeLabel)
        data.scaleCubeLabel = uicontrol(data.scatterFig, 'Style', 'text', ...
            'String', sizeStr, ...
            'Units', 'normalized', ...
            'Position', labelPos, ...
            'HorizontalAlignment', 'center', ...
            'BackgroundColor', [1 1 1], ...
            'ForegroundColor', [0 0 0], ...
            'FontWeight', 'bold', ...
            'FontSize', 11, ...
            'Tag', 'scatterPlotPosWidgetScaleCubeLabel', ...
            'Enable', 'inactive');
    else
        data.scaleCubeLabel.String = sizeStr;
        data.scaleCubeLabel.Position = labelPos;
        data.scaleCubeLabel.Visible = 'on';
    end

    try
        uistack(data.scaleCubeLabel, 'top');
    catch
    end
end

function labelPos = getScaleCubeLabelPosition(scaleCubeAx)
    axPos = getAxesInnerPositionSafe(scaleCubeAx);
    labelHeight = max(0.03, min(0.05, axPos(4) * 0.16));
    xPad = axPos(3) * 0.06;
    labelPos = [axPos(1) + xPad, max(0.002, axPos(2) - labelHeight * 0.95), ...
        max(0.10, axPos(3) - 2 * xPad), labelHeight];
end

function data = applyScaleCubeAxisLimits(data)
    if ~isValidGraphics(data.ax) || ~isValidGraphics(data.scaleCubeAx)
        return;
    end

    mainRange = [max(diff(data.ax.XLim), eps), ...
        max(diff(data.ax.YLim), eps), ...
        max(diff(data.ax.ZLim), eps)];

    mainPos = getAxesInnerPositionSafe(data.ax);
    cubePos = getAxesInnerPositionSafe(data.scaleCubeAx);

    ratioW = max(cubePos(3), eps) / max(mainPos(3), eps);
    ratioH = max(cubePos(4), eps) / max(mainPos(4), eps);
    sizeRatio = min(ratioW, ratioH);
    if ~isfinite(sizeRatio) || sizeRatio <= 0
        sizeRatio = 0.25;
    end

    cubeRange = mainRange .* sizeRatio;
    cubeRange = max(cubeRange, eps);
    halfRange = cubeRange / 2;

    xlim(data.scaleCubeAx, [-halfRange(1), halfRange(1)]);
    ylim(data.scaleCubeAx, [-halfRange(2), halfRange(2)]);
    zlim(data.scaleCubeAx, [-halfRange(3), halfRange(3)]);
end

function applyScaleCubeCamera(data)
    if ~isValidGraphics(data.ax) || ~isValidGraphics(data.scaleCubeAx)
        return;
    end

    try
        data.scaleCubeAx.XDir = data.ax.XDir;
        data.scaleCubeAx.YDir = data.ax.YDir;
        data.scaleCubeAx.ZDir = data.ax.ZDir;
        data.scaleCubeAx.Projection = data.ax.Projection;
        data.scaleCubeAx.CameraViewAngle = data.ax.CameraViewAngle;
        data.scaleCubeAx.CameraViewAngleMode = 'manual';
    catch
        % Fallback handled below.
    end

    didCopyCamera = true;
    try
        mainCenter = [mean(data.ax.XLim), mean(data.ax.YLim), mean(data.ax.ZLim)];
        mainRange = [max(diff(data.ax.XLim), eps), ...
            max(diff(data.ax.YLim), eps), ...
            max(diff(data.ax.ZLim), eps)];
        cubeRange = [max(diff(data.scaleCubeAx.XLim), eps), ...
            max(diff(data.scaleCubeAx.YLim), eps), ...
            max(diff(data.scaleCubeAx.ZLim), eps)];

        camPosMain = data.ax.CameraPosition;
        camTargetMain = data.ax.CameraTarget;
        camUpMain = data.ax.CameraUpVector;

        camPosCube = ((camPosMain - mainCenter) ./ mainRange) .* cubeRange;
        camTargetCube = ((camTargetMain - mainCenter) ./ mainRange) .* cubeRange;
        camUpCube = (camUpMain ./ mainRange) .* cubeRange;
        if norm(camUpCube) <= eps
            camUpCube = [0 0 1];
        else
            camUpCube = camUpCube ./ norm(camUpCube);
        end

        data.scaleCubeAx.CameraPositionMode = 'manual';
        data.scaleCubeAx.CameraTargetMode = 'manual';
        data.scaleCubeAx.CameraUpVectorMode = 'manual';
        data.scaleCubeAx.CameraPosition = camPosCube;
        data.scaleCubeAx.CameraTarget = camTargetCube;
        data.scaleCubeAx.CameraUpVector = camUpCube;
    catch
        didCopyCamera = false;
    end

    if ~didCopyCamera
        [az, el] = view(data.ax);
        view(data.scaleCubeAx, az, el);
    end
end

function pos = getAxesInnerPositionSafe(ax)
    try
        pos = ax.InnerPosition;
    catch
        pos = ax.Position;
    end
    if numel(pos) ~= 4
        pos = [0 0 1 1];
    end
end

function listenerHandle = tryAddAxesPropertyListener(ax, propertyName, callbackFn)
    listenerHandle = [];
    if ~isValidGraphics(ax)
        return;
    end
    try
        listenerHandle = addlistener(ax, propertyName, 'PostSet', callbackFn);
    catch
        listenerHandle = [];
    end
end

%% Species Reordering (Z-Order)
function onMoveSpecies(controlFig, direction)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || isempty(data.selectedRow)
            return;
        end

        % Get actual index from filtered view
        if data.selectedRow > numel(data.filteredIndices)
            return;
        end

        actualIdx = data.filteredIndices(data.selectedRow);
        orderIdx = find(data.speciesOrder == actualIdx);

        if isempty(orderIdx)
            return;
        end

        newOrderIdx = orderIdx + direction;
        if newOrderIdx < 1 || newOrderIdx > numel(data.speciesOrder)
            return;
        end

        pushUndo(controlFig);

        % Swap positions in order
        temp = data.speciesOrder(orderIdx);
        data.speciesOrder(orderIdx) = data.speciesOrder(newOrderIdx);
        data.speciesOrder(newOrderIdx) = temp;

        % Reorder scatter handles (affects z-order)
        reorderScatterHandles(data);

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
    catch ME
        warning('scatterPlotPosWidget:moveError', 'Move species error: %s', ME.message);
    end
end

function reorderScatterHandles(data)
    if isempty(data.scatterHandles) || ~isValidGraphics(data.ax)
        return;
    end

    % Reorder children of axes based on speciesOrder
    % Lower index = drawn first = appears behind
    for i = numel(data.speciesOrder):-1:1
        idx = data.speciesOrder(i);
        if idx <= numel(data.scatterHandles) && isvalid(data.scatterHandles(idx))
            uistack(data.scatterHandles(idx), 'top');
        end
    end
end

%% Color Customization
function onChangeColor(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || isempty(data.selectedRow)
            msgbox('Please select a species in the table first.', 'No Selection', 'warn');
            return;
        end

        if data.selectedRow > numel(data.filteredIndices)
            return;
        end

        actualIdx = data.filteredIndices(data.selectedRow);
        currentColor = data.colors(actualIdx, :);

        % Open color picker
        newColor = uisetcolor(currentColor, sprintf('Color for %s', data.displayNames(actualIdx)));

        if isequal(newColor, 0)
            return; % User cancelled
        end

        pushUndo(controlFig);

        data.colors(actualIdx, :) = newColor;

        % Update scatter handle
        if actualIdx <= numel(data.scatterHandles) && isvalid(data.scatterHandles(actualIdx))
            set(data.scatterHandles(actualIdx), 'MarkerFaceColor', newColor, 'MarkerEdgeColor', newColor);
        end

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:colorError', 'Color change error: %s', ME.message);
    end
end

function onTableSelection(~, evd, controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end

    if ~isempty(evd.Indices)
        data.selectedRow = evd.Indices(1);
    else
        data.selectedRow = [];
    end

    setappdata(controlFig, 'scatterPlotPosWidget', data);
end

%% Reset View
function onResetView(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data) || ~isValidGraphics(data.ax)
            return;
        end

        % Reset to default 3D view
        view(data.ax, 3);

        % Reset axis limits to bounding box
        if isfield(data, 'bboxOriginal')
            data.bbox = data.bboxOriginal;
            data.bboxUserEdited = false;
            setappdata(controlFig, 'scatterPlotPosWidget', data);
            applyBoundingBox(data.ax, data.bbox, data.fixBoundingBox);
            updateBoundingBoxControls(controlFig);
        end
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:resetViewError', 'Reset view error: %s', ME.message);
    end
end

function onExportProfileToWorkspace(controlFig)
    try
        profile = visualisationProfileFromWidget(controlFig);
        varName = "visualisationProfile";
        if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
            timestamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            varName = "visualisationProfile_" + timestamp;
        end
        assignin('base', char(varName), profile);
        msgbox(sprintf('Exported visualisation profile to workspace variable ''%s''.', varName), ...
            'Export Complete');
    catch ME
        warning('scatterPlotPosWidget:exportProfileError', ...
            'Export profile error: %s', ME.message);
    end
end

function onImportProfileFromWorkspace(controlFig)
    try
        [~, loaded, sourceName] = scatterPlotPosWidgetLoadState(controlFig);
        if loaded
            msgbox(sprintf('Loaded visualisation profile/state from workspace variable ''%s''.', char(sourceName)), ...
                'Load Complete');
        end
    catch ME
        warning('scatterPlotPosWidget:importProfileError', ...
            'Import profile error: %s', ME.message);
    end
end

%% Status Bar
function updateStatusBar(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(data) || isempty(ctrls) || ~isfield(ctrls, 'statusBar')
        return;
    end

    % Calculate statistics
    totalPoints = height(data.pos);
    visibleCount = 0;

    for i = 1:numel(data.speciesNames)
        if data.visible(i)
            idx = data.speciesIndices{i};
            if data.clipToBoundingBox
                inside = isInsideBBox(data.coords(idx, :), data.bbox);
                idx = idx(inside);
            end
            n = computeSampleCount(numel(idx), data.sampleFractions(i));
            visibleCount = visibleCount + n;
        end
    end

    visibleSpecies = sum(data.visible);
    totalSpecies = numel(data.speciesNames);

    % Estimate memory (rough: 8 bytes per double, 3 coords per point)
    memoryMB = (visibleCount * 3 * 8) / (1024 * 1024);

    statusText = sprintf('Points: %s/%s | Species: %d/%d | ~%.1f MB', ...
        formatNumber(visibleCount), formatNumber(totalPoints), ...
        visibleSpecies, totalSpecies, memoryMB);

    if isgraphics(ctrls.statusBar)
        ctrls.statusBar.String = statusText;
    end
end

function str = formatNumber(n)
    if n >= 1e6
        str = sprintf('%.1fM', n/1e6);
    elseif n >= 1e3
        str = sprintf('%.1fK', n/1e3);
    else
        str = sprintf('%d', n);
    end
end

%% Update Functions
function updateScatter(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        data = ensureMarkerSizeLinkState(data);
        setappdata(controlFig, 'scatterPlotPosWidget', data);

        % Use pre-computed sampling if cache is valid, otherwise compute fresh
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
                    % Check if we can use cached sampling
                    if isfield(data, 'samplingCache') && i <= numel(data.samplingCache) && ...
                            data.samplingCacheValid(i) && numel(data.samplingCache{i}) == n
                        sampleIdx = data.samplingCache{i};
                    else
                        sampleIdx = idx(randperm(numel(idx), n));
                        % Cache the result
                        if isfield(data, 'samplingCache')
                            data.samplingCache{i} = sampleIdx;
                            data.samplingCacheValid(i) = true;
                        end
                    end
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

        setappdata(controlFig, 'scatterPlotPosWidget', data);
        applyBoundingBox(data.ax, data.bbox, data.fixBoundingBox);
        updateStatusBar(controlFig);
        persistStateToAxis(controlFig);
    catch ME
        warning('scatterPlotPosWidget:updateScatterError', 'Update scatter error: %s', ME.message);
    end
end

function n = computeSampleCount(total, sampleValue)
    if sampleValue <= 1
        n = round(total * sampleValue);
    else
        n = round(sampleValue);
    end
    n = max(0, min(total, n));
end

function invalidateSamplingCache(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data.samplingCacheValid(:) = false;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
end

function onSlider(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        value = max(0, min(1, src.Value));
        data.sampleValue = value;
        data.sampleFractions = initSampleFractions(data.speciesCounts, value);
        invalidateSamplingCache(controlFig);
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        data.sampleValue = value;
        data.sampleFractions = initSampleFractions(data.speciesCounts, value);
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
        if ~isempty(ctrls) && isfield(ctrls, 'editBox')
            ctrls.editBox.String = num2str(value, '%.3f');
        end
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:sliderError', 'Slider error: %s', ME.message);
    end
end

function onEdit(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        value = str2double(src.String);
        if isnan(value) || value < 0
            value = data.sampleValue;
        end
        data.sampleValue = value;
        data.sampleFractions = initSampleFractions(data.speciesCounts, value);
        invalidateSamplingCache(controlFig);
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        data.sampleValue = value;
        data.sampleFractions = initSampleFractions(data.speciesCounts, value);
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
        if ~isempty(ctrls) && isfield(ctrls, 'slider')
            ctrls.slider.Value = min(value, 1);
        end
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:editError', 'Edit error: %s', ME.message);
    end
end

function onMarkerSizeSlider(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        value = max(1, src.Value);
        data = applyGlobalMarkerSize(data, value);
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
        if ~isempty(ctrls) && isfield(ctrls, 'markerEdit')
            ctrls.markerEdit.String = num2str(value, '%.2f');
        end
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:markerSliderError', 'Marker slider error: %s', ME.message);
    end
end

function onMarkerSizeEdit(src, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        value = str2double(src.String);
        if isnan(value) || value <= 0
            value = data.markerSizeGlobal;
        end
        data = applyGlobalMarkerSize(data, value);
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
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:markerEditError', 'Marker edit error: %s', ME.message);
    end
end

function onRandomSeedEdit(src, controlFig)
    try
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
        invalidateSamplingCache(controlFig);
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        data.randomSeed = value;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        if data.liveUpdate
            updateScatter(controlFig);
        end
    catch ME
        warning('scatterPlotPosWidget:seedEditError', 'Seed edit error: %s', ME.message);
    end
end

function onRotationSlider(src, controlFig)
    try
        angleDeg = max(0, min(360, src.Value));
        applyRotation(controlFig, angleDeg);
        ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
        if ~isempty(ctrls) && isfield(ctrls, 'rotationEdit') && isgraphics(ctrls.rotationEdit)
            ctrls.rotationEdit.String = num2str(angleDeg, '%.1f');
        end
    catch ME
        warning('scatterPlotPosWidget:rotationSliderError', 'Rotation slider error: %s', ME.message);
    end
end

function onRotationEdit(src, controlFig)
    try
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
    catch ME
        warning('scatterPlotPosWidget:rotationEditError', 'Rotation edit error: %s', ME.message);
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
    invalidateSamplingCache(controlFig);
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    data.rotationAngle = angleDeg;
    if isfield(data, 'coordsOriginal')
        data.coords = rotateCoordsZ(data.coordsOriginal, angleDeg);
    end
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    updateBoundingBoxControls(controlFig);
    conditionalUpdate(controlFig);
end

function onExportImages(controlFig)
    try
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
        baseName = string(baseName);
        ext = lower(string(ext));
        if ~ismember(ext, [".png", ".tif", ".tiff", ".jpg", ".jpeg", ".pdf", ".eps", ".svg"])
            ext = resolveExportExtension(filterIndex);
        end

        visibleIdx = find(data.visible);
        if isempty(visibleIdx)
            return;
        end

        originalVisible = data.visible;
        originalSuspend = data.suspendStatePersistence;
        data.suspendStatePersistence = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        suspendCleanup = onCleanup(@() restoreStatePersistence(controlFig, originalSuspend)); %#ok<NASGU>
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
        data.suspendStatePersistence = originalSuspend;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
        updateScatter(controlFig);
    catch ME
        warning('scatterPlotPosWidget:exportImagesError', 'Export images error: %s', ME.message);
    end
end

function onExportTurntable(controlFig)
    try
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
        baseName = string(baseName);
        ext = string(ext);
        if isempty(ext)
            ext = ".avi";
        end

        visibleIdx = find(data.visible);
        if isempty(visibleIdx)
            return;
        end

        originalVisible = data.visible;
        originalSuspend = data.suspendStatePersistence;
        data.suspendStatePersistence = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        suspendCleanup = onCleanup(@() restoreStatePersistence(controlFig, originalSuspend)); %#ok<NASGU>
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

            exportTurntableFromWidget(controlFig, outFile, 0.5, 30);
        end

        data.visible = originalVisible;
        data.suspendStatePersistence = originalSuspend;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        updateTable(controlFig);
        updateScatter(controlFig);
    catch ME
        warning('scatterPlotPosWidget:exportTurntableError', 'Export turntable error: %s', ME.message);
    end
end

function exportTurntableFromWidget(controlFig, outFile, stepDeg, frameRate)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data) || ~isfield(data, 'ax') || ~isValidGraphics(data.ax)
        return;
    end

    ax = data.ax;
    fig = ancestor(ax, 'figure');
    if ~isValidGraphics(fig)
        return;
    end

    if nargin < 4 || isempty(frameRate)
        frameRate = 30;
    end
    if nargin < 3 || isempty(stepDeg) || stepDeg <= 0
        stepDeg = 0.5;
    end

    nFrames = max(1, round(360 / stepDeg));
    stepDeg = 360 / nFrames;

    state = captureTurntableState(data);
    stateCleanup = onCleanup(@() restoreTurntableState(state));

    lockAxesForTurntable(data);

    [~, ~, ext] = fileparts(outFile);
    if strcmpi(ext, '.mp4')
        writer = VideoWriter(outFile, 'MPEG-4');
    else
        writer = VideoWriter(outFile);
    end
    writer.FrameRate = frameRate;
    writer.Quality = 100;
    open(writer);
    writerCleanup = onCleanup(@() closeVideoWriterSafe(writer));

    for frameIdx = 1:nFrames
        if frameIdx > 1
            rotateAxesAroundZ(ax, stepDeg);
        end

        dataNow = getappdata(controlFig, 'scatterPlotPosWidget');
        if ~isempty(dataNow) && isfield(dataNow, 'showScaleCube') && dataNow.showScaleCube
            syncScaleCubeView(dataNow);
        end
        drawnow;

        writeVideo(writer, getframe(fig));
    end
end

function state = captureTurntableState(data)
    state = struct();
    state.mainAx = [];
    state.main = struct();
    state.hasCube = false;
    state.cubeAx = [];
    state.cube = struct();

    if isfield(data, 'ax') && isValidGraphics(data.ax)
        state.mainAx = data.ax;
        state.main = captureAxesState(data.ax);
    end

    if isfield(data, 'showScaleCube') && data.showScaleCube && ...
            isfield(data, 'scaleCubeAx') && isValidGraphics(data.scaleCubeAx)
        state.hasCube = true;
        state.cubeAx = data.scaleCubeAx;
        state.cube = captureAxesState(data.scaleCubeAx);
    end
end

function stateAx = captureAxesState(ax)
    stateAx = struct();
    props = {'XLim','YLim','ZLim', ...
        'XLimMode','YLimMode','ZLimMode', ...
        'DataAspectRatio','DataAspectRatioMode', ...
        'PlotBoxAspectRatio','PlotBoxAspectRatioMode', ...
        'CameraPosition','CameraPositionMode', ...
        'CameraTarget','CameraTargetMode', ...
        'CameraUpVector','CameraUpVectorMode', ...
        'CameraViewAngle','CameraViewAngleMode', ...
        'Projection'};
    for i = 1:numel(props)
        p = props{i};
        if isprop(ax, p)
            stateAx.(p) = ax.(p);
        end
    end
end

function restoreTurntableState(state)
    if isfield(state, 'mainAx') && isValidGraphics(state.mainAx)
        restoreAxesState(state.mainAx, state.main);
    end
    if isfield(state, 'hasCube') && state.hasCube && ...
            isfield(state, 'cubeAx') && isValidGraphics(state.cubeAx)
        restoreAxesState(state.cubeAx, state.cube);
    end
end

function restoreAxesState(ax, stateAx)
    if ~isValidGraphics(ax) || isempty(stateAx)
        return;
    end

    numericProps = {'XLim','YLim','ZLim','DataAspectRatio','PlotBoxAspectRatio', ...
        'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle','Projection'};
    modeProps = {'XLimMode','YLimMode','ZLimMode', ...
        'DataAspectRatioMode','PlotBoxAspectRatioMode', ...
        'CameraPositionMode','CameraTargetMode', ...
        'CameraUpVectorMode','CameraViewAngleMode'};

    for i = 1:numel(numericProps)
        p = numericProps{i};
        if isfield(stateAx, p) && isprop(ax, p)
            try
                ax.(p) = stateAx.(p);
            catch
            end
        end
    end

    for i = 1:numel(modeProps)
        p = modeProps{i};
        if isfield(stateAx, p) && isprop(ax, p)
            try
                ax.(p) = stateAx.(p);
            catch
            end
        end
    end
end

function lockAxesForTurntable(data)
    if ~isfield(data, 'ax') || ~isValidGraphics(data.ax)
        return;
    end

    ax = data.ax;
    axis(ax, 'vis3d');
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
    ax.DataAspectRatioMode = 'manual';
    ax.PlotBoxAspectRatioMode = 'manual';
    ax.CameraPositionMode = 'manual';
    ax.CameraTargetMode = 'manual';
    ax.CameraUpVectorMode = 'manual';
    ax.CameraViewAngleMode = 'manual';

    if isfield(data, 'showScaleCube') && data.showScaleCube && ...
            isfield(data, 'scaleCubeAx') && isValidGraphics(data.scaleCubeAx)
        axis(data.scaleCubeAx, 'vis3d');
        data.scaleCubeAx.XLimMode = 'manual';
        data.scaleCubeAx.YLimMode = 'manual';
        data.scaleCubeAx.ZLimMode = 'manual';
        data.scaleCubeAx.DataAspectRatioMode = 'manual';
        data.scaleCubeAx.PlotBoxAspectRatioMode = 'manual';
        syncScaleCubeView(data);
    end
end

function rotateAxesAroundZ(ax, stepDeg)
    if ~isValidGraphics(ax)
        return;
    end
    try
        camorbit(ax, stepDeg, 0, 'data', [0 0 1]);
    catch
        axes(ax);
        camorbit(stepDeg, 0);
    end
end

function closeVideoWriterSafe(writer)
    if isempty(writer)
        return;
    end
    try
        close(writer);
    catch
    end
end

function onTableEdit(src, evd, controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end

        row = evd.Indices(1);
        col = evd.Indices(2);

        % Map filtered row to actual index
        if row > numel(data.filteredIndices)
            return;
        end
        actualRow = data.filteredIndices(row);
        data = ensureMarkerSizeLinkState(data);

        pushUndo(controlFig);

        if col == 1  % Show checkbox
            data.visible(actualRow) = logical(evd.NewData);
        elseif col == 5  % Number shown
            newCount = max(0, round(evd.NewData));
            total = data.speciesCounts(actualRow);
            data.sampleFractions(actualRow) = min(1, newCount / max(total, 1));
            invalidateSamplingCache(controlFig);
            data = getappdata(controlFig, 'scatterPlotPosWidget');
            data.sampleFractions(actualRow) = min(1, newCount / max(total, 1));
        elseif col == 6  % Fraction
            newFrac = max(0, min(1, evd.NewData));
            data.sampleFractions(actualRow) = newFrac;
            invalidateSamplingCache(controlFig);
            data = getappdata(controlFig, 'scatterPlotPosWidget');
            data.sampleFractions(actualRow) = newFrac;
        elseif col == 7  % Marker size
            newSize = evd.NewData;
            if ~isnumeric(newSize) || isnan(newSize)
                newSize = data.markerSizes(actualRow);
            end
            newSize = max(0.1, newSize);
            data.markerSizes(actualRow) = newSize;
            data.markerSizeLinked(actualRow) = false;
        elseif col == 8  % Link marker size to global
            isLinked = logical(evd.NewData);
            data.markerSizeLinked(actualRow) = isLinked;
            if isLinked
                data.markerSizes(actualRow) = data.markerSizeGlobal;
            end
        end
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:tableEditError', 'Table edit error: %s', ME.message);
    end
end

function onShowAll(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        data.visible(:) = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:showAllError', 'Show all error: %s', ME.message);
    end
end

function onHideAll(controlFig)
    try
        data = getappdata(controlFig, 'scatterPlotPosWidget');
        if isempty(data)
            return;
        end
        pushUndo(controlFig);
        data.visible(:) = false;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:hideAllError', 'Hide all error: %s', ME.message);
    end
end

function onGroupingChanged(controlFig, ~, modeGroup)
    try
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
        data.markerSizeLinked = true(numel(speciesNames), 1);

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
        data.speciesOrder = (1:numel(speciesNames))';
        data.filteredIndices = (1:numel(speciesNames))';
        data.searchFilter = "";
        data.sortMode = "none";

        % Reset sampling cache
        data.samplingCache = cell(numel(speciesNames), 1);
        data.samplingCacheValid = false(numel(speciesNames), 1);

        % Clear undo/redo stacks on grouping change
        data.undoStack = {};
        data.redoStack = {};

        setappdata(controlFig, 'scatterPlotPosWidget', data);

        % Clear search box and reset sort dropdown
        if isfield(ctrls, 'searchBox') && isgraphics(ctrls.searchBox)
            ctrls.searchBox.String = '';
        end
        if isfield(ctrls, 'sortDropdown') && isgraphics(ctrls.sortDropdown)
            ctrls.sortDropdown.Value = 1;  % Reset to 'None'
        end

        updateTable(controlFig);
        updateScatter(controlFig);
        updateUndoRedoButtons(controlFig);
    catch ME
        errLoc = "";
        if ~isempty(ME.stack)
            errLoc = sprintf(' (%s:%d)', ME.stack(1).name, ME.stack(1).line);
        end
        warning('scatterPlotPosWidget:groupingChangedError', ...
            'Grouping changed error%s: %s', errLoc, ME.message);
    end
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
            'Tag', 'scatterPlotPosWidgetScatter', ...
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
    panelBg = get(parent, 'BackgroundColor');
    uicontrol(parent, 'Style', 'text', ...
        'String', labelText, ...
        'Units', 'normalized', ...
        'Position', [0.02 y 0.05 0.2], ...
        'HorizontalAlignment', 'left', ...
        'BackgroundColor', panelBg);
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
    try
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
        conditionalUpdate(controlFig);
    catch ME
        warning('scatterPlotPosWidget:bboxChangedError', 'Bounding box changed error: %s', ME.message);
    end
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
                chargePart = string(repmat('+', 1, n));
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

function data = ensureMarkerSizeLinkState(data)
    nSpecies = numel(data.speciesNames);
    if ~isfield(data, 'markerSizes') || numel(data.markerSizes) ~= nSpecies
        data.markerSizes = initMarkerSizes(nSpecies, data.markerSizeGlobal);
    end

    if ~isfield(data, 'markerSizeLinked') || numel(data.markerSizeLinked) ~= nSpecies
        data.markerSizeLinked = true(nSpecies, 1);
        if ~isempty(data.markerSizes)
            tol = max(1e-9, 1e-6 * max(1, data.markerSizeGlobal));
            data.markerSizeLinked = abs(data.markerSizes(:) - data.markerSizeGlobal) <= tol;
        end
    else
        data.markerSizeLinked = logical(data.markerSizeLinked(:));
    end
end

function data = applyGlobalMarkerSize(data, markerSizeGlobal)
    if isempty(markerSizeGlobal) || ~isfinite(markerSizeGlobal) || markerSizeGlobal <= 0
        markerSizeGlobal = data.markerSizeGlobal;
    end
    data = ensureMarkerSizeLinkState(data);
    data.markerSizeGlobal = markerSizeGlobal;
    data.markerSize = markerSizeGlobal;
    linkedIdx = data.markerSizeLinked;
    if any(linkedIdx)
        data.markerSizes(linkedIdx) = markerSizeGlobal;
    end
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
    % Build table data for filtered species only
    data = ensureMarkerSizeLinkState(data);
    nFiltered = numel(data.filteredIndices);
    tableData = cell(nFiltered, 8);

    for j = 1:nFiltered
        i = data.filteredIndices(j);
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

        tableData{j, 1} = data.visible(i);

        % Color swatch as HTML (for display purposes, show as colored text)
        color = data.colors(i, :);
        tableData{j, 2} = sprintf('<html><div style="background-color:rgb(%d,%d,%d);width:20px;height:12px;"></div></html>', ...
            round(color(1)*255), round(color(2)*255), round(color(3)*255));

        if isfield(data, 'displayNames') && numel(data.displayNames) >= i
            tableData{j, 3} = char(data.displayNames(i));
        else
            tableData{j, 3} = char(data.speciesNames(i));
        end
        tableData{j, 4} = total;
        tableData{j, 5} = shown;
        tableData{j, 6} = frac;
        if isfield(data, 'markerSizes') && numel(data.markerSizes) >= i
            tableData{j, 7} = data.markerSizes(i);
        else
            tableData{j, 7} = data.markerSizeGlobal;
        end
        if isfield(data, 'markerSizeLinked') && numel(data.markerSizeLinked) >= i
            tableData{j, 8} = logical(data.markerSizeLinked(i));
        else
            tableData{j, 8} = true;
        end
    end
end

function updateTable(controlFig)
    data = getappdata(controlFig, 'scatterPlotPosWidget');
    if isempty(data)
        return;
    end
    data = ensureMarkerSizeLinkState(data);
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    ctrls = getappdata(controlFig, 'scatterPlotPosWidgetControls');
    if isempty(ctrls) || ~isfield(ctrls, 'table')
        return;
    end
    ctrls.table.Data = buildTableData(data);
    updateBoundingBoxControls(controlFig);
    updateStatusBar(controlFig);
    drawnow;
end
