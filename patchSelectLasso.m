function [selectedFaces, selectedVertices, info] = patchSelectLasso(varargin)
% PATCHSELECTLASSO Select mesh vertices or faces using a screen-space lasso.
%
% [selectedFaces, selectedVertices, info] = patchSelectLasso(mesh)
% [selectedFaces, selectedVertices, info] = patchSelectLasso(mesh, 'axes', ax)
% [selectedFaces, selectedVertices, info] = patchSelectLasso(patchHandle)
% patchSelectLasso('clear')
% patchSelectLasso('clear', ax)
% patchSelectLasso(patchHandle, 'clear')
%
% Draws a lasso in the current view and selects mesh elements whose
% projections fall inside the lasso polygon.
%
% INPUT:
%   mesh - patch handle or struct with fields:
%          .vertices (Nx3) and .faces (Mx3)
%
% OPTIONS:
%   'axes'             - Target axes (default: axes of patchHandle or gca)
%   'patchHandle'      - Patch handle to highlight selection (default: [])
%   'highlight'        - Highlight selected faces on patchHandle (default: true)
%   'highlightColor'   - RGB highlight color (default: [1 0.6 0.2])
%   'select'           - 'vertices' (default) or 'faces'
%   'useVisibleOnly'   - For faces: only select front-facing faces.
%                        For vertices: apply a screen-space depth test.
%                        (default: false)
%   'depthResolution'  - Depth buffer resolution for vertex selection
%                        (default: 512)
%   'depthTolerance'   - Depth tolerance for visibility (default: 1e-6)
%   'returnMask'       - Return logical masks instead of index lists (default: false)
%   'selectionMode'    - 'auto' | 'replace' | 'add' | 'remove' | 'toggle'
%                        (default: 'auto')
%   'useStoredSelection' - Use stored selection in figure appdata when
%                        combining (default: true)
%   'existingSelection'  - Struct with optional fields:
%                          .facesMask (Mx1) and .verticesMask (Nx1)
%   'highlightHandle'  - Existing highlight handle to update/replace
%   'forceOrthographic'- Force orthographic projection during selection
%                        (default: true)
%   'expandSteps'      - Expand selection by adjacency steps (default: 0)
%   'expandDistance'   - Expand selection by geodesic distance (default: 0)
%   'expandIncludeOriginal' - Keep original selection when expanding (default: true)
%
% OUTPUT:
%   selectedFaces    - Indices of faces inside the lasso
%   selectedVertices - Unique vertex indices of selected faces
%   info             - Struct with diagnostic data and original patch state
%
% Note: This function is currently in the "untested" folder.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin >= 1
    if (ischar(varargin{1}) || (isstring(varargin{1}) && isscalar(varargin{1}))) && ...
            strcmpi(string(varargin{1}), "clear")
        args = varargin(2:end);
        if ~isempty(args) && isgraphics(args{1}, 'axes')
            ax = args{1};
            args = args(2:end);
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", 'axes', ax, args{:});
        else
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", args{:});
        end
        return;
    end

    if nargin >= 2 && (ischar(varargin{2}) || (isstring(varargin{2}) && isscalar(varargin{2}))) && ...
            strcmpi(string(varargin{2}), "clear")
        target = varargin{1};
        args = varargin(3:end);
        if isgraphics(target, 'patch')
            ax = ancestor(target, 'axes');
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", 'axes', ax, args{:});
        elseif isgraphics(target, 'axes')
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", 'axes', target, args{:});
        elseif isgraphics(target, 'figure')
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", args{:});
        else
            [selectedFaces, selectedVertices, info] = patchSelectLassoImpl("clear", args{:});
        end
        return;
    end
end

[selectedFaces, selectedVertices, info] = patchSelectLassoImpl(varargin{:});
end

function [selectedFaces, selectedVertices, info] = patchSelectLassoImpl(mesh, options)
% PATCHSELECTLASSOIMPL Implementation for patchSelectLasso.

arguments
    mesh
    options.axes = []
    options.patchHandle = []
    options.highlight (1,1) logical = true
    options.highlightColor (1,3) double = [1 0.6 0.2]
    options.select (1,1) string {mustBeMember(options.select, ["vertices","faces"])} = "vertices"
    options.useVisibleOnly (1,1) logical = false
    options.depthResolution (1,1) double {mustBePositive} = 512
    options.depthTolerance (1,1) double {mustBeNonnegative} = 1e-6
    options.returnMask (1,1) logical = false
    options.selectionMode (1,1) string {mustBeMember(options.selectionMode, ["auto","replace","add","remove","toggle"])} = "auto"
    options.useStoredSelection (1,1) logical = true
    options.existingSelection = struct()
    options.highlightHandle = []
    options.forceOrthographic (1,1) logical = true
    options.expandSteps (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.expandDistance (1,1) double {mustBeNonnegative} = 0
    options.expandIncludeOriginal (1,1) logical = true
end

if ischar(mesh) || (isstring(mesh) && isscalar(mesh))
    if strcmpi(string(mesh), "clear")
        clearSelection(options.axes);
        selectedFaces = [];
        selectedVertices = [];
        info = struct('cleared', true);
        return;
    end
end

if isgraphics(mesh, 'patch')
    if isempty(options.axes)
        options.axes = ancestor(mesh, 'axes');
    end
    if isempty(options.patchHandle)
        options.patchHandle = mesh;
    end
    mesh = struct('vertices', mesh.Vertices, 'faces', mesh.Faces);
end

if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
    error('patchSelectLasso:invalidMesh', ...
        'mesh must be a patch handle or a struct with ''vertices'' and ''faces''.');
end

if isempty(options.axes)
    options.axes = gca;
end

if exist('drawpolygon', 'file') ~= 2
    error('patchSelectLasso:missingDependency', ...
        'drawpolygon is required for lasso selection.');
end

ax = options.axes;
fig = ancestor(ax, 'figure');

oldProjection = ax.Projection;
if options.forceOrthographic
    ax.Projection = 'orthographic';
end
cleanupProjection = onCleanup(@() set(ax, 'Projection', oldProjection));

overlay = createOverlayAxes(fig, ax);
cleanupOverlay = onCleanup(@() deleteOverlay(overlay));

modeState = disableInteractionModes(fig, ax);
cleanupModes = onCleanup(@() restoreInteractionModes(fig, ax, modeState));

roi = drawpolygon(overlay, 'LineWidth', 1.5, 'Color', options.highlightColor);
lassoXY = roi.Position;
delete(roi);

if isempty(lassoXY)
    selectedFaces = [];
    selectedVertices = [];
    info = struct('lasso', lassoXY);
    return;
end

function modeState = disableInteractionModes(fig, ax)
    modeState = struct();
    if isempty(fig) || ~ishandle(fig)
        return;
    end

    modeState.rotate = rotate3d(fig);
    modeState.zoom = zoom(fig);
    modeState.pan = pan(fig);
    modeState.dcm = datacursormode(fig);
    modeState.brush = brush(fig);
    if exist('isplotedit', 'file') == 2
        modeState.plotEditOn = isplotedit(fig);
    else
        modeState.plotEditOn = false;
    end

    if isprop(modeState.rotate, 'Enable')
        modeState.rotateEnable = modeState.rotate.Enable;
        modeState.rotate.Enable = 'off';
    end
    if isprop(modeState.zoom, 'Enable')
        modeState.zoomEnable = modeState.zoom.Enable;
        modeState.zoom.Enable = 'off';
    end
    if isprop(modeState.pan, 'Enable')
        modeState.panEnable = modeState.pan.Enable;
        modeState.pan.Enable = 'off';
    end
    if isprop(modeState.dcm, 'Enable')
        modeState.dcmEnable = modeState.dcm.Enable;
        modeState.dcm.Enable = 'off';
    end
    if isprop(modeState.brush, 'Enable')
        modeState.brushEnable = modeState.brush.Enable;
        modeState.brush.Enable = 'off';
    end
    if modeState.plotEditOn && exist('plotedit', 'file') == 2
        plotedit(fig, 'off');
    end

    modeState.usedDisableDefaultInteractivity = false;
    if nargin >= 2 && ~isempty(ax)
        if isprop(ax, 'Interactions')
            modeState.interactions = ax.Interactions;
        end
        if exist('disableDefaultInteractivity', 'file') == 2
            try
                disableDefaultInteractivity(ax);
                modeState.usedDisableDefaultInteractivity = true;
            catch
            end
        end
        if isprop(ax, 'Interactions')
            try
                ax.Interactions = [];
            catch
            end
        end
    end
end

function restoreInteractionModes(fig, ax, modeState)
    if isempty(fig) || ~ishandle(fig) || isempty(modeState)
        return;
    end
    if isfield(modeState, 'rotate') && isprop(modeState.rotate, 'Enable')
        modeState.rotate.Enable = modeState.rotateEnable;
    end
    if isfield(modeState, 'zoom') && isprop(modeState.zoom, 'Enable')
        modeState.zoom.Enable = modeState.zoomEnable;
    end
    if isfield(modeState, 'pan') && isprop(modeState.pan, 'Enable')
        modeState.pan.Enable = modeState.panEnable;
    end
    if isfield(modeState, 'dcm') && isprop(modeState.dcm, 'Enable')
        modeState.dcm.Enable = modeState.dcmEnable;
    end
    if isfield(modeState, 'brush') && isprop(modeState.brush, 'Enable')
        modeState.brush.Enable = modeState.brushEnable;
    end
    if isfield(modeState, 'plotEditOn') && modeState.plotEditOn && exist('plotedit', 'file') == 2
        plotedit(fig, 'on');
    end

    if nargin >= 2 && ~isempty(ax)
        if isfield(modeState, 'usedDisableDefaultInteractivity') && modeState.usedDisableDefaultInteractivity
            if exist('enableDefaultInteractivity', 'file') == 2
                try
                    enableDefaultInteractivity(ax);
                catch
                end
            end
        end
        if isfield(modeState, 'interactions') && isprop(ax, 'Interactions')
            try
                ax.Interactions = modeState.interactions;
            catch
            end
        end
    end
end

faces = mesh.faces;
vertices = mesh.vertices;

if options.select == "faces"
    centroids = (vertices(faces(:, 1), :) + vertices(faces(:, 2), :) + vertices(faces(:, 3), :)) / 3;
    [x2d, y2d, projInfo] = projectToAxesNormalized(ax, centroids);
    inMask = inpolygon(x2d, y2d, lassoXY(:, 1), lassoXY(:, 2));

    if options.useVisibleOnly
        faceNormals = cross( ...
            vertices(faces(:, 2), :) - vertices(faces(:, 1), :), ...
            vertices(faces(:, 3), :) - vertices(faces(:, 1), :), 2);
        nLen = sqrt(sum(faceNormals.^2, 2));
        faceNormals = faceNormals ./ max(nLen, eps);
        viewDir = ax.CameraTarget - ax.CameraPosition;
        viewDir = viewDir ./ norm(viewDir);
        facingMask = (faceNormals * viewDir') < 0;
        inMask = inMask & facingMask;
    end

    newFacesMask = inMask;
    newVerticesMask = false(size(vertices, 1), 1);
    newVerticesMask(unique(faces(newFacesMask, :))) = true;
else
    [x2d, y2d, projInfo] = projectToAxesNormalized(ax, vertices);
    inMask = inpolygon(x2d, y2d, lassoXY(:, 1), lassoXY(:, 2));
    if options.useVisibleOnly
        visibleMask = vertexVisibilityMask(vertices, faces, ax.CameraPosition, ...
            x2d, y2d, projInfo.viewDepth, options, inMask);
        inMask = inMask & visibleMask;
    end
    newVerticesMask = inMask;
    newFacesMask = false(size(faces, 1), 1);
end

selectionMode = resolveSelectionMode(fig, options.selectionMode);
[facesMask, verticesMask] = combineSelectionMasks( ...
    newFacesMask, newVerticesMask, faces, vertices, options, selectionMode);

if options.expandSteps > 0 || options.expandDistance > 0
    sel = struct('facesMask', facesMask, 'verticesMask', verticesMask);
    [facesMask, verticesMask, expandInfo] = patchExpandSelection(mesh, sel, ...
        'select', options.select, ...
        'steps', options.expandSteps, ...
        'maxDistance', options.expandDistance, ...
        'includeOriginal', options.expandIncludeOriginal, ...
        'returnMask', true);
else
    expandInfo = struct();
end

selectedFaces = find(facesMask);
selectedVertices = find(verticesMask);

info = struct();
info.lasso = lassoXY;
info.projectedPoints = [x2d, y2d];
info.projection = projInfo;
info.selectedFaces = selectedFaces;
info.selectedVertices = selectedVertices;
info.selectionMode = selectionMode;
info.expand = expandInfo;
info.patch = struct('handle', options.patchHandle, 'originalFaceColor', [], ...
    'originalCData', [], 'originalFaceVertexCData', []);
info.highlightHandle = [];

if options.highlight
    if options.select == "faces"
        if ~isempty(options.patchHandle) && isvalid(options.patchHandle)
            info.patch.originalFaceColor = options.patchHandle.FaceColor;
            info.patch.originalCData = options.patchHandle.CData;
            info.patch.originalFaceVertexCData = options.patchHandle.FaceVertexCData;
            highlightSelection(options.patchHandle, size(faces, 1), facesMask, options.highlightColor);
            storeFaceRestore(fig, options.patchHandle, info.patch);
        end
    else
        info.highlightHandle = updateVertexHighlight(ax, vertices, verticesMask, ...
            options.highlightColor, options.highlightHandle);
    end
end

info.selectedFacesMask = facesMask;
info.selectedVerticesMask = verticesMask;
info.selectedFacesIndices = selectedFaces;
info.selectedVerticesIndices = selectedVertices;

if options.returnMask
    selectedFaces = facesMask;
    selectedVertices = verticesMask;
end

end

function overlay = createOverlayAxes(fig, ax)
    overlay = axes('Parent', fig, ...
        'Units', ax.Units, ...
        'Position', ax.Position, ...
        'Color', 'none', ...
        'XLim', [0 1], ...
        'YLim', [0 1], ...
        'Visible', 'off', ...
        'HitTest', 'on', ...
        'PickableParts', 'all');
    overlay.Tag = 'patchSelectLassoOverlay';
end

function deleteOverlay(overlay)
    if ~isempty(overlay) && isvalid(overlay)
        delete(overlay);
    end
end

function [x2d, y2d, info] = projectToAxesNormalized(ax, points)
    camPos = ax.CameraPosition;
    camTarget = ax.CameraTarget;
    camUp = ax.CameraUpVector;

    forward = camTarget - camPos;
    forward = forward / norm(forward);
    right = cross(forward, camUp);
    right = right / norm(right);
    up = cross(right, forward);

    rel = points - camTarget;
    xView = rel * right';
    yView = rel * up';
    zView = (points - camPos) * forward';

    corners = axisCorners(ax);
    relCorners = corners - camTarget;
    xCorners = relCorners * right';
    yCorners = relCorners * up';
    minX = min(xCorners);
    maxX = max(xCorners);
    minY = min(yCorners);
    maxY = max(yCorners);

    x2d = (xView - minX) / (maxX - minX);
    y2d = (yView - minY) / (maxY - minY);

    info = struct('minX', minX, 'maxX', maxX, 'minY', minY, 'maxY', maxY, ...
        'viewDepth', zView, 'viewRight', right, 'viewUp', up, 'viewForward', forward);
end

function corners = axisCorners(ax)
    xl = ax.XLim;
    yl = ax.YLim;
    zl = ax.ZLim;
    corners = [ ...
        xl(1) yl(1) zl(1); ...
        xl(1) yl(1) zl(2); ...
        xl(1) yl(2) zl(1); ...
        xl(1) yl(2) zl(2); ...
        xl(2) yl(1) zl(1); ...
        xl(2) yl(1) zl(2); ...
        xl(2) yl(2) zl(1); ...
        xl(2) yl(2) zl(2)];
end

function highlightSelection(patchHandle, numFaces, selectedFaces, color)
    baseColor = [0.8 0.8 0.8];
    if isnumeric(patchHandle.FaceColor) && numel(patchHandle.FaceColor) == 3
        baseColor = patchHandle.FaceColor;
    end
    if islogical(selectedFaces)
        mask = selectedFaces;
    else
        mask = false(numFaces, 1);
        mask(selectedFaces) = true;
    end
    faceColors = repmat(baseColor, numFaces, 1);
    faceColors(mask, :) = repmat(color, sum(mask), 1);

    patchHandle.FaceColor = 'flat';
    patchHandle.CData = faceColors;
end

function visibleMask = screenDepthMask(x2d, y2d, depth, resolution, tol)
    valid = x2d >= 0 & x2d <= 1 & y2d >= 0 & y2d <= 1 & depth > 0;
    x2d = x2d(valid);
    y2d = y2d(valid);
    depth = depth(valid);

    ix = min(resolution, max(1, floor(x2d * (resolution - 1)) + 1));
    iy = min(resolution, max(1, floor(y2d * (resolution - 1)) + 1));

    idx = sub2ind([resolution, resolution], iy, ix);
    minDepth = accumarray(idx, depth, [resolution * resolution, 1], @min, inf);

    visible = false(numel(valid), 1);
    depthCell = minDepth(idx);
    visible(valid) = depth <= (depthCell + tol);

    visibleMask = visible;
end

function visibleMask = vertexVisibilityMask(vertices, faces, camPos, x2d, y2d, depth, options, candidateMask)
    visibleMask = false(size(vertices, 1), 1);
    candidates = x2d >= 0 & x2d <= 1 & y2d >= 0 & y2d <= 1 & depth > 0;
    if nargin >= 8 && ~isempty(candidateMask)
        candidates = candidates & candidateMask;
    end
    idx = find(candidates);
    if isempty(idx)
        return;
    end

    if exist('TriangleRayIntersection', 'file') ~= 2
        toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
        triPath = fullfile(toolboxRoot, 'analysis', 'utilities_patch_normals', 'TriangleRayIntersection');
        if isfolder(triPath)
            addpath(triPath);
        end
    end

    if exist('TriangleRayIntersection', 'file') ~= 2
        visibleMask = screenDepthMask(x2d, y2d, depth, options.depthResolution, options.depthTolerance);
        return;
    end

    numFaces = size(faces, 1);
    if numFaces == 0
        visibleMask(candidates) = true;
        return;
    end

    vert1 = vertices(faces(:, 1), :);
    vert2 = vertices(faces(:, 2), :);
    vert3 = vertices(faces(:, 3), :);
    optionsTri = struct('ray', 'segment', 'border', 'inclusive');

    tol = max(options.depthTolerance, 1e-6);
    for k = 1:numel(idx)
        vi = idx(k);
        dir = vertices(vi, :) - camPos;
        if all(dir == 0)
            visibleMask(vi) = true;
            continue;
        end
        orig = repmat(camPos, [numFaces, 1]);
        dirs = repmat(dir, [numFaces, 1]);
        [hit, t] = TriangleRayIntersection(orig, dirs, vert1, vert2, vert3, optionsTri);
        if ~any(hit)
            visibleMask(vi) = true;
            continue;
        end
        tHit = t(hit);
        tMin = min(tHit);
        if tMin >= (1 - tol)
            visibleMask(vi) = true;
        end
    end
end

function selectionMode = resolveSelectionMode(fig, selectionMode)
    if selectionMode ~= "auto"
        return;
    end
    mods = get(fig, 'CurrentModifier');
    if any(strcmpi(mods, 'shift'))
        selectionMode = "add";
    elseif any(strcmpi(mods, 'control')) || any(strcmpi(mods, 'command'))
        selectionMode = "remove";
    else
        selectionMode = "replace";
    end
end

function [facesMask, verticesMask] = combineSelectionMasks( ...
    newFacesMask, newVerticesMask, faces, vertices, options, selectionMode)

    [prevFacesMask, prevVerticesMask] = getExistingMasks(options, size(faces, 1), size(vertices, 1));

    switch selectionMode
        case "replace"
            facesMask = newFacesMask;
            verticesMask = newVerticesMask;
        case "add"
            facesMask = prevFacesMask | newFacesMask;
            verticesMask = prevVerticesMask | newVerticesMask;
        case "remove"
            facesMask = prevFacesMask & ~newFacesMask;
            verticesMask = prevVerticesMask & ~newVerticesMask;
        case "toggle"
            facesMask = xor(prevFacesMask, newFacesMask);
            verticesMask = xor(prevVerticesMask, newVerticesMask);
        otherwise
            facesMask = newFacesMask;
            verticesMask = newVerticesMask;
    end

    if options.useStoredSelection
        fig = ancestor(options.axes, 'figure');
        setappdata(fig, 'patchSelectFacesMask', facesMask);
        setappdata(fig, 'patchSelectVerticesMask', verticesMask);
    end
end

function [facesMask, verticesMask] = getExistingMasks(options, numFaces, numVertices)
    facesMask = false(numFaces, 1);
    verticesMask = false(numVertices, 1);

    if isfield(options.existingSelection, 'facesMask')
        facesMask = logical(options.existingSelection.facesMask);
    end
    if isfield(options.existingSelection, 'verticesMask')
        verticesMask = logical(options.existingSelection.verticesMask);
    end

    if options.useStoredSelection
        fig = ancestor(options.axes, 'figure');
        storedFaces = getappdata(fig, 'patchSelectFacesMask');
        storedVertices = getappdata(fig, 'patchSelectVerticesMask');
        if ~isempty(storedFaces)
            facesMask = logical(storedFaces);
        end
        if ~isempty(storedVertices)
            verticesMask = logical(storedVertices);
        end
    end
end

function h = updateVertexHighlight(ax, vertices, mask, color, hPrev)
    if ~isempty(hPrev) && isvalid(hPrev)
        delete(hPrev);
    else
        fig = ancestor(ax, 'figure');
        hStored = getappdata(fig, 'patchSelectHighlight');
        if ~isempty(hStored) && isvalid(hStored)
            delete(hStored);
        end
    end

    wasHold = ishold(ax);
    hold(ax, 'on');
    if any(mask)
        h = scatter3(ax, vertices(mask, 1), vertices(mask, 2), vertices(mask, 3), ...
            25, color, 'filled', 'MarkerEdgeColor', 'k');
    else
        h = [];
    end
    if ~wasHold
        hold(ax, 'off');
    end

    fig = ancestor(ax, 'figure');
    setappdata(fig, 'patchSelectHighlight', h);
end

function storeFaceRestore(fig, patchHandle, patchInfo)
    data = struct();
    data.handle = patchHandle;
    data.originalFaceColor = patchInfo.originalFaceColor;
    data.originalCData = patchInfo.originalCData;
    data.originalFaceVertexCData = patchInfo.originalFaceVertexCData;
    setappdata(fig, 'patchSelectFaceRestore', data);
end

function clearSelection(target)
    if nargin < 1 || isempty(target)
        fig = gcf;
    elseif isgraphics(target, 'figure')
        fig = target;
    elseif isgraphics(target, 'axes') || isgraphics(target, 'patch')
        fig = ancestor(target, 'figure');
    else
        fig = gcf;
    end

    h = getappdata(fig, 'patchSelectHighlight');
    if ~isempty(h) && isvalid(h)
        delete(h);
    end
    setappdata(fig, 'patchSelectHighlight', []);

    hExpand = getappdata(fig, 'patchExpandHighlight');
    if ~isempty(hExpand) && isvalid(hExpand)
        delete(hExpand);
    end
    setappdata(fig, 'patchExpandHighlight', []);

    restore = getappdata(fig, 'patchSelectFaceRestore');
    if ~isempty(restore) && isfield(restore, 'handle') && isvalid(restore.handle)
        restore.handle.FaceColor = restore.originalFaceColor;
        restore.handle.CData = restore.originalCData;
        restore.handle.FaceVertexCData = restore.originalFaceVertexCData;
    end
    setappdata(fig, 'patchSelectFaceRestore', []);

    restoreExpand = getappdata(fig, 'patchExpandFaceRestore');
    if ~isempty(restoreExpand) && isfield(restoreExpand, 'handle') && isvalid(restoreExpand.handle)
        restoreExpand.handle.FaceColor = restoreExpand.originalFaceColor;
        restoreExpand.handle.CData = restoreExpand.originalCData;
        restoreExpand.handle.FaceVertexCData = restoreExpand.originalFaceVertexCData;
    end
    setappdata(fig, 'patchExpandFaceRestore', []);

    setappdata(fig, 'patchSelectFacesMask', []);
    setappdata(fig, 'patchSelectVerticesMask', []);

    overlays = findobj(fig, 'Type', 'axes', 'Tag', 'patchSelectLassoOverlay');
    if ~isempty(overlays)
        delete(overlays);
    end
end
