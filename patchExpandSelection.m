function [selectedFaces, selectedVertices, info] = patchExpandSelection(mesh, selection, options)
% PATCHEXPANDSELECTION Expand a vertex or face selection on a patch mesh.
%
% [selectedFaces, selectedVertices, info] = patchExpandSelection(mesh, selection)
% [selectedFaces, selectedVertices, info] = patchExpandSelection(mesh, selection, 'select', 'faces')
% [selectedFaces, selectedVertices, info] = patchExpandSelection(mesh, selection, 'maxDistance', 5)
%
% Expands a selection by a given number of adjacency steps. For vertices,
% adjacency is defined by shared edges. For faces, adjacency is defined by
% shared edges (edge-neighbors).
%
% INPUT:
%   mesh      - patch handle or struct with fields:
%               .vertices (Nx3) and .faces (Mx3)
%   selection - indices or logical mask for vertices/faces (see options.select)
%               or struct with fields .verticesMask / .facesMask
%
% OPTIONS:
%   'select'        - 'vertices' (default) or 'faces'
%   'steps'         - number of expansion steps (default: 1). Negative
%                     values shrink the selection by removing boundary rings.
%   'maxDistance'   - geodesic expansion distance along mesh edges
%                     (default: 0, disabled)
%   'includeOriginal' - keep original selection in output (default: true)
%   'returnMask'    - return logical masks (default: false)
%   'axes'          - Target axes for highlighting (default: axes of patchHandle or gca)
%   'patchHandle'   - Patch handle to highlight selection (default: [])
%   'highlight'     - Highlight selection (default: true)
%   'highlightColor'- RGB highlight color (default: [1 0.6 0.2])
%   'highlightHandle' - Existing highlight handle to update/replace
%
% OUTPUT:
%   selectedFaces    - logical mask or indices (see returnMask)
%   selectedVertices - logical mask or indices (see returnMask)
%   info            - struct with fields:
%                     .facesMask, .verticesMask
%                     .facesIndices, .verticesIndices
%
% Note: This function is currently in the "untested" folder.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    mesh
    selection
    options.select (1,1) string {mustBeMember(options.select, ["vertices","faces"])} = "vertices"
    options.steps (1,1) double {mustBeInteger} = 1
    options.maxDistance (1,1) double {mustBeNonnegative} = 0
    options.includeOriginal (1,1) logical = true
    options.returnMask (1,1) logical = false
    options.axes = []
    options.patchHandle = []
    options.highlight = []
    options.highlightColor (1,3) double = [1 0.6 0.2]
    options.highlightHandle = []
end

isPatchInput = isgraphics(mesh, 'patch');
if isPatchInput
    if isempty(options.axes)
        options.axes = ancestor(mesh, 'axes');
    end
    if isempty(options.patchHandle)
        options.patchHandle = mesh;
    end
    mesh = struct('vertices', mesh.Vertices, 'faces', mesh.Faces);
end
if ~isstruct(mesh) || ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
    error('patchExpandSelection:invalidMesh', ...
        'mesh must be a patch handle or a struct with ''vertices'' and ''faces''.');
end

faces = mesh.faces;
vertices = mesh.vertices;
numFaces = size(faces, 1);
numVertices = size(vertices, 1);

if isempty(options.axes)
    options.axes = gca;
end

if isempty(options.highlight)
    options.highlight = ~isempty(options.patchHandle) || isPatchInput;
end

facesMask = false(numFaces, 1);
verticesMask = false(numVertices, 1);

if isstruct(selection)
    if isfield(selection, 'facesMask')
        facesMask = logical(selection.facesMask);
    end
    if isfield(selection, 'verticesMask')
        verticesMask = logical(selection.verticesMask);
    end
else
    if options.select == "faces"
        facesMask = toMask(selection, numFaces);
    else
        verticesMask = toMask(selection, numVertices);
    end
end

useDistance = options.maxDistance > 0;
if options.steps < 0 && useDistance
    error('patchExpandSelection:invalidOptions', ...
        'Negative steps are not supported together with maxDistance.');
end
if options.select == "faces"
    if useDistance
        [facesMask, verticesMask] = expandFacesByDistance(vertices, faces, facesMask, options);
    else
        facesMask = expandFaces(faces, facesMask, options.steps, options.includeOriginal);
        verticesMask = false(numVertices, 1);
        if any(facesMask)
            verticesMask(unique(faces(facesMask, :))) = true;
        end
    end
else
    if useDistance
        [facesMask, verticesMask] = expandVerticesByDistance(vertices, faces, verticesMask, options);
    else
        verticesMask = expandVertices(faces, verticesMask, options.steps, options.includeOriginal, numVertices);
        facesMask = false(numFaces, 1);
        if any(verticesMask)
            usedFaces = any(verticesMask(faces), 2);
            facesMask = usedFaces;
        end
    end
end

info = struct();
info.facesMask = facesMask;
info.verticesMask = verticesMask;
info.facesIndices = find(facesMask);
info.verticesIndices = find(verticesMask);
info.highlightHandle = [];
info.patch = struct('handle', options.patchHandle, 'originalFaceColor', [], ...
    'originalCData', [], 'originalFaceVertexCData', []);

if options.returnMask
    selectedFaces = facesMask;
    selectedVertices = verticesMask;
else
    selectedFaces = info.facesIndices;
    selectedVertices = info.verticesIndices;
end

if options.highlight
    if options.select == "faces"
        if ~isempty(options.patchHandle) && isvalid(options.patchHandle)
            info.patch.originalFaceColor = options.patchHandle.FaceColor;
            info.patch.originalCData = options.patchHandle.CData;
            info.patch.originalFaceVertexCData = options.patchHandle.FaceVertexCData;
            highlightSelection(options.patchHandle, numFaces, facesMask, options.highlightColor);
            storeFaceRestore(options.axes, options.patchHandle, info.patch);
        end
    else
        info.highlightHandle = updateVertexHighlight(options.axes, vertices, verticesMask, ...
            options.highlightColor, options.highlightHandle);
    end
end

end

function mask = toMask(selection, n)
    if isempty(selection)
        mask = false(n, 1);
        return;
    end
    if islogical(selection)
        mask = selection(:);
        if numel(mask) ~= n
            error('patchExpandSelection:invalidSelection', 'Logical mask has wrong length.');
        end
        return;
    end
    mask = false(n, 1);
    idx = selection(:);
    idx = idx(idx >= 1 & idx <= n);
    mask(idx) = true;
end

function facesMask = expandFaces(faces, facesMask, steps, includeOriginal)
    if steps == 0
        return;
    end
    if steps < 0
        facesMask = shrinkFaces(faces, facesMask, abs(steps));
        return;
    end
    adj = faceAdjacency(faces);
    current = facesMask;
    expanded = facesMask;
    for i = 1:steps
        neighbors = adj * current > 0;
        current = neighbors;
        expanded = expanded | neighbors;
    end
    if includeOriginal
        facesMask = expanded;
    else
        facesMask = expanded & ~facesMask;
    end
end

function verticesMask = expandVertices(faces, verticesMask, steps, includeOriginal, numVertices)
    if steps == 0
        return;
    end
    if steps < 0
        verticesMask = shrinkVertices(faces, verticesMask, abs(steps), numVertices);
        return;
    end
    adj = vertexAdjacency(faces, numVertices);
    current = verticesMask;
    expanded = verticesMask;
    for i = 1:steps
        neighbors = adj * current > 0;
        current = neighbors;
        expanded = expanded | neighbors;
    end
    if includeOriginal
        verticesMask = expanded;
    else
        verticesMask = expanded & ~verticesMask;
    end
end

function verticesMask = shrinkVertices(faces, verticesMask, steps, numVertices)
    if steps == 0
        return;
    end
    adj = vertexAdjacency(faces, numVertices);
    totalDegree = full(sum(adj, 2));
    current = verticesMask;
    for i = 1:steps
        if ~any(current)
            break;
        end
        inDegree = full(adj * current);
        boundary = current & (inDegree < totalDegree | totalDegree == 0);
        current = current & ~boundary;
    end
    verticesMask = current;
end

function facesMask = shrinkFaces(faces, facesMask, steps)
    if steps == 0
        return;
    end
    adj = faceAdjacency(faces);
    totalDegree = full(sum(adj, 2));
    current = facesMask;
    for i = 1:steps
        if ~any(current)
            break;
        end
        inDegree = full(adj * current);
        boundary = current & (inDegree < totalDegree | totalDegree == 0);
        current = current & ~boundary;
    end
    facesMask = current;
end

function [facesMask, verticesMask] = expandVerticesByDistance(vertices, faces, seedVerticesMask, options)
    numVertices = size(vertices, 1);
    if ~any(seedVerticesMask)
        facesMask = false(size(faces, 1), 1);
        verticesMask = false(numVertices, 1);
        return;
    end
    dist = vertexGeodesicDistance(vertices, faces, seedVerticesMask);
    verticesMask = dist <= options.maxDistance;
    if ~options.includeOriginal
        verticesMask = verticesMask & ~seedVerticesMask;
    end
    facesMask = any(verticesMask(faces), 2);
end

function [facesMask, verticesMask] = expandFacesByDistance(vertices, faces, seedFacesMask, options)
    numVertices = size(vertices, 1);
    seedVerticesMask = false(numVertices, 1);
    if any(seedFacesMask)
        seedVerticesMask(unique(faces(seedFacesMask, :))) = true;
    end
    if ~any(seedVerticesMask)
        facesMask = false(size(faces, 1), 1);
        verticesMask = false(numVertices, 1);
        return;
    end
    dist = vertexGeodesicDistance(vertices, faces, seedVerticesMask);
    verticesMask = dist <= options.maxDistance;
    if ~options.includeOriginal
        verticesMask = verticesMask & ~seedVerticesMask;
    end
    facesMask = any(verticesMask(faces), 2);
end

function dist = vertexGeodesicDistance(vertices, faces, seedMask)
    numVertices = size(vertices, 1);
    edges = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [1 3])];
    edges = sort(edges, 2);
    edges = unique(edges, 'rows');
    i = edges(:, 1);
    j = edges(:, 2);
    w = sqrt(sum((vertices(i, :) - vertices(j, :)).^2, 2));
    if exist('graph', 'file') ~= 2
        dist = inf(numVertices, 1);
        dist(seedMask) = 0;
        return;
    end
    G = graph(i, j, w, numVertices);
    seedIdx = find(seedMask);
    d = distances(G, seedIdx);
    dist = min(d, [], 1)';
end

function adj = vertexAdjacency(faces, numVertices)
    edges = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [1 3])];
    edges = sort(edges, 2);
    i = edges(:, 1);
    j = edges(:, 2);
    adj = sparse([i; j], [j; i], 1, numVertices, numVertices);
end

function adj = faceAdjacency(faces)
    numFaces = size(faces, 1);
    edges = [faces(:, [1 2]); faces(:, [2 3]); faces(:, [1 3])];
    edges = sort(edges, 2);
    faceIdx = repelem((1:numFaces)', 3, 1);
    [~, ~, ic] = unique(edges, 'rows');
    faceGroups = accumarray(ic, faceIdx, [], @(x) {x});

    adj = sparse(numFaces, numFaces);
    for i = 1:numel(faceGroups)
        f = faceGroups{i};
        if numel(f) < 2
            continue;
        end
        if numel(f) == 2
            adj(f(1), f(2)) = 1;
            adj(f(2), f(1)) = 1;
        else
            for a = 1:numel(f)
                for b = a+1:numel(f)
                    adj(f(a), f(b)) = 1;
                    adj(f(b), f(a)) = 1;
                end
            end
        end
    end
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

function h = updateVertexHighlight(ax, vertices, mask, color, hPrev)
    if ~isempty(hPrev) && isvalid(hPrev)
        delete(hPrev);
    else
        fig = ancestor(ax, 'figure');
        hStored = getappdata(fig, 'patchExpandHighlight');
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
    setappdata(fig, 'patchExpandHighlight', h);
end

function storeFaceRestore(ax, patchHandle, patchInfo)
    fig = ancestor(ax, 'figure');
    data = struct();
    data.handle = patchHandle;
    data.originalFaceColor = patchInfo.originalFaceColor;
    data.originalCData = patchInfo.originalCData;
    data.originalFaceVertexCData = patchInfo.originalFaceVertexCData;
    setappdata(fig, 'patchExpandFaceRestore', data);
end
