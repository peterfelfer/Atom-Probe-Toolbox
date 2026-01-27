function [fv, info] = removeTrappedVertex(fv)
% REMOVETRAPPEDVERTEX Remove degree-3 interior vertices from a triangular mesh.
%
% [fv, info] = removeTrappedVertex(fv)
%
% A "trapped" vertex is defined here as an interior vertex with exactly
% 3 unique neighbors (degree 3) and exactly 3 incident faces. The vertex
% is removed and replaced by a single triangle connecting its neighbors.
%
% Note: This function is currently in the "untested" folder.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

validateattributes(fv, {'struct'}, {'scalar'});
assert(isfield(fv, 'vertices') && isfield(fv, 'faces'), ...
    'removeTrappedVertex:invalidMesh', 'fv must have vertices and faces.');

vertices = fv.vertices;
faces = fv.faces;

bndVerts = boundaryVertices(faces);

numVerts = size(vertices, 1);
vertFaces = accumarray(faces(:), repelem((1:size(faces, 1))', 3, 1), ...
    [numVerts, 1], @(x) {x});

removedVerts = [];
removedFaces = [];
newFaces = [];

for v = 1:numVerts
    if ismember(v, bndVerts)
        continue;
    end
    attachedFaces = vertFaces{v};
    if numel(attachedFaces) ~= 3
        continue;
    end

    neigh = unique(faces(attachedFaces, :));
    neigh(neigh == v) = [];
    if numel(neigh) ~= 3
        continue;
    end

    tri = neigh(:)';

    if hasDuplicateFace(faces, tri)
        continue;
    end

    tri = orientTriangle(tri, faces(attachedFaces, :), vertices);

    removedVerts(end+1, 1) = v; %#ok<AGROW>
    removedFaces = [removedFaces; attachedFaces(:)]; %#ok<AGROW>
    newFaces = [newFaces; tri]; %#ok<AGROW>
end

if isempty(removedVerts)
    info = struct('removedVertices', [], 'removedFaces', [], 'newFaces', []);
    return;
end

keepMask = true(size(faces, 1), 1);
keepMask(unique(removedFaces)) = false;
faces = faces(keepMask, :);
faces = [faces; newFaces];

[vertices, faces, oldToNew] = compactMesh(vertices, faces);

fv.vertices = vertices;
fv.faces = faces;

info = struct();
info.removedVertices = removedVerts;
info.removedFaces = unique(removedFaces);
info.newFaces = newFaces;
info.oldToNew = oldToNew;

end

function bnd = boundaryVertices(faces)
    edges = [faces(:, [1 2]); faces(:, [1 3]); faces(:, [2 3])];
    edges = sort(edges, 2);
    [edgesU, ~, c] = unique(edges, 'rows');
    cnt = accumarray(c, 1);
    bndEdges = edgesU(cnt == 1, :);
    bnd = unique(bndEdges(:));
end

function tf = hasDuplicateFace(faces, tri)
    tri = sort(tri);
    facesSorted = sort(faces, 2);
    tf = any(all(facesSorted == tri, 2));
end

function tri = orientTriangle(tri, oldFaces, vertices)
    a = vertices(tri(1), :);
    b = vertices(tri(2), :);
    c = vertices(tri(3), :);
    newNormal = cross(b - a, c - a);

    oldNormals = zeros(size(oldFaces, 1), 3);
    for i = 1:size(oldFaces, 1)
        f = oldFaces(i, :);
        p1 = vertices(f(1), :);
        p2 = vertices(f(2), :);
        p3 = vertices(f(3), :);
        oldNormals(i, :) = cross(p2 - p1, p3 - p1);
    end
    avgNormal = sum(oldNormals, 1);

    if norm(newNormal) == 0 || norm(avgNormal) == 0
        return;
    end
    if dot(newNormal, avgNormal) < 0
        tri = tri([1 3 2]);
    end
end

function [verticesOut, facesOut, oldToNew] = compactMesh(vertices, faces)
    used = unique(faces(:));
    oldToNew = zeros(size(vertices, 1), 1);
    oldToNew(used) = 1:numel(used);
    verticesOut = vertices(used, :);
    facesOut = oldToNew(faces);
end
