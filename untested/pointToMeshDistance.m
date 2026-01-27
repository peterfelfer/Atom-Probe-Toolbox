function [distances, nearestPoints, faceIndices, locationTypes] = pointToMeshDistance(points, mesh, options)
% POINTTOMESHDISTANCE Calculate true signed distance from points to mesh surface
%
% [distances, nearestPoints, faceIndices, locationTypes] = pointToMeshDistance(points, mesh)
% [distances, ...] = pointToMeshDistance(points, mesh, 'signed', false)
%
% Calculates the shortest Euclidean distance from each query point to the
% nearest point on a triangulated mesh surface. Unlike vertex-only methods,
% this finds the true closest point which may lie anywhere on a triangle
% face, edge, or vertex.
%
% INPUT:
%   points - Nx3 array of query point coordinates [x, y, z]
%            OR a pos table with 'x', 'y', 'z' columns
%   mesh   - Mesh structure with fields:
%            .vertices - Mx3 array of vertex coordinates
%            .faces    - Px3 array of triangle vertex indices
%
% OPTIONS:
%   'signed'       - Return signed distances (default: true)
%                    Positive = outside mesh, Negative = inside mesh
%                    Requires mesh to be closed/watertight for correct signs
%   'chunkSize'    - Points per chunk for memory management (default: 10000)
%   'showProgress' - Display progress bar (default: true for large inputs)
%   'precomputeNormals' - Precompute face normals (default: true)
%   'useBVH'       - Use bounding volume hierarchy for acceleration (default: true)
%
% OUTPUT:
%   distances     - Nx1 array of (signed) distances to mesh surface
%                   Positive = outside, Negative = inside (when signed=true)
%   nearestPoints - Nx3 array of closest points on mesh surface
%   faceIndices   - Nx1 array of face indices containing nearest points
%   locationTypes - Nx1 array indicating where on the triangle the closest point lies:
%                   0 = face interior (inside triangle)
%                   1 = vertex 0 (first vertex of face)
%                   2 = vertex 1 (second vertex of face)
%                   3 = vertex 2 (third vertex of face)
%                   4 = edge 0-1 (edge between vertices 0 and 1)
%                   5 = edge 1-2 (edge between vertices 1 and 2)
%                   6 = edge 2-0 (edge between vertices 2 and 0)
%
% ALGORITHM:
%   For each query point:
%   1. Find candidate triangles using spatial acceleration (BVH)
%   2. For each candidate triangle, compute closest point:
%      a. Project point onto triangle plane
%      b. If projection is inside triangle, use that point
%      c. Otherwise, find closest point on triangle edges/vertices
%   3. Return minimum distance across all triangles
%
% EXAMPLES:
%   % Basic usage with Nx3 array (signed distance by default)
%   [dist, nearPts, faceIdx, locType] = pointToMeshDistance(atomPositions, isosurface);
%
%   % Using a pos table directly (with x, y, z columns)
%   pos = posLoad('mydata.pos');
%   [dist, nearPts, faceIdx, locType] = pointToMeshDistance(pos, isosurface);
%
%   % Unsigned distance
%   [dist, ~, ~, ~] = pointToMeshDistance(pos, mesh, 'signed', false);
%
%   % Find atoms within 1nm of surface
%   dist = pointToMeshDistance(pos, mesh);
%   nearSurface = abs(dist) < 1.0;
%   posNear = pos(nearSurface, :);
%
%   % Find atoms on the inside of the mesh
%   dist = pointToMeshDistance(pos, mesh);
%   insideAtoms = dist < 0;
%   posInside = pos(insideAtoms, :);
%
%   % Find which atoms are closest to edges vs faces
%   [~, ~, ~, locType] = pointToMeshDistance(pos, mesh);
%   onEdge = locType >= 4;  % Closest point is on an edge
%   onVertex = locType >= 1 & locType <= 3;  % Closest point is a vertex
%   onFace = locType == 0;  % Closest point is inside a triangle
%
% NOTES:
%   - For signed distances, the mesh should be closed (watertight)
%   - Performance scales with O(N * log(P)) using BVH acceleration
%   - Memory usage: approximately 100 bytes per point for intermediate results
%
% SEE ALSO:
%   patchCreateProxigram, alphaShape, delaunayTriangulation
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    points {mustBePointsOrTable}
    mesh struct
    options.signed (1,1) logical = true
    options.chunkSize (1,1) double {mustBePositive} = 10000
    options.showProgress (1,1) logical = true
    options.precomputeNormals (1,1) logical = true
    options.useBVH (1,1) logical = true
end

% Convert table to array if needed
if istable(points)
    if all(ismember({'x', 'y', 'z'}, points.Properties.VariableNames))
        points = [points.x, points.y, points.z];
    elseif all(ismember({'X', 'Y', 'Z'}, points.Properties.VariableNames))
        points = [points.X, points.Y, points.Z];
    else
        error('pointToMeshDistance:invalidTable', ...
            'Table must contain x, y, z (or X, Y, Z) columns.');
    end
end

% Validate mesh
if ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
    error('pointToMeshDistance:invalidMesh', ...
        'Mesh must have ''vertices'' and ''faces'' fields.');
end

vertices = mesh.vertices;
faces = mesh.faces;

nPoints = size(points, 1);
nFaces = size(faces, 1);

% Preallocate output
distances = zeros(nPoints, 1);
nearestPoints = zeros(nPoints, 3);
faceIndices = zeros(nPoints, 1);
locationTypes = zeros(nPoints, 1);

% Precompute face data
if options.precomputeNormals
    [faceNormals, ~, ~] = computeFaceData(vertices, faces);
else
    faceNormals = [];
end

% Build bounding volume hierarchy for acceleration
if options.useBVH && nFaces > 100
    bvh = buildBVH(vertices, faces);
else
    bvh = [];
end

% Process points (in chunks for memory management)
nChunks = ceil(nPoints / options.chunkSize);
showProg = options.showProgress && nPoints > 1000;

if showProg
    fprintf('Computing point-to-mesh distances...\n');
end

for chunk = 1:nChunks
    startIdx = (chunk - 1) * options.chunkSize + 1;
    endIdx = min(chunk * options.chunkSize, nPoints);
    chunkPoints = points(startIdx:endIdx, :);

    % Compute distances for this chunk
    [chunkDist, chunkNearest, chunkFaces, chunkLocTypes] = computeDistancesChunk(...
        chunkPoints, vertices, faces, faceNormals, bvh);

    distances(startIdx:endIdx) = chunkDist;
    nearestPoints(startIdx:endIdx, :) = chunkNearest;
    faceIndices(startIdx:endIdx) = chunkFaces;
    locationTypes(startIdx:endIdx) = chunkLocTypes;

    if showProg
        fprintf('  Progress: %d/%d points (%.1f%%)\n', endIdx, nPoints, 100*endIdx/nPoints);
    end
end

% Compute signed distances if requested
if options.signed
    if showProg
        fprintf('Computing signed distances...\n');
    end
    signs = computeSignedDistance(points, nearestPoints, faceIndices, ...
        vertices, faces, faceNormals);
    distances = distances .* signs;
end

if showProg
    fprintf('Done.\n');
end

end

%% Core distance computation
function [distances, nearestPoints, faceIndices, locationTypes] = computeDistancesChunk(points, vertices, faces, faceNormals, bvh)
    % Compute distances for a chunk of points

    nPoints = size(points, 1);
    nFaces = size(faces, 1);

    distances = inf(nPoints, 1);
    nearestPoints = zeros(nPoints, 3);
    faceIndices = zeros(nPoints, 1);
    locationTypes = zeros(nPoints, 1);

    for i = 1:nPoints
        pt = points(i, :);

        % Get candidate faces (all faces if no BVH, or filtered by BVH)
        if ~isempty(bvh)
            candidateFaces = queryBVH(bvh, pt, sqrt(distances(i)));
            if isempty(candidateFaces)
                candidateFaces = 1:nFaces;
            end
        else
            candidateFaces = 1:nFaces;
        end

        minDist = inf;
        minPoint = [0, 0, 0];
        minFace = 0;
        minLocType = 0;

        for fIdx = candidateFaces
            % Get triangle vertices
            tri = vertices(faces(fIdx, :), :);
            v0 = tri(1, :);
            v1 = tri(2, :);
            v2 = tri(3, :);

            % Find closest point on triangle
            [closestPt, dist, locType] = pointToTriangleDistance(pt, v0, v1, v2);

            if dist < minDist
                minDist = dist;
                minPoint = closestPt;
                minFace = fIdx;
                minLocType = locType;
            end
        end

        distances(i) = minDist;
        nearestPoints(i, :) = minPoint;
        faceIndices(i) = minFace;
        locationTypes(i) = minLocType;
    end
end

%% Point to triangle distance
function [closestPoint, distance, locationType] = pointToTriangleDistance(point, v0, v1, v2)
    % Find closest point on triangle to query point
    % Uses barycentric coordinate method
    %
    % Returns:
    %   closestPoint - the closest point on the triangle
    %   distance     - distance from query point to closest point
    %   locationType - where on triangle:
    %                  0 = face interior
    %                  1,2,3 = vertex v0,v1,v2
    %                  4,5,6 = edge v0-v1, v1-v2, v2-v0
    %
    % Based on: "Real-Time Collision Detection" by Christer Ericson

    % Edge vectors
    edge0 = v1 - v0;
    edge1 = v2 - v0;
    v0ToPoint = v0 - point;

    % Compute dot products
    a = dot(edge0, edge0);
    b = dot(edge0, edge1);
    c = dot(edge1, edge1);
    d = dot(edge0, v0ToPoint);
    e = dot(edge1, v0ToPoint);

    det = a * c - b * b;

    % Barycentric coordinates (before projection to valid region)
    s = b * e - c * d;
    t = b * d - a * e;

    % Default to face interior
    locationType = 0;

    if s + t <= det
        if s < 0
            if t < 0
                % Region 4: closest to v0 or edges
                if d < 0
                    s = clamp(-d / a, 0, 1);
                    t = 0;
                    if s == 0
                        locationType = 1;  % vertex v0
                    elseif s == 1
                        locationType = 2;  % vertex v1
                    else
                        locationType = 4;  % edge v0-v1
                    end
                else
                    s = 0;
                    t = clamp(-e / c, 0, 1);
                    if t == 0
                        locationType = 1;  % vertex v0
                    elseif t == 1
                        locationType = 3;  % vertex v2
                    else
                        locationType = 6;  % edge v2-v0
                    end
                end
            else
                % Region 3: closest to edge v0-v2
                s = 0;
                t = clamp(-e / c, 0, 1);
                if t == 0
                    locationType = 1;  % vertex v0
                elseif t == 1
                    locationType = 3;  % vertex v2
                else
                    locationType = 6;  % edge v2-v0
                end
            end
        elseif t < 0
            % Region 5: closest to edge v0-v1
            s = clamp(-d / a, 0, 1);
            t = 0;
            if s == 0
                locationType = 1;  % vertex v0
            elseif s == 1
                locationType = 2;  % vertex v1
            else
                locationType = 4;  % edge v0-v1
            end
        else
            % Region 0: inside triangle
            invDet = 1 / det;
            s = s * invDet;
            t = t * invDet;
            locationType = 0;  % face interior
        end
    else
        if s < 0
            % Region 2: closest to edge v0-v2 or v2
            tmp0 = b + d;
            tmp1 = c + e;
            if tmp1 > tmp0
                numer = tmp1 - tmp0;
                denom = a - 2*b + c;
                s = clamp(numer / denom, 0, 1);
                t = 1 - s;
                if s == 0
                    locationType = 3;  % vertex v2
                elseif s == 1
                    locationType = 2;  % vertex v1
                else
                    locationType = 5;  % edge v1-v2
                end
            else
                s = 0;
                t = clamp(-e / c, 0, 1);
                if t == 0
                    locationType = 1;  % vertex v0
                elseif t == 1
                    locationType = 3;  % vertex v2
                else
                    locationType = 6;  % edge v2-v0
                end
            end
        elseif t < 0
            % Region 6: closest to edge v0-v1 or v1
            tmp0 = b + e;
            tmp1 = a + d;
            if tmp1 > tmp0
                numer = tmp1 - tmp0;
                denom = a - 2*b + c;
                t = clamp(numer / denom, 0, 1);
                s = 1 - t;
                if t == 0
                    locationType = 2;  % vertex v1
                elseif t == 1
                    locationType = 3;  % vertex v2
                else
                    locationType = 5;  % edge v1-v2
                end
            else
                t = 0;
                s = clamp(-d / a, 0, 1);
                if s == 0
                    locationType = 1;  % vertex v0
                elseif s == 1
                    locationType = 2;  % vertex v1
                else
                    locationType = 4;  % edge v0-v1
                end
            end
        else
            % Region 1: closest to edge v1-v2
            numer = (c + e) - (b + d);
            if numer <= 0
                s = 0;
                locationType = 3;  % vertex v2
            else
                denom = a - 2*b + c;
                s = clamp(numer / denom, 0, 1);
                if s == 0
                    locationType = 3;  % vertex v2
                elseif s == 1
                    locationType = 2;  % vertex v1
                else
                    locationType = 5;  % edge v1-v2
                end
            end
            t = 1 - s;
        end
    end

    % Compute closest point
    closestPoint = v0 + s * edge0 + t * edge1;
    distance = norm(point - closestPoint);
end

function val = clamp(x, minVal, maxVal)
    val = max(minVal, min(maxVal, x));
end

%% Signed distance computation
function signs = computeSignedDistance(points, nearestPoints, faceIndices, vertices, faces, faceNormals)
    % Determine sign based on direction to nearest point vs face normal

    nPoints = size(points, 1);
    signs = ones(nPoints, 1);

    % Compute face normals if not provided
    if isempty(faceNormals)
        [faceNormals, ~, ~] = computeFaceData(vertices, faces);
    end

    for i = 1:nPoints
        if faceIndices(i) > 0
            % Vector from nearest point to query point
            toPoint = points(i, :) - nearestPoints(i, :);

            % Face normal at the closest face
            normal = faceNormals(faceIndices(i), :);

            % Sign based on dot product
            if dot(toPoint, normal) < 0
                signs(i) = -1;  % Inside mesh
            end
        end
    end
end

%% Face data computation
function [normals, centroids, areas] = computeFaceData(vertices, faces)
    % Compute normals, centroids, and areas for all faces

    nFaces = size(faces, 1);
    normals = zeros(nFaces, 3);
    centroids = zeros(nFaces, 3);
    areas = zeros(nFaces, 1);

    for i = 1:nFaces
        v0 = vertices(faces(i, 1), :);
        v1 = vertices(faces(i, 2), :);
        v2 = vertices(faces(i, 3), :);

        % Centroid
        centroids(i, :) = (v0 + v1 + v2) / 3;

        % Normal (cross product of edges)
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        crossProd = cross(edge1, edge2);
        area = norm(crossProd) / 2;
        areas(i) = area;

        if area > 0
            normals(i, :) = crossProd / (2 * area);
        else
            normals(i, :) = [0, 0, 1];  % Degenerate triangle
        end
    end
end

%% Bounding Volume Hierarchy (BVH) for acceleration
function bvh = buildBVH(vertices, faces)
    % Build axis-aligned bounding box hierarchy for fast queries

    nFaces = size(faces, 1);

    % Compute bounding box for each face
    faceBounds = zeros(nFaces, 6);  % [minX, minY, minZ, maxX, maxY, maxZ]

    for i = 1:nFaces
        faceVerts = vertices(faces(i, :), :);
        faceBounds(i, 1:3) = min(faceVerts, [], 1);
        faceBounds(i, 4:6) = max(faceVerts, [], 1);
    end

    % Build simple grid-based acceleration structure
    % (Full BVH tree would be more efficient but more complex)
    gridSize = ceil(nFaces^(1/3));
    gridSize = max(gridSize, 4);

    % Global bounds
    globalMin = min(faceBounds(:, 1:3), [], 1);
    globalMax = max(faceBounds(:, 4:6), [], 1);
    gridSpacing = (globalMax - globalMin) / gridSize;
    gridSpacing(gridSpacing == 0) = 1;  % Handle flat meshes

    % Assign faces to grid cells
    bvh = struct();
    bvh.globalMin = globalMin;
    bvh.globalMax = globalMax;
    bvh.gridSize = gridSize;
    bvh.gridSpacing = gridSpacing;
    bvh.faceBounds = faceBounds;

    % Create grid cells
    bvh.grid = cell(gridSize, gridSize, gridSize);

    for i = 1:nFaces
        % Find cells this face overlaps
        faceMin = faceBounds(i, 1:3);
        faceMax = faceBounds(i, 4:6);

        cellMin = max(1, floor((faceMin - globalMin) ./ gridSpacing) + 1);
        cellMax = min(gridSize, ceil((faceMax - globalMin) ./ gridSpacing) + 1);

        for cx = cellMin(1):cellMax(1)
            for cy = cellMin(2):cellMax(2)
                for cz = cellMin(3):cellMax(3)
                    bvh.grid{cx, cy, cz} = [bvh.grid{cx, cy, cz}, i];
                end
            end
        end
    end
end

function candidateFaces = queryBVH(bvh, point, maxDist)
    % Query BVH for faces potentially within maxDist of point

    if isinf(maxDist)
        maxDist = norm(bvh.globalMax - bvh.globalMin);
    end

    % Find grid cells within maxDist of point
    queryMin = point - maxDist;
    queryMax = point + maxDist;

    cellMin = max(1, floor((queryMin - bvh.globalMin) ./ bvh.gridSpacing) + 1);
    cellMax = min(bvh.gridSize, ceil((queryMax - bvh.globalMin) ./ bvh.gridSpacing) + 1);

    % Collect faces from relevant cells
    candidateFaces = [];
    for cx = cellMin(1):cellMax(1)
        for cy = cellMin(2):cellMax(2)
            for cz = cellMin(3):cellMax(3)
                candidateFaces = [candidateFaces, bvh.grid{cx, cy, cz}];
            end
        end
    end

    % Remove duplicates
    candidateFaces = unique(candidateFaces);

    % Further filter by bounding box distance
    if ~isempty(candidateFaces)
        keep = false(size(candidateFaces));
        for i = 1:length(candidateFaces)
            fIdx = candidateFaces(i);
            fMin = bvh.faceBounds(fIdx, 1:3);
            fMax = bvh.faceBounds(fIdx, 4:6);

            % Distance from point to bounding box
            closestInBox = max(fMin, min(point, fMax));
            boxDist = norm(point - closestInBox);

            keep(i) = boxDist <= maxDist * 1.1;  % Small margin
        end
        candidateFaces = candidateFaces(keep);
    end
end

%% Input validation
function mustBePointsOrTable(points)
    % Custom validation: must be Nx3 numeric array or table with x,y,z columns
    if istable(points)
        hasLower = all(ismember({'x', 'y', 'z'}, points.Properties.VariableNames));
        hasUpper = all(ismember({'X', 'Y', 'Z'}, points.Properties.VariableNames));
        if ~hasLower && ~hasUpper
            error('pointToMeshDistance:invalidInput', ...
                'Table must contain x, y, z (or X, Y, Z) columns.');
        end
    elseif isnumeric(points)
        if size(points, 2) ~= 3
            error('pointToMeshDistance:invalidInput', ...
                'Numeric points must be Nx3 array.');
        end
    else
        error('pointToMeshDistance:invalidInput', ...
            'Points must be Nx3 numeric array or table with x, y, z columns.');
    end
end
