function [fvOut, info] = centroidalVoronoiTessellation(fvIn, stepsOrTolerance, angleLimitDeg, varargin)
% CENTROIDALVORONOITESSELLATION Centroidal Voronoi relaxation for an fv mesh.
%
% [fvOut, info] = centroidalVoronoiTessellation(fv)
% [fvOut, info] = centroidalVoronoiTessellation(fv, steps, angleLimitDeg)
% [fvOut, info] = centroidalVoronoiTessellation(fv, steps, angleLimitDeg, ...)
%
% Inputs
%   fv        struct with fields .vertices (Nx3) and .faces (Mx3)
%   steps     number of iterations. If < 1, used as relative tolerance.
%   angleLimitDeg    boundary angle limit in degrees (default 30)
%
% Name-value options
%   'maxIterations'  maximum iterations when steps < 1 (default 50)
%   'tolerance'      fallback tolerance if steps is empty (default 1e-3)
%   'showProgress'   show waitbar (default true)
%   'debug'          verbose debug output (default false)
%   'plot'           plot before/after (default false)
%   'snapToOriginal' snap vertices to original mesh each iteration (default false)
%   'snapMethod'     'nearest-point' or 'nearest-vertex' (default 'nearest-point')
%   'snapEvery'      apply snapping every N iterations (default 1)
%
% Note: This function is currently in the "untested" folder.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 2 || ischar(stepsOrTolerance) || isstring(stepsOrTolerance)
    varargin = [{stepsOrTolerance} {angleLimitDeg} varargin];
    stepsOrTolerance = [];
    angleLimitDeg = [];
end

if nargin < 3 || isempty(angleLimitDeg)
    angleLimitDeg = 30;
end

if isempty(stepsOrTolerance)
    stepsOrTolerance = 3;
end

p = inputParser;
p.addParameter('maxIterations', 50, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('tolerance', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('showProgress', true, @(x) islogical(x) && isscalar(x));
p.addParameter('debug', false, @(x) islogical(x) && isscalar(x));
p.addParameter('plot', false, @(x) islogical(x) && isscalar(x));
p.addParameter('snapToOriginal', false, @(x) islogical(x) && isscalar(x));
p.addParameter('snapMethod', 'nearest-point', @(x) ischar(x) || isstring(x));
p.addParameter('snapEvery', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
options = p.Results;
options.snapMethod = char(lower(string(options.snapMethod)));

validateattributes(fvIn, {'struct'}, {'scalar'});
assert(isfield(fvIn, 'vertices') && isfield(fvIn, 'faces'), ...
    'centroidalVoronoiTessellation:invalidMesh', 'Mesh must have vertices and faces.');

fvOut.vertices = fvIn.vertices;
fvOut.faces = fvIn.faces;

numVertices = size(fvOut.vertices, 1);
edgePairs = getEdgePairs(fvOut.faces);
boundaryVertexIndices = boundaryVertices(fvOut.faces);

useTolerance = false;
if stepsOrTolerance < 1
    useTolerance = true;
    relativeTolerance = stepsOrTolerance;
    if isempty(relativeTolerance) || relativeTolerance <= 0
        relativeTolerance = options.tolerance;
    end
    maxIterations = options.maxIterations;
else
    maxIterations = round(stepsOrTolerance);
    relativeTolerance = NaN;
end

info = struct();
info.iterations = 0;
info.meanShift = [];
info.boundaryVertices = boundaryVertexIndices;
info.options = options;
info.snapIterations = [];

if options.debug
    fprintf('CVT: %d vertices, boundary %d\n', numVertices, numel(boundaryVertexIndices));
end

for iterIndex = 1:maxIterations
    vertexNormalVectors = vertexNormals(fvOut.vertices, fvOut.faces);
    newVertices = fvOut.vertices;

    if options.showProgress
        wb = waitbar(0, sprintf('CVT iteration %d', iterIndex));
    else
        wb = [];
    end

    for vertexIndex = 1:numVertices
        if options.showProgress && (mod(vertexIndex, 500) == 0 || vertexIndex == numVertices)
            waitbar(vertexIndex / numVertices, wb);
        end

        currentVertex = fvOut.vertices(vertexIndex, :);
        neighbors = edgeNeighbors(edgePairs, vertexIndex);

        if ismember(vertexIndex, boundaryVertexIndices)
            neighbors = neighbors(ismember(neighbors, boundaryVertexIndices));
            if numel(neighbors) < 2
                continue;
            end
            if numel(neighbors) > 2
                neighborDistances = vecnorm3(fvOut.vertices(neighbors, :) - currentVertex);
                [~, order] = sort(neighborDistances, 'ascend');
                neighbors = neighbors(order(1:2));
            end

            neighborVertices = fvOut.vertices(neighbors, :);
            edgeVec1 = currentVertex - neighborVertices(1, :);
            edgeVec2 = neighborVertices(2, :) - currentVertex;

            len1 = norm(edgeVec1);
            len2 = norm(edgeVec2);
            if len1 == 0 || len2 == 0
                continue;
            end

            midLength = (len1 + len2) / 2;
            if len1 >= len2
                targetPoint = neighborVertices(1, :) + edgeVec1 / len1 * midLength;
            else
                targetPoint = currentVertex + edgeVec2 / len2 * (midLength - len1);
            end

            angleDeg = acosd(max(-1, min(1, dot(edgeVec1, edgeVec2) / (len1 * len2))));
            if angleDeg < angleLimitDeg
                newVertices(vertexIndex, :) = targetPoint;
            end
        else
            if numel(neighbors) < 3
                continue;
            end

            neighborVertices = fvOut.vertices(neighbors, :);
            vertexNormal = vertexNormalVectors(vertexIndex, :);
            if all(vertexNormal == 0)
                continue;
            end

            basis = null(vertexNormal);
            neighborUv = neighborVertices * basis;

            try
                hullIndices = convhull(neighborUv(:, 1), neighborUv(:, 2));
            catch
                if options.debug
                    warning('CVT:badNeighbourhood', 'Bad neighbourhood at vertex %d', vertexIndex);
                end
                continue;
            end

            hullIndices(end) = [];
            hullUv = neighborUv(hullIndices, :);
            centroid2d = polygonCentroid(hullUv);
            currentUv = currentVertex * basis;
            shift2d = centroid2d - currentUv;
            shift3d = (basis * shift2d')';

            newVertices(vertexIndex, :) = currentVertex + shift3d;
        end
    end

    if options.showProgress && ~isempty(wb)
        close(wb);
    end

    shiftVec = newVertices - fvOut.vertices;
    shiftMagnitude = vecnorm3(shiftVec);
    meanShift = mean(shiftMagnitude);

    fvOut.vertices = newVertices;

    if options.snapToOriginal && mod(iterIndex, options.snapEvery) == 0
        fvOut.vertices = snapToMesh(fvOut.vertices, fvIn, options.snapMethod);
        info.snapIterations(end+1, 1) = iterIndex;
    end

    info.iterations = iterIndex;
    info.meanShift(iterIndex, 1) = meanShift;

    if useTolerance
        if iterIndex == 1
            baselineShift = meanShift;
            if baselineShift == 0
                break;
            end
        else
            if meanShift <= baselineShift * relativeTolerance
                break;
            end
        end
    end
end

if options.plot
    figure('Name', 'CVT before');
    patch(fvIn, 'FaceColor', [0 1 1], 'FaceAlpha', 0.8, 'EdgeColor', [0 0 0]);
    axis equal; rotate3d on;

    figure('Name', sprintf('CVT after %d iterations', info.iterations));
    patch(fvOut, 'FaceColor', [0 1 0], 'FaceAlpha', 0.8, 'EdgeColor', [0 0 0]);
    axis equal; rotate3d on;
end

end

function points = snapToMesh(points, mesh, method)
    switch method
        case 'nearest-vertex'
            points = snapToNearestVertex(points, mesh.vertices);
        case 'nearest-point'
            points = snapToNearestPoint(points, mesh.vertices, mesh.faces);
        otherwise
            error('centroidalVoronoiTessellation:badSnapMethod', ...
                'Unknown snap method: %s', method);
    end
end

function points = snapToNearestVertex(points, vertices)
    for i = 1:size(points, 1)
        diffs = vertices - points(i, :);
        d2 = sum(diffs.^2, 2);
        [~, idx] = min(d2);
        points(i, :) = vertices(idx, :);
    end
end

function points = snapToNearestPoint(points, vertices, faces)
    triA = vertices(faces(:, 1), :);
    triB = vertices(faces(:, 2), :);
    triC = vertices(faces(:, 3), :);

    for i = 1:size(points, 1)
        p = points(i, :);
        bestDist2 = inf;
        bestPoint = p;
        for t = 1:size(triA, 1)
            [q, d2] = closestPointOnTriangle(p, triA(t, :), triB(t, :), triC(t, :));
            if d2 < bestDist2
                bestDist2 = d2;
                bestPoint = q;
            end
        end
        points(i, :) = bestPoint;
    end
end

function [q, dist2] = closestPointOnTriangle(p, a, b, c)
    ab = b - a;
    ac = c - a;
    ap = p - a;

    d1 = dot(ab, ap);
    d2 = dot(ac, ap);
    if d1 <= 0 && d2 <= 0
        q = a;
        dist2 = sum((p - q).^2);
        return;
    end

    bp = p - b;
    d3 = dot(ab, bp);
    d4 = dot(ac, bp);
    if d3 >= 0 && d4 <= d3
        q = b;
        dist2 = sum((p - q).^2);
        return;
    end

    vc = d1 * d4 - d3 * d2;
    if vc <= 0 && d1 >= 0 && d3 <= 0
        v = d1 / (d1 - d3);
        q = a + v * ab;
        dist2 = sum((p - q).^2);
        return;
    end

    cp = p - c;
    d5 = dot(ab, cp);
    d6 = dot(ac, cp);
    if d6 >= 0 && d5 <= d6
        q = c;
        dist2 = sum((p - q).^2);
        return;
    end

    vb = d5 * d2 - d1 * d6;
    if vb <= 0 && d2 >= 0 && d6 <= 0
        w = d2 / (d2 - d6);
        q = a + w * ac;
        dist2 = sum((p - q).^2);
        return;
    end

    va = d3 * d6 - d5 * d4;
    if va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        q = b + w * (c - b);
        dist2 = sum((p - q).^2);
        return;
    end

    denom = 1 / (va + vb + vc);
    v = vb * denom;
    w = vc * denom;
    q = a + ab * v + ac * w;
    dist2 = sum((p - q).^2);
end

function neighbors = edgeNeighbors(edgePairs, vertexIndex)
    mask = any(edgePairs == vertexIndex, 2);
    neighbors = unique(edgePairs(mask, :));
    neighbors = neighbors(neighbors ~= vertexIndex);
end

function edgePairs = getEdgePairs(faces)
    if exist('tri2edgeList', 'file') == 2
        edgePairs = tri2edgeList(faces);
        return;
    end
    edges = [faces(:, [1 2]); faces(:, [1 3]); faces(:, [2 3])];
    edgePairs = unique(sort(edges, 2), 'rows');
end

function bnd = boundaryVertices(faces)
    edges = [faces(:, [1 2]); faces(:, [1 3]); faces(:, [2 3])];
    edges = sort(edges, 2);
    [edgesU, ~, c] = unique(edges, 'rows');
    cnt = accumarray(c, 1);
    bndEdges = edgesU(cnt == 1, :);
    bnd = unique(bndEdges(:));
end

function normals = vertexNormals(vertices, faces)
    normals = zeros(size(vertices));
    v1 = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
    faceNormals = cross(v1, v2, 2);

    for i = 1:size(faces, 1)
        f = faces(i, :);
        normals(f, :) = normals(f, :) + faceNormals(i, :);
    end

    n = sqrt(sum(normals.^2, 2));
    n(n == 0) = 1;
    normals = normals ./ n;
end

function n = vecnorm3(v)
    n = sqrt(sum(v.^2, 2));
end
