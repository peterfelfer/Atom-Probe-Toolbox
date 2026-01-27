function results = test_pointToMeshDistance()
% TEST_POINTTOMESHDISTANCE Test suite for pointToMeshDistance function
%
% results = test_pointToMeshDistance()
%
% Runs various tests to verify correct behavior of the point-to-mesh
% distance calculation.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

fprintf('Running pointToMeshDistance tests...\n\n');

results = struct();
results.passed = 0;
results.failed = 0;
results.tests = {};

% Test 1: Distance to single triangle
results = runTest(results, @test_singleTriangle, 'Single triangle distance');

% Test 2: Distance to triangle vertices
results = runTest(results, @test_triangleVertices, 'Distance to vertices');

% Test 3: Distance to triangle edges
results = runTest(results, @test_triangleEdges, 'Distance to edges');

% Test 4: Distance inside triangle
results = runTest(results, @test_insideTriangle, 'Point projection inside triangle');

% Test 5: Cube mesh
results = runTest(results, @test_cubeMesh, 'Distance to cube mesh');

% Test 6: Sphere mesh (approximate)
results = runTest(results, @test_sphereMesh, 'Distance to sphere mesh');

% Test 7: Signed distance (now default)
results = runTest(results, @test_signedDistance, 'Signed distance');

% Test 8: Location types output
results = runTest(results, @test_locationTypes, 'Location types output');

% Test 9: Table input (pos variable)
results = runTest(results, @test_tableInput, 'Table input (pos variable)');

% Test 10: Performance with many points
results = runTest(results, @test_performance, 'Performance test');

% Summary
fprintf('\n========================================\n');
fprintf('Test Summary: %d passed, %d failed\n', results.passed, results.failed);
fprintf('========================================\n');

end

function results = runTest(results, testFunc, testName)
    try
        success = testFunc();
        if success
            fprintf('[PASS] %s\n', testName);
            results.passed = results.passed + 1;
        else
            fprintf('[FAIL] %s\n', testName);
            results.failed = results.failed + 1;
        end
        results.tests{end+1} = struct('name', testName, 'passed', success, 'error', '');
    catch ME
        fprintf('[ERROR] %s: %s\n', testName, ME.message);
        results.failed = results.failed + 1;
        results.tests{end+1} = struct('name', testName, 'passed', false, 'error', ME.message);
    end
end

%% Test Functions

function success = test_singleTriangle()
    % Test distance to a simple triangle in the XY plane

    mesh.vertices = [0 0 0; 1 0 0; 0 1 0];
    mesh.faces = [1 2 3];

    % Point directly above triangle centroid
    centroid = mean(mesh.vertices, 1);
    testPoint = centroid + [0 0 1];

    [dist, nearPt, ~, locType] = pointToMeshDistance(testPoint, mesh, 'showProgress', false);

    % Distance should be 1 (height above plane), positive since outside
    success = abs(dist - 1) < 1e-10;

    % Nearest point should be the centroid
    success = success && norm(nearPt - centroid) < 1e-10;

    % Location type should be face interior (0)
    success = success && locType == 0;
end

function success = test_triangleVertices()
    % Test that vertices are correctly identified as nearest points

    mesh.vertices = [0 0 0; 2 0 0; 1 2 0];
    mesh.faces = [1 2 3];

    % Points directly above each vertex
    testPoints = mesh.vertices + [0 0 1; 0 0 1; 0 0 1];

    [dist, nearPts, ~, locTypes] = pointToMeshDistance(testPoints, mesh, 'showProgress', false);

    success = true;

    % All distances should be 1
    success = success && all(abs(dist - 1) < 1e-10);

    % Nearest points should be the vertices
    for i = 1:3
        success = success && norm(nearPts(i,:) - mesh.vertices(i,:)) < 1e-10;
    end

    % Location types should be vertices (1, 2, 3)
    success = success && all(locTypes >= 1 & locTypes <= 3);
end

function success = test_triangleEdges()
    % Test distance to triangle edges

    mesh.vertices = [0 0 0; 2 0 0; 1 2 0];
    mesh.faces = [1 2 3];

    % Point outside triangle, nearest to edge midpoint
    edgeMid = (mesh.vertices(1,:) + mesh.vertices(2,:)) / 2;  % [1, 0, 0]
    testPoint = edgeMid + [0 -1 0];  % Below edge in -Y direction

    [dist, nearPt, ~, locType] = pointToMeshDistance(testPoint, mesh, 'showProgress', false);

    success = abs(dist - 1) < 1e-10;
    success = success && norm(nearPt - edgeMid) < 1e-10;

    % Location type should be an edge (4, 5, or 6)
    success = success && locType >= 4 && locType <= 6;
end

function success = test_insideTriangle()
    % Test point projection inside triangle (not at vertex or edge)

    mesh.vertices = [0 0 0; 3 0 0; 0 3 0];
    mesh.faces = [1 2 3];

    % Point directly above interior point
    interiorPt = [1 1 0];
    testPoint = [1 1 2];

    [dist, nearPt, ~, locType] = pointToMeshDistance(testPoint, mesh, 'showProgress', false);

    success = abs(dist - 2) < 1e-10;
    success = success && norm(nearPt - interiorPt) < 1e-10;

    % Location type should be face interior (0)
    success = success && locType == 0;
end

function success = test_cubeMesh()
    % Test distance to a unit cube mesh

    % Create unit cube centered at origin
    mesh = createCubeMesh(1);

    % Test points (signed distances: inside=negative, outside=positive)
    testPoints = [
        0 0 0;      % Center - should be -0.5 (inside)
        0.5 0 0;    % On +X face - should be 0
        1 0 0;      % Outside +X face - should be +0.5
        0.25 0.25 0.25  % Interior point - should be -0.25 (inside)
    ];

    expectedDist = [-0.5; 0; 0.5; -0.25];  % Signed distances

    [dist, ~, ~, ~] = pointToMeshDistance(testPoints, mesh, 'showProgress', false);

    success = all(abs(dist - expectedDist) < 1e-6);
end

function success = test_sphereMesh()
    % Test distance to an approximate sphere mesh

    % Create icosphere approximation
    mesh = createSphereMesh(1, 2);  % Radius 1, 2 subdivisions

    % Test points inside and outside
    testPoints = [
        2 0 0;    % Outside at radius 2 - should be ~+1
        0.5 0 0;  % Inside at radius 0.5 - should be ~-0.5
    ];

    [dist, ~, ~, ~] = pointToMeshDistance(testPoints, mesh, 'showProgress', false);

    % Allow some tolerance due to mesh approximation
    success = abs(dist(1) - 1) < 0.1;  % Outside, positive
    success = success && abs(dist(2) + 0.5) < 0.1;  % Inside, negative
end

function success = test_signedDistance()
    % Test signed distance (inside vs outside)
    % Note: signed=true is now the default

    mesh = createCubeMesh(2);  % 2x2x2 cube centered at origin

    testPoints = [
        0 0 0;      % Inside - should be negative
        2 0 0;      % Outside - should be positive
    ];

    % Test with default (signed=true)
    [dist, ~, ~, ~] = pointToMeshDistance(testPoints, mesh, 'showProgress', false);
    success = dist(1) < 0 && dist(2) > 0;

    % Test explicitly unsigned
    [distUnsigned, ~, ~, ~] = pointToMeshDistance(testPoints, mesh, 'signed', false, 'showProgress', false);
    success = success && all(distUnsigned > 0);  % Both should be positive when unsigned
end

function success = test_locationTypes()
    % Test that locationTypes correctly identifies vertex, edge, and face hits

    mesh.vertices = [0 0 0; 3 0 0; 0 3 0];
    mesh.faces = [1 2 3];

    % Test points that should hit different parts of the triangle
    testPoints = [
        1 1 1;      % Above interior - should be face (0)
        0 0 1;      % Above vertex 0 - should be vertex 1
        3 0 1;      % Above vertex 1 - should be vertex 2
        0 3 1;      % Above vertex 2 - should be vertex 3
        1.5 0 1;    % Above edge v0-v1 midpoint - should be edge 4
        1.5 1.5 1;  % Above edge v1-v2 midpoint - should be edge 5
        0 1.5 1;    % Above edge v2-v0 midpoint - should be edge 6
    ];

    [~, ~, ~, locTypes] = pointToMeshDistance(testPoints, mesh, 'showProgress', false);

    success = true;

    % Check face interior (locType == 0)
    success = success && locTypes(1) == 0;

    % Check vertices (locType == 1, 2, or 3)
    success = success && locTypes(2) == 1;  % vertex 0
    success = success && locTypes(3) == 2;  % vertex 1
    success = success && locTypes(4) == 3;  % vertex 2

    % Check edges (locType == 4, 5, or 6)
    success = success && locTypes(5) == 4;  % edge v0-v1
    success = success && locTypes(6) == 5;  % edge v1-v2
    success = success && locTypes(7) == 6;  % edge v2-v0

    if ~success
        fprintf('  Location types: expected [0,1,2,3,4,5,6], got [%s]\n', ...
            num2str(locTypes'));
    end
end

function success = test_tableInput()
    % Test that table input (pos variable) works correctly

    mesh.vertices = [0 0 0; 3 0 0; 0 3 0];
    mesh.faces = [1 2 3];

    % Create test points as array
    pointsArray = [1 1 1; 0 0 1; 1.5 0 1];

    % Create same points as table (lowercase x, y, z)
    posTable = table(pointsArray(:,1), pointsArray(:,2), pointsArray(:,3), ...
        'VariableNames', {'x', 'y', 'z'});

    % Compute distances with both
    [distArray, nearPtsArray, ~, locTypesArray] = pointToMeshDistance(pointsArray, mesh, 'showProgress', false);
    [distTable, nearPtsTable, ~, locTypesTable] = pointToMeshDistance(posTable, mesh, 'showProgress', false);

    % Results should be identical
    success = all(abs(distArray - distTable) < 1e-10);
    success = success && all(all(abs(nearPtsArray - nearPtsTable) < 1e-10));
    success = success && all(locTypesArray == locTypesTable);

    % Also test uppercase X, Y, Z
    posTableUpper = table(pointsArray(:,1), pointsArray(:,2), pointsArray(:,3), ...
        'VariableNames', {'X', 'Y', 'Z'});
    [distUpper, ~, ~, ~] = pointToMeshDistance(posTableUpper, mesh, 'showProgress', false);
    success = success && all(abs(distArray - distUpper) < 1e-10);
end

function success = test_performance()
    % Test performance with many points

    mesh = createSphereMesh(1, 3);  % Higher resolution sphere

    % Generate random points around sphere
    nPoints = 1000;
    points = randn(nPoints, 3);
    points = points ./ vecnorm(points, 2, 2) * 1.5;  % On sphere of radius 1.5 (outside)

    tic;
    [dist, ~, ~, ~] = pointToMeshDistance(points, mesh, 'showProgress', false);
    elapsed = toc;

    % Should complete in reasonable time and all distances should be ~+0.5 (outside, positive)
    success = elapsed < 30 && all(abs(dist - 0.5) < 0.2);

    fprintf('  (%.3f seconds for %d points, %d faces)\n', elapsed, nPoints, size(mesh.faces, 1));
end

%% Helper functions to create test meshes

function mesh = createCubeMesh(size)
    % Create a cube mesh centered at origin

    s = size / 2;

    mesh.vertices = [
        -s -s -s;  % 1
         s -s -s;  % 2
         s  s -s;  % 3
        -s  s -s;  % 4
        -s -s  s;  % 5
         s -s  s;  % 6
         s  s  s;  % 7
        -s  s  s;  % 8
    ];

    % Two triangles per face
    mesh.faces = [
        1 2 3; 1 3 4;  % -Z face
        5 7 6; 5 8 7;  % +Z face
        1 5 6; 1 6 2;  % -Y face
        4 3 7; 4 7 8;  % +Y face
        1 4 8; 1 8 5;  % -X face
        2 6 7; 2 7 3;  % +X face
    ];
end

function mesh = createSphereMesh(radius, subdivisions)
    % Create an icosphere mesh

    % Start with icosahedron
    t = (1 + sqrt(5)) / 2;

    vertices = [
        -1  t  0;
         1  t  0;
        -1 -t  0;
         1 -t  0;
         0 -1  t;
         0  1  t;
         0 -1 -t;
         0  1 -t;
         t  0 -1;
         t  0  1;
        -t  0 -1;
        -t  0  1;
    ];

    % Normalize to unit sphere
    vertices = vertices ./ vecnorm(vertices, 2, 2);

    faces = [
        1 12 6;
        1 6 2;
        1 2 8;
        1 8 11;
        1 11 12;
        2 6 10;
        6 12 5;
        12 11 3;
        11 8 7;
        8 2 9;
        4 10 5;
        4 5 3;
        4 3 7;
        4 7 9;
        4 9 10;
        5 10 6;
        3 5 12;
        7 3 11;
        9 7 8;
        10 9 2;
    ];

    % Subdivide
    for sub = 1:subdivisions
        newFaces = [];
        edgeMidpoints = containers.Map('KeyType', 'char', 'ValueType', 'double');

        for i = 1:size(faces, 1)
            v1 = faces(i, 1);
            v2 = faces(i, 2);
            v3 = faces(i, 3);

            % Get or create midpoints
            m12 = getOrCreateMidpoint(v1, v2, vertices, edgeMidpoints);
            m23 = getOrCreateMidpoint(v2, v3, vertices, edgeMidpoints);
            m31 = getOrCreateMidpoint(v3, v1, vertices, edgeMidpoints);

            % Update vertices list
            vertices = [vertices; m12; m23; m31];

            nv = size(vertices, 1);
            i12 = nv - 2;
            i23 = nv - 1;
            i31 = nv;

            % Create 4 new faces
            newFaces = [newFaces;
                v1 i12 i31;
                v2 i23 i12;
                v3 i31 i23;
                i12 i23 i31;
            ];
        end

        faces = newFaces;

        % Remove duplicate vertices (from shared edges)
        [vertices, ~, ic] = unique(vertices, 'rows', 'stable');
        faces = ic(faces);
    end

    % Scale to desired radius
    vertices = vertices * radius;

    mesh.vertices = vertices;
    mesh.faces = faces;
end

function midpoint = getOrCreateMidpoint(v1, v2, vertices, ~)
    % Get midpoint between two vertices, normalized to unit sphere
    midpoint = (vertices(v1, :) + vertices(v2, :)) / 2;
    midpoint = midpoint / norm(midpoint);
end
