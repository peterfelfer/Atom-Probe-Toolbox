function [charData, ax] = boundaryCharacter(fv, ori1, ori2, varargin)
% BOUNDARYCHARACTER Analyze and visualize grain boundary character
%
%   [CHARDATA, AX] = BOUNDARYCHARACTER(FV, ORI1, ORI2) analyzes the character
%   of a grain boundary between two crystals with orientations ORI1 and ORI2.
%   Displays the boundary mesh colored by misorientation angle by default.
%
%   The boundary character is determined by:
%   - Misorientation: The rotation relating the two crystal orientations
%   - Boundary plane: The interface normal expressed in each crystal's frame
%
%   Input:
%       fv   - Boundary mesh structure with fields:
%              fv.vertices - Nx3 array of vertex coordinates
%              fv.faces    - Mx3 array of face indices (triangles)
%
%       ori1 - crystalOrientation object or 3x3 rotation matrix for crystal 1
%       ori2 - crystalOrientation object or 3x3 rotation matrix for crystal 2
%
%   Optional Name-Value Pairs:
%       'colorBy'      - What to color the mesh by:
%                        'misorientation' (default) - misorientation angle
%                        'plane1' - boundary plane IPF color in crystal 1 frame
%                        'plane2' - boundary plane IPF color in crystal 2 frame
%                        'csl'    - CSL sigma value (special boundaries)
%       'property'     - Nx1 per-vertex property to weight (optional)
%       'showIPF'      - Logical, show side-by-side IPF triangles (default: true)
%       'showCSL'      - Logical, identify CSL boundaries (default: true)
%       'cslTolerance' - Angular tolerance for CSL detection in degrees (default: 2)
%       'colormap'     - Colormap for scalar values (default: 'parula')
%       'plot'         - Logical, create visualization (default: true)
%
%   Output:
%       charData - Structure with fields:
%                  .misorientation     - Misorientation struct with angle, axis
%                  .boundaryPlane1     - Nx3 boundary normals in crystal 1 frame
%                  .boundaryPlane2     - Nx3 boundary normals in crystal 2 frame
%                  .cslSigma           - CSL sigma value (NaN if not CSL)
%                  .faceColors         - Mx3 RGB colors based on colorBy mode
%                  .vertexNormals      - Nx3 boundary normals in world frame
%
%       ax - Axes handles (struct with .mesh, .ipf1, .ipf2 if showIPF is true)
%
%   Example:
%       % Create boundary mesh (e.g., a flat plane)
%       [X, Y] = meshgrid(-1:0.1:1, -1:0.1:1);
%       Z = zeros(size(X));
%       fv = surf2patch(X, Y, Z, 'triangles');
%
%       % Define two crystal orientations
%       ori1 = crystalOrientation(eye(3), 'cubic');
%       ori2 = crystalOrientation(rotationMatrix3D([1,1,1], 60), 'cubic');
%
%       % Visualize boundary character
%       [charData, ax] = boundaryCharacter(fv, ori1, ori2);
%
%       % Color by boundary plane in crystal 1
%       [charData, ax] = boundaryCharacter(fv, ori1, ori2, 'colorBy', 'plane1');
%
%   See also: crystalOrientation, ipfMesh, ipfHistogram

    %% Parse inputs
    % Extract rotation matrices
    if isa(ori1, 'crystalOrientation')
        R1 = ori1.getRotationMatrix();
    elseif isnumeric(ori1) && all(size(ori1) == [3, 3])
        R1 = ori1;
    else
        error('ori1 must be a 3x3 rotation matrix or crystalOrientation object');
    end

    if isa(ori2, 'crystalOrientation')
        R2 = ori2.getRotationMatrix();
    elseif isnumeric(ori2) && all(size(ori2) == [3, 3])
        R2 = ori2;
    else
        error('ori2 must be a 3x3 rotation matrix or crystalOrientation object');
    end

    p = inputParser;
    addRequired(p, 'fv', @(x) isstruct(x) && isfield(x,'vertices') && isfield(x,'faces'));
    addRequired(p, 'ori1');
    addRequired(p, 'ori2');
    addParameter(p, 'colorBy', 'misorientation', @(x) ismember(lower(x), {'misorientation', 'plane1', 'plane2', 'csl'}));
    addParameter(p, 'property', [], @(x) isempty(x) || isnumeric(x));
    addParameter(p, 'showIPF', true, @islogical);
    addParameter(p, 'showCSL', true, @islogical);
    addParameter(p, 'cslTolerance', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'colormap', 'parula', @(x) ischar(x) || (isnumeric(x) && size(x,2)==3));
    addParameter(p, 'plot', true, @islogical);
    parse(p, fv, ori1, ori2, varargin{:});

    colorBy = lower(p.Results.colorBy);
    showIPF = p.Results.showIPF;
    showCSL = p.Results.showCSL;
    cslTolerance = p.Results.cslTolerance;
    cmapName = p.Results.colormap;
    doPlot = p.Results.plot;

    vertices = fv.vertices;
    faces = fv.faces;
    nVertices = size(vertices, 1);
    nFaces = size(faces, 1);

    %% Step 1: Calculate misorientation
    misori = calculateMisorientation(R1, R2);

    %% Step 2: Calculate boundary normals in world frame
    vertexNormals = calculateVertexNormals(vertices, faces);

    %% Step 3: Transform boundary normals to both crystal frames
    % R1 transforms crystal1 -> world, so R1' transforms world -> crystal1
    boundaryPlane1 = (R1' * vertexNormals')';
    boundaryPlane2 = (R2' * vertexNormals')';

    %% Step 4: Reduce to fundamental zone for IPF coloring
    reducedPlane1 = reduceToFundamentalZone(boundaryPlane1);
    reducedPlane2 = reduceToFundamentalZone(boundaryPlane2);

    %% Step 5: Detect CSL character
    if showCSL
        [cslSigma, cslName] = detectCSL(misori, cslTolerance);
    else
        cslSigma = NaN;
        cslName = '';
    end

    %% Step 6: Calculate colors based on colorBy mode
    switch colorBy
        case 'misorientation'
            % Uniform color based on misorientation angle
            vertexColors = repmat(misori.angle / 62.8, nVertices, 1);  % Normalize to max cubic misorientation
            vertexColors = repmat([vertexColors, vertexColors, vertexColors], 1, 1);
            scalarValues = repmat(misori.angle, nVertices, 1);
            colorLabel = 'Misorientation Angle (°)';

        case 'plane1'
            % IPF color based on boundary plane in crystal 1 frame
            vertexColors = zeros(nVertices, 3);
            for i = 1:nVertices
                vertexColors(i, :) = calculateIPFColor(reducedPlane1(i, :));
            end
            scalarValues = [];
            colorLabel = 'Boundary Plane (Crystal 1)';

        case 'plane2'
            % IPF color based on boundary plane in crystal 2 frame
            vertexColors = zeros(nVertices, 3);
            for i = 1:nVertices
                vertexColors(i, :) = calculateIPFColor(reducedPlane2(i, :));
            end
            scalarValues = [];
            colorLabel = 'Boundary Plane (Crystal 2)';

        case 'csl'
            % Color by CSL sigma value
            if ~isnan(cslSigma)
                % Use a specific color for the CSL type
                cslColors = getCSLColors();
                if isfield(cslColors, sprintf('sigma%d', cslSigma))
                    col = cslColors.(sprintf('sigma%d', cslSigma));
                else
                    col = [0.5, 0.5, 0.5];  % Gray for unknown CSL
                end
                vertexColors = repmat(col, nVertices, 1);
            else
                vertexColors = repmat([0.7, 0.7, 0.7], nVertices, 1);  % Light gray for non-CSL
            end
            scalarValues = repmat(cslSigma, nVertices, 1);
            colorLabel = 'CSL Sigma Value';
    end

    %% Build output structure
    charData.misorientation = misori;
    charData.boundaryPlane1 = boundaryPlane1;
    charData.boundaryPlane2 = boundaryPlane2;
    charData.reducedPlane1 = reducedPlane1;
    charData.reducedPlane2 = reducedPlane2;
    charData.cslSigma = cslSigma;
    charData.cslName = cslName;
    charData.vertexColors = vertexColors;
    charData.vertexNormals = vertexNormals;
    charData.R1 = R1;
    charData.R2 = R2;

    %% Step 7: Visualization
    if doPlot
        ax = struct();

        if showIPF
            % Create figure with mesh and two IPF triangles
            figure('Color', 'w', 'Position', [100, 100, 1400, 500]);

            % Mesh plot
            ax.mesh = subplot(1, 3, 1);
            plotBoundaryMesh(ax.mesh, fv, vertexColors, scalarValues, cmapName, colorLabel, misori, cslSigma, cslName);

            % IPF triangle for crystal 1
            ax.ipf1 = subplot(1, 3, 2);
            plotBoundaryIPF(ax.ipf1, reducedPlane1, 'Crystal 1 Frame');

            % IPF triangle for crystal 2
            ax.ipf2 = subplot(1, 3, 3);
            plotBoundaryIPF(ax.ipf2, reducedPlane2, 'Crystal 2 Frame');
        else
            figure('Color', 'w');
            ax.mesh = axes;
            plotBoundaryMesh(ax.mesh, fv, vertexColors, scalarValues, cmapName, colorLabel, misori, cslSigma, cslName);
        end
    else
        ax = [];
    end
end

%% Helper Functions

function normals = calculateVertexNormals(vertices, faces)
% Calculate vertex normals by averaging face normals

    nVertices = size(vertices, 1);
    nFaces = size(faces, 1);
    normals = zeros(nVertices, 3);

    for i = 1:nFaces
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);

        p1 = vertices(v1, :);
        p2 = vertices(v2, :);
        p3 = vertices(v3, :);

        edge1 = p2 - p1;
        edge2 = p3 - p1;
        faceNormal = cross(edge1, edge2);

        normals(v1, :) = normals(v1, :) + faceNormal;
        normals(v2, :) = normals(v2, :) + faceNormal;
        normals(v3, :) = normals(v3, :) + faceNormal;
    end

    for i = 1:nVertices
        n = norm(normals(i, :));
        if n > 0
            normals(i, :) = normals(i, :) / n;
        end
    end
end

function misori = calculateMisorientation(R1, R2)
% Calculate misorientation between two orientations
% For cubic symmetry, finds the minimum angle among all equivalent misorientations

    % Misorientation: rotation from crystal 1 to crystal 2
    % deltaR = R2 * R1'  (in world frame)
    % or equivalently in crystal frame: R1' * R2
    deltaR = R1' * R2;

    % Get cubic symmetry operators
    symOps = getCubicSymmetryOperators();

    % Find minimum misorientation angle
    minAngle = 180;
    minAxis = [0, 0, 1];
    minR = deltaR;

    for i = 1:size(symOps, 3)
        for j = 1:size(symOps, 3)
            % Apply symmetry: Si * deltaR * Sj'
            Rtest = symOps(:,:,i) * deltaR * symOps(:,:,j)';

            % Extract angle from rotation matrix
            traceR = trace(Rtest);
            traceR = max(-1, min(3, traceR));  % Clamp for numerical stability
            angle = acosd((traceR - 1) / 2);

            if angle < minAngle
                minAngle = angle;
                minR = Rtest;

                % Extract axis
                if angle > 0.01 && angle < 179.99
                    axis = [Rtest(3,2) - Rtest(2,3);
                            Rtest(1,3) - Rtest(3,1);
                            Rtest(2,1) - Rtest(1,2)];
                    axis = axis / norm(axis);
                    minAxis = axis';
                end
            end
        end
    end

    misori.angle = minAngle;
    misori.axis = minAxis;
    misori.matrix = minR;
end

function symOps = getCubicSymmetryOperators()
% Return the 24 cubic symmetry operators

    symOps = zeros(3, 3, 24);

    % Identity
    symOps(:,:,1) = eye(3);

    % 90° rotations about <100>
    symOps(:,:,2) = [1 0 0; 0 0 -1; 0 1 0];   % 90° about [100]
    symOps(:,:,3) = [1 0 0; 0 -1 0; 0 0 -1];  % 180° about [100]
    symOps(:,:,4) = [1 0 0; 0 0 1; 0 -1 0];   % 270° about [100]

    symOps(:,:,5) = [0 0 1; 0 1 0; -1 0 0];   % 90° about [010]
    symOps(:,:,6) = [-1 0 0; 0 1 0; 0 0 -1];  % 180° about [010]
    symOps(:,:,7) = [0 0 -1; 0 1 0; 1 0 0];   % 270° about [010]

    symOps(:,:,8) = [0 -1 0; 1 0 0; 0 0 1];   % 90° about [001]
    symOps(:,:,9) = [-1 0 0; 0 -1 0; 0 0 1];  % 180° about [001]
    symOps(:,:,10) = [0 1 0; -1 0 0; 0 0 1];  % 270° about [001]

    % 120° rotations about <111>
    symOps(:,:,11) = [0 0 1; 1 0 0; 0 1 0];   % 120° about [111]
    symOps(:,:,12) = [0 1 0; 0 0 1; 1 0 0];   % 240° about [111]

    symOps(:,:,13) = [0 0 -1; -1 0 0; 0 1 0]; % 120° about [-1,1,1]
    symOps(:,:,14) = [0 -1 0; 0 0 1; -1 0 0]; % 240° about [-1,1,1]

    symOps(:,:,15) = [0 0 1; -1 0 0; 0 -1 0]; % 120° about [1,-1,1]
    symOps(:,:,16) = [0 -1 0; 0 0 -1; 1 0 0]; % 240° about [1,-1,1]

    symOps(:,:,17) = [0 0 -1; 1 0 0; 0 -1 0]; % 120° about [1,1,-1]
    symOps(:,:,18) = [0 1 0; 0 0 -1; -1 0 0]; % 240° about [1,1,-1]

    % 180° rotations about <110>
    symOps(:,:,19) = [0 1 0; 1 0 0; 0 0 -1];  % 180° about [110]
    symOps(:,:,20) = [0 -1 0; -1 0 0; 0 0 -1];% 180° about [-110]
    symOps(:,:,21) = [-1 0 0; 0 0 1; 0 1 0];  % 180° about [011]
    symOps(:,:,22) = [-1 0 0; 0 0 -1; 0 -1 0];% 180° about [0-11]
    symOps(:,:,23) = [0 0 1; 0 -1 0; 1 0 0];  % 180° about [101]
    symOps(:,:,24) = [0 0 -1; 0 -1 0; -1 0 0];% 180° about [-101]
end

function reduced = reduceToFundamentalZone(directions)
% Reduce directions to standard triangle

    nDirs = size(directions, 1);
    reduced = zeros(nDirs, 3);

    for i = 1:nDirs
        dir = directions(i, :);
        n = norm(dir);
        if n > 0
            dir = dir / n;
        end
        dirAbs = abs(dir);
        dirSorted = sort(dirAbs, 'ascend');
        reduced(i, :) = dirSorted;
    end
end

function [sigma, name] = detectCSL(misori, tolerance)
% Detect if misorientation corresponds to a CSL boundary
% Returns sigma value and name, or NaN if not CSL

    % Common CSL boundaries for cubic crystals
    % Format: [sigma, angle, axis_h, axis_k, axis_l]
    cslTable = [
        3,   60.00,  1, 1, 1;    % Sigma 3 (twin)
        5,   36.87,  1, 0, 0;    % Sigma 5
        7,   38.21,  1, 1, 1;    % Sigma 7
        9,   38.94,  1, 1, 0;    % Sigma 9
        11,  50.48,  1, 1, 0;    % Sigma 11
        13,  22.62,  1, 0, 0;    % Sigma 13a
        13,  27.80,  1, 1, 1;    % Sigma 13b
        15,  48.19,  2, 1, 0;    % Sigma 15
        17,  28.07,  1, 0, 0;    % Sigma 17a
        17,  61.93,  2, 2, 1;    % Sigma 17b
        19,  26.53,  1, 1, 0;    % Sigma 19a
        19,  46.83,  1, 1, 1;    % Sigma 19b
        21,  21.79,  1, 1, 1;    % Sigma 21a
        21,  44.42,  2, 1, 1;    % Sigma 21b
        23,  40.45,  3, 1, 1;    % Sigma 23
        25,  16.26,  1, 0, 0;    % Sigma 25a
        25,  51.68,  3, 3, 1;    % Sigma 25b
        27,  31.59,  1, 1, 0;    % Sigma 27a
        27,  35.43,  2, 1, 0;    % Sigma 27b
        29,  43.60,  1, 0, 0;    % Sigma 29a
        29,  46.40,  2, 2, 1;    % Sigma 29b
    ];

    sigma = NaN;
    name = '';

    angle = misori.angle;
    axis = misori.axis;

    % Normalize axis
    axis = abs(axis);
    axis = sort(axis, 'descend');

    for i = 1:size(cslTable, 1)
        cslAngle = cslTable(i, 2);
        cslAxis = cslTable(i, 3:5);
        cslAxis = cslAxis / norm(cslAxis);
        cslAxis = sort(abs(cslAxis), 'descend');

        % Check angle match
        if abs(angle - cslAngle) <= tolerance
            % Check axis match (allow some tolerance)
            axisDiff = norm(axis - cslAxis);
            if axisDiff < 0.1
                sigma = cslTable(i, 1);
                name = sprintf('Sigma %d', sigma);
                return;
            end
        end
    end
end

function rgb = calculateIPFColor(direction)
% Calculate IPF color for a direction in fundamental zone

    h = direction(1);
    k = direction(2);
    l = direction(3);

    % Normalize so l = 1
    if l > 0
        h = h / l;
        k = k / l;
    else
        rgb = [0, 0, 0];
        return;
    end

    % Color components
    red = h;
    green = k - h;
    blue = 1 - k;

    % Normalize to max = 1
    maxVal = max([red, green, blue]);
    if maxVal > 0
        rgb = [red, green, blue] / maxVal;
    else
        rgb = [0, 0, 1];
    end

    rgb = max(0, min(1, rgb));
end

function colors = getCSLColors()
% Define colors for common CSL types

    colors.sigma3 = [1, 0, 0];      % Red - twin
    colors.sigma5 = [0, 0.7, 0];    % Green
    colors.sigma7 = [0, 0, 1];      % Blue
    colors.sigma9 = [1, 0.5, 0];    % Orange
    colors.sigma11 = [0.5, 0, 0.5]; % Purple
    colors.sigma13 = [0, 0.7, 0.7]; % Cyan
    colors.sigma15 = [0.7, 0.7, 0]; % Yellow
    colors.sigma17 = [1, 0, 0.5];   % Pink
    colors.sigma19 = [0.5, 0.5, 0]; % Olive
    colors.sigma21 = [0, 0.5, 0.5]; % Teal
    colors.sigma25 = [0.5, 0, 0];   % Dark red
    colors.sigma27 = [0, 0.5, 0];   % Dark green
    colors.sigma29 = [0, 0, 0.5];   % Dark blue
end

function plotBoundaryMesh(ax, fv, vertexColors, scalarValues, cmapName, colorLabel, misori, cslSigma, cslName)
% Plot the boundary mesh with colors

    axes(ax);
    hold(ax, 'on');

    % Check if colors are RGB or scalar
    if size(vertexColors, 2) == 3 && isempty(scalarValues)
        % RGB colors (IPF)
        patch(ax, 'Faces', fv.faces, 'Vertices', fv.vertices, ...
              'FaceVertexCData', vertexColors, ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');
    else
        % Scalar values with colormap
        if ~isempty(scalarValues)
            patch(ax, 'Faces', fv.faces, 'Vertices', fv.vertices, ...
                  'FaceVertexCData', scalarValues, ...
                  'FaceColor', 'interp', ...
                  'EdgeColor', 'none');
            colormap(ax, cmapName);
            cb = colorbar(ax);
            cb.Label.String = colorLabel;
        else
            patch(ax, 'Faces', fv.faces, 'Vertices', fv.vertices, ...
                  'FaceVertexCData', vertexColors, ...
                  'FaceColor', 'interp', ...
                  'EdgeColor', 'none');
        end
    end

    axis(ax, 'equal');
    xlabel(ax, 'X'); ylabel(ax, 'Y'); zlabel(ax, 'Z');

    % Title with misorientation info
    titleStr = sprintf('Boundary: %.1f° about [%.2f, %.2f, %.2f]', ...
        misori.angle, misori.axis(1), misori.axis(2), misori.axis(3));
    if ~isnan(cslSigma)
        titleStr = sprintf('%s\n%s', titleStr, cslName);
    end
    title(ax, titleStr, 'FontSize', 11);

    view(ax, 3);
    rotate3d(ax, 'on');
    hold(ax, 'off');
end

function plotBoundaryIPF(ax, reducedNormals, titleStr)
% Plot boundary plane distribution in IPF triangle

    axes(ax);
    hold(ax, 'on');

    % Draw the colored standard triangle background
    drawStandardTriangle(ax);

    % Plot the boundary plane orientations as points
    [x, y] = direction2stereo(reducedNormals);
    scatter(ax, x, y, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);

    % Draw triangle boundary
    drawTriangleBoundary(ax);

    % Labels
    [x001, y001] = direction2stereoSingle(0, 0, 1);
    [x011, y011] = direction2stereoSingle(0, 1, 1);
    [x111, y111] = direction2stereoSingle(1, 1, 1);

    text(ax, x001 - 0.03, y001 - 0.02, '[001]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(ax, x011 + 0.03, y011 - 0.02, '[011]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(ax, x111, y111 + 0.03, '[111]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    axis(ax, 'equal');
    axis(ax, 'off');
    title(ax, titleStr, 'FontSize', 12);
    hold(ax, 'off');
end

function drawStandardTriangle(ax)
% Draw the IPF-colored standard triangle

    nDiv = 50;
    vertices = [];
    colors = [];
    indexMap = zeros(nDiv+1, nDiv+1);
    idx = 0;

    for i = 0:nDiv
        for j = i:nDiv
            h = i / nDiv;
            k = j / nDiv;
            l = 1;

            [x, y] = direction2stereoSingle(h, k, l);

            idx = idx + 1;
            vertices(idx, :) = [x, y];
            indexMap(i+1, j+1) = idx;

            rgb = calculateIPFColor([h, k, l]);
            colors(idx, :) = rgb;
        end
    end

    % Build faces
    faces = [];
    for i = 0:(nDiv-1)
        for j = i:(nDiv-1)
            idx1 = indexMap(i+1, j+1);
            idx2 = indexMap(i+1, j+2);
            idx3 = indexMap(i+2, j+2);
            faces = [faces; idx1, idx2, idx3];

            if j > i
                idx4 = indexMap(i+2, j+1);
                faces = [faces; idx1, idx3, idx4];
            end
        end
    end

    patch(ax, 'Faces', faces, 'Vertices', vertices, ...
          'FaceVertexCData', colors, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none', ...
          'FaceAlpha', 0.3);
end

function drawTriangleBoundary(ax)
% Draw the triangle boundary lines

    nBoundary = 100;
    t = linspace(0, 1, nBoundary);

    bx1 = zeros(1, nBoundary); by1 = zeros(1, nBoundary);
    bx2 = zeros(1, nBoundary); by2 = zeros(1, nBoundary);
    bx3 = zeros(1, nBoundary); by3 = zeros(1, nBoundary);

    for i = 1:nBoundary
        [bx1(i), by1(i)] = direction2stereoSingle(0, t(i), 1);
        [bx2(i), by2(i)] = direction2stereoSingle(t(i), 1, 1);
        [bx3(i), by3(i)] = direction2stereoSingle(1-t(i), 1-t(i), 1);
    end

    plot(ax, bx1, by1, 'k-', 'LineWidth', 1.5);
    plot(ax, bx2, by2, 'k-', 'LineWidth', 1.5);
    plot(ax, bx3, by3, 'k-', 'LineWidth', 1.5);
end

function [x, y] = direction2stereo(dirs)
% Convert multiple directions to stereographic coordinates

    n = size(dirs, 1);
    x = zeros(n, 1);
    y = zeros(n, 1);

    for i = 1:n
        [x(i), y(i)] = direction2stereoSingle(dirs(i,1), dirs(i,2), dirs(i,3));
    end
end

function [x, y] = direction2stereoSingle(h, k, l)
% Convert single direction to stereographic coordinates

    norm_d = sqrt(h^2 + k^2 + l^2);
    if norm_d > 0
        h = h / norm_d;
        k = k / norm_d;
        l = l / norm_d;
    end

    theta = acos(l);
    if (h^2 + k^2) > 0
        phi = atan2(k, h);
    else
        phi = 0;
    end

    r = tan(theta / 2);
    x_stereo = r * cos(phi);
    y_stereo = r * sin(phi);

    x = y_stereo;
    y = x_stereo;
end
