function [fvColored, ax] = ipfMesh(fv, orientationInput, varargin)
% IPFMESH Colors a mesh based on crystallographic surface normal orientation
%
%   [FVCOLORED, AX] = IPFMESH(FV, ORIENTATIONINPUT) takes a mesh structure
%   and a crystal orientation, calculates the crystallographic direction of
%   each vertex normal, and displays the mesh colored according to the IPF
%   coloring scheme for cubic crystals.
%
%   Input:
%       fv - Mesh structure with fields:
%            fv.vertices - Nx3 array of vertex coordinates
%            fv.faces    - Mx3 array of face indices (triangles)
%
%       orientationInput - Either:
%                          - 3x3 rotation matrix describing the orientation
%                            of the crystal coordinate system relative to the
%                            world coordinate system. The matrix transforms
%                            vectors from crystal to world coordinates:
%                            world_vector = R * crystal_vector
%                          - crystalOrientation object
%
%   Optional Name-Value Pairs:
%       'plot'       - logical, if true displays the mesh (default: true)
%       'showLegend' - logical, if true shows IPF triangle legend (default: true)
%       'lighting'   - logical, if true adds lighting (default: false)
%
%   Output:
%       fvColored - Structure with fields:
%                   .vertices - same as input
%                   .faces    - same as input
%                   .faceVertexCData - Nx3 RGB colors for each vertex
%                   .vertexNormals   - Nx3 vertex normals in world coords
%                   .crystalNormals  - Nx3 vertex normals in crystal coords
%
%       ax - Axes handle (if plot is true)
%
%   Example:
%       % Create a sphere mesh
%       [X,Y,Z] = sphere(50);
%       fv = surf2patch(X,Y,Z,'triangles');
%
%       % Define crystal orientation (identity = crystal aligned with world)
%       R = eye(3);
%
%       % Color and display
%       [fvColored, ax] = ipfMesh(fv, R);
%
%       % Using crystalOrientation object
%       ori = crystalOrientation(eye(3), 'cubic');
%       [fvColored, ax] = ipfMesh(fv, ori);
%
%   See also: crystalOrientation, ipfColor, patch, surf2patch

    % Determine input type and extract rotation matrix
    if isa(orientationInput, 'crystalOrientation')
        R = orientationInput.getRotationMatrix();
    elseif isnumeric(orientationInput) && all(size(orientationInput) == [3, 3])
        R = orientationInput;
    else
        error('orientationInput must be a 3x3 rotation matrix or a crystalOrientation object');
    end

    % Parse inputs
    p = inputParser;
    addRequired(p, 'fv', @(x) isstruct(x) && isfield(x,'vertices') && isfield(x,'faces'));
    addRequired(p, 'orientationInput');
    addParameter(p, 'plot', true, @islogical);
    addParameter(p, 'showLegend', true, @islogical);
    addParameter(p, 'lighting', false, @islogical);
    parse(p, fv, orientationInput, varargin{:});

    doPlot = p.Results.plot;
    showLegend = p.Results.showLegend;
    useLighting = p.Results.lighting;

    % Extract mesh data
    vertices = fv.vertices;
    faces = fv.faces;
    nVertices = size(vertices, 1);

    % Calculate vertex normals
    vertexNormals = calculateVertexNormals(vertices, faces);

    % Transform normals from world to crystal coordinates
    % If R transforms crystal->world, then R' (transpose) transforms world->crystal
    R_world2crystal = R';
    crystalNormals = (R_world2crystal * vertexNormals')';

    % Calculate IPF color for each vertex normal
    vertexColors = zeros(nVertices, 3);
    for i = 1:nVertices
        vertexColors(i, :) = ipfColor(crystalNormals(i, :));
    end

    % Create output structure
    fvColored.vertices = vertices;
    fvColored.faces = faces;
    fvColored.faceVertexCData = vertexColors;
    fvColored.vertexNormals = vertexNormals;
    fvColored.crystalNormals = crystalNormals;

    % Plot if requested
    if doPlot
        figure('Color', 'w');

        if showLegend
            % Create subplot for mesh and legend
            ax = subplot(1, 4, [1 2 3]);
        else
            ax = axes;
        end

        % Display the mesh with vertex colors
        patch(ax, 'Faces', faces, 'Vertices', vertices, ...
              'FaceVertexCData', vertexColors, ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');

        axis(ax, 'equal');
        xlabel(ax, 'X');
        ylabel(ax, 'Y');
        zlabel(ax, 'Z');
        title(ax, 'IPF Colored Mesh');
        view(ax, 3);

        if useLighting
            camlight('headlight');
            material('dull');
        end

        % Add IPF legend
        if showLegend
            axLegend = subplot(1, 4, 4);
            plotIPFLegend(axLegend);
        end
    else
        ax = [];
    end
end

function normals = calculateVertexNormals(vertices, faces)
% CALCULATEVERTEXNORMALS Compute vertex normals by averaging face normals
%
%   NORMALS = CALCULATEVERTEXNORMALS(VERTICES, FACES) computes the normal
%   at each vertex by averaging the normals of all faces sharing that vertex.

    nVertices = size(vertices, 1);
    nFaces = size(faces, 1);

    % Initialize vertex normals
    normals = zeros(nVertices, 3);

    % Calculate face normals and accumulate at vertices
    for i = 1:nFaces
        % Get vertex indices for this face
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);

        % Get vertex positions
        p1 = vertices(v1, :);
        p2 = vertices(v2, :);
        p3 = vertices(v3, :);

        % Calculate face normal using cross product
        edge1 = p2 - p1;
        edge2 = p3 - p1;
        faceNormal = cross(edge1, edge2);

        % Accumulate to vertex normals (weighted by face area, implicit in cross product magnitude)
        normals(v1, :) = normals(v1, :) + faceNormal;
        normals(v2, :) = normals(v2, :) + faceNormal;
        normals(v3, :) = normals(v3, :) + faceNormal;
    end

    % Normalize vertex normals
    for i = 1:nVertices
        n = norm(normals(i, :));
        if n > 0
            normals(i, :) = normals(i, :) / n;
        end
    end
end

function plotIPFLegend(ax)
% PLOTIPFLEGEND Plot the IPF standard triangle as a legend

    axes(ax);
    hold(ax, 'on');

    % Create triangular mesh for continuous coloring
    nDiv = 50;

    vertices = [];
    colors = [];
    idx = 0;
    indexMap = zeros(nDiv+1, nDiv+1);

    for i = 0:nDiv
        for j = i:nDiv
            h = i / nDiv;
            k = j / nDiv;
            l = 1;

            [x, y] = legendCoords(h, k, l);

            idx = idx + 1;
            vertices(idx, :) = [x, y];
            indexMap(i+1, j+1) = idx;

            rgb = ipfColor([h, k, l]);
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

    % Plot mesh
    patch(ax, 'Faces', faces, 'Vertices', vertices, ...
          'FaceVertexCData', colors, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    % Plot boundary
    nBoundary = 100;
    t = linspace(0, 1, nBoundary);

    bx1 = zeros(1, nBoundary); by1 = zeros(1, nBoundary);
    bx2 = zeros(1, nBoundary); by2 = zeros(1, nBoundary);
    bx3 = zeros(1, nBoundary); by3 = zeros(1, nBoundary);

    for i = 1:nBoundary
        [bx1(i), by1(i)] = legendCoords(0, t(i), 1);
        [bx2(i), by2(i)] = legendCoords(t(i), 1, 1);
        [bx3(i), by3(i)] = legendCoords(1-t(i), 1-t(i), 1);
    end

    plot(ax, bx1, by1, 'k-', 'LineWidth', 1.5);
    plot(ax, bx2, by2, 'k-', 'LineWidth', 1.5);
    plot(ax, bx3, by3, 'k-', 'LineWidth', 1.5);

    % Labels
    [x001, y001] = legendCoords(0, 0, 1);
    [x011, y011] = legendCoords(0, 1, 1);
    [x111, y111] = legendCoords(1, 1, 1);

    text(ax, x001 - 0.03, y001 - 0.02, '[001]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(ax, x011 + 0.03, y011 - 0.02, '[011]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(ax, x111, y111 + 0.03, '[111]', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    axis(ax, 'equal');
    axis(ax, 'off');
    title(ax, 'IPF Legend', 'FontSize', 12);
    hold(ax, 'off');
end

function [x, y] = legendCoords(h, k, l)
% Convert crystallographic direction to legend plot coordinates
% Uses stereographic projection with [001] bottom-left, [011] bottom-right, [111] top

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

    % Swap for desired orientation
    x = y_stereo;
    y = x_stereo;
end
