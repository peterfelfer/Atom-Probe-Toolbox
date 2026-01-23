function [rgb, varargout] = ipfColor(directionOrOrientation, varargin)
% IPFCOLOR Computes IPF color for crystallographic direction (cubic system)
%
%   RGB = IPFCOLOR(DIRECTION) returns the RGB color for a crystallographic
%   direction in a cubic crystal system using the standard IPF coloring
%   scheme with the [001]-[011]-[111] standard triangle.
%
%   RGB = IPFCOLOR(CRYSTALORIENTATION, WORLDDIRECTION) returns the IPF color
%   for a world direction transformed into the crystal frame using the given
%   crystalOrientation object. This is useful for coloring surfaces based on
%   their normal directions.
%
%   [RGB, AX] = IPFCOLOR(..., 'plot', true) additionally plots the
%   standard triangle with the direction marked and returns the axes handle.
%
%   Input:
%       direction - 3-component vector [h k l] specifying crystallographic direction
%
%       OR
%
%       crystalOrientation - crystalOrientation object
%       worldDirection     - 3-component vector in world coordinates
%
%   Optional Name-Value Pairs:
%       'plot'      - logical, if true creates a plot (default: false)
%       'MarkerSize'- size of the marker in the plot (default: 100)
%
%   Output:
%       rgb - 1x3 RGB color vector (values 0-1)
%       ax  - axes handle (only if 'plot' is true)
%
%   Color scheme:
%       [001] = Blue  (bottom-left)
%       [011] = Green (bottom-right)
%       [111] = Red   (top)
%
%   Example:
%       % Direct crystallographic direction
%       rgb = ipfColor([1 2 3]);
%       [rgb, ax] = ipfColor([1 1 0], 'plot', true);
%
%       % Using crystalOrientation with world direction
%       ori = crystalOrientation(eye(3), 'cubic');
%       worldNormal = [0 0 1];  % World Z direction
%       rgb = ipfColor(ori, worldNormal);
%
%   See also: crystalOrientation, ipfMesh, plotStandardTriangle

    % Determine input type
    if isa(directionOrOrientation, 'crystalOrientation')
        % Two-argument form: crystalOrientation + world direction
        if nargin < 2
            error('When passing a crystalOrientation, a world direction must also be provided');
        end
        worldDirection = varargin{1};
        if ~isnumeric(worldDirection) || numel(worldDirection) ~= 3
            error('World direction must be a 3-component vector');
        end

        % Transform world direction to crystal coordinates
        R = directionOrOrientation.getRotationMatrix();
        direction = (R' * worldDirection(:))';  % R' transforms world->crystal

        % Remove the world direction from varargin for parsing
        varargin = varargin(2:end);
    elseif isnumeric(directionOrOrientation) && numel(directionOrOrientation) == 3
        % Single-argument form: direct crystallographic direction
        direction = directionOrOrientation;
    else
        error('First argument must be a 3-component direction vector or a crystalOrientation object');
    end

    % Parse remaining inputs
    p = inputParser;
    addParameter(p, 'plot', false, @islogical);
    addParameter(p, 'MarkerSize', 100, @isnumeric);
    parse(p, varargin{:});

    doPlot = p.Results.plot;
    markerSize = p.Results.MarkerSize;

    % Ensure direction is a row vector
    direction = direction(:)';

    % Normalize direction
    direction = direction / norm(direction);

    % Reduce to fundamental zone (standard triangle) using cubic symmetry
    % Take absolute values and sort to get equivalent direction in standard triangle
    dirAbs = abs(direction);
    dirSorted = sort(dirAbs, 'ascend');

    % Now dirSorted = [min, mid, max] which maps to the standard triangle
    % with vertices at [001], [011], [111]
    h = dirSorted(1);  % smallest component
    k = dirSorted(2);  % middle component
    l = dirSorted(3);  % largest component

    % Calculate RGB using the standard IPF coloring scheme
    rgb = calculateIPFColor(h, k, l);

    % Ensure RGB values are in valid range
    rgb = max(0, min(1, rgb));

    % Plot if requested
    if doPlot
        ax = plotStandardTriangle();
        hold(ax, 'on');

        % Get stereographic projection coordinates
        [x, y] = direction2stereo(h, k, l);

        scatter(ax, x, y, markerSize, rgb, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold(ax, 'off');
        varargout{1} = ax;
    elseif nargout > 1
        varargout{1} = [];
    end
end

function rgb = calculateIPFColor(h, k, l)
% Calculate IPF color based on position in standard triangle
% Using the scheme: [001]=Blue, [011]=Green, [111]=Red

    % Normalize so that l (largest component) = 1
    if l > 0
        h = h / l;
        k = k / l;
        l = 1;
    else
        rgb = [0 0 0];
        return;
    end

    % Calculate color components based on proximity to vertices
    % [001]: h=0, k=0, l=1 -> Blue
    % [011]: h=0, k=1, l=1 -> Green
    % [111]: h=1, k=1, l=1 -> Red

    % The color is determined by the relative magnitudes
    % Red component: proportional to h (distance toward [111])
    % Green component: proportional to k-h (distance toward [011])
    % Blue component: proportional to 1-k (distance toward [001])

    red = h;
    green = k - h;
    blue = 1 - k;

    % Normalize to maximum = 1 for vivid colors
    maxVal = max([red, green, blue]);
    if maxVal > 0
        rgb = [red, green, blue] / maxVal;
    else
        rgb = [0 0 1];  % Default to blue for [001]
    end
end

function [x, y] = direction2stereo(h, k, l)
% Convert crystallographic direction (in fundamental zone) to stereographic
% projection coordinates with orientation:
% [001] bottom-left, [011] bottom-right, [111] top

    % Normalize
    norm_d = sqrt(h^2 + k^2 + l^2);
    h = h / norm_d;
    k = k / norm_d;
    l = l / norm_d;

    % Stereographic projection
    % theta = polar angle from [001]
    % phi = azimuthal angle in h-k plane
    theta = acos(l);

    if (h^2 + k^2) > 0
        phi = atan2(k, h);
    else
        phi = 0;
    end

    r = tan(theta / 2);
    x_stereo = r * cos(phi);
    y_stereo = r * sin(phi);

    % Rotate to desired orientation: swap x and y
    % This puts [001] at origin (bottom-left), [011] to the right, [111] at top
    x = y_stereo;
    y = x_stereo;
end

function ax = plotStandardTriangle()
% PLOTSTANDARDTRIANGLE Creates a plot of the standard stereographic triangle
% with continuous IPF coloring using a triangular mesh.
% Orientation: [001] bottom-left, [011] bottom-right, [111] top

    % Create figure
    figure('Color', 'w');
    ax = axes;
    hold(ax, 'on');

    % Create triangular mesh for continuous coloring
    % Sample the fundamental zone and create mesh vertices
    nDiv = 80;  % Number of divisions for smooth coloring

    % Generate vertices by sampling the standard triangle
    % Parameterize using (h, k) where 0 <= h <= k <= 1 (with l=1)
    vertices = [];
    colors = [];
    hkValues = [];  % Store h,k values for face generation

    % Create a grid in the (h,k) parameter space
    % The standard triangle in (h,k) space (with l=1 normalized) is:
    % 0 <= h <= k <= 1

    idx = 0;
    indexMap = zeros(nDiv+1, nDiv+1);  % Map from (i,j) to vertex index

    for i = 0:nDiv
        for j = i:nDiv
            h = i / nDiv;
            k = j / nDiv;
            l = 1;

            % Get stereographic coordinates
            [x, y] = direction2stereo(h, k, l);

            idx = idx + 1;
            vertices(idx, :) = [x, y];
            indexMap(i+1, j+1) = idx;

            % Get color
            rgb = calculateIPFColor(h, k, l);
            colors(idx, :) = rgb;
            hkValues(idx, :) = [h, k];
        end
    end

    % Build triangular faces
    faces = [];

    for i = 0:(nDiv-1)
        for j = i:(nDiv-1)
            % Get vertex indices
            idx1 = indexMap(i+1, j+1);      % (i, j)
            idx2 = indexMap(i+1, j+2);      % (i, j+1)
            idx3 = indexMap(i+2, j+2);      % (i+1, j+1)

            % Lower triangle
            faces = [faces; idx1, idx2, idx3];

            % Upper triangle (if valid)
            if j > i
                idx4 = indexMap(i+2, j+1);  % (i+1, j)
                faces = [faces; idx1, idx3, idx4];
            end
        end
    end

    % Plot the mesh with interpolated colors
    patch(ax, 'Faces', faces, 'Vertices', vertices, ...
          'FaceVertexCData', colors, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    % Plot triangle boundary (the curved edges of stereographic projection)
    nBoundary = 100;

    % Edge [001] to [011]: h=0, k varies from 0 to 1
    t = linspace(0, 1, nBoundary);
    bx1 = zeros(1, nBoundary);
    by1 = zeros(1, nBoundary);
    for i = 1:nBoundary
        [bx1(i), by1(i)] = direction2stereo(0, t(i), 1);
    end
    plot(ax, bx1, by1, 'k-', 'LineWidth', 1.5);

    % Edge [011] to [111]: h varies from 0 to 1, k=1
    bx2 = zeros(1, nBoundary);
    by2 = zeros(1, nBoundary);
    for i = 1:nBoundary
        [bx2(i), by2(i)] = direction2stereo(t(i), 1, 1);
    end
    plot(ax, bx2, by2, 'k-', 'LineWidth', 1.5);

    % Edge [111] to [001]: h=k, both vary from 1 to 0
    bx3 = zeros(1, nBoundary);
    by3 = zeros(1, nBoundary);
    for i = 1:nBoundary
        val = 1 - t(i);
        [bx3(i), by3(i)] = direction2stereo(val, val, 1);
    end
    plot(ax, bx3, by3, 'k-', 'LineWidth', 1.5);

    % Get vertex positions for labels
    [x001, y001] = direction2stereo(0, 0, 1);
    [x011, y011] = direction2stereo(0, 1, 1);
    [x111, y111] = direction2stereo(1, 1, 1);

    % Add vertex labels
    text(ax, x001 - 0.03, y001 - 0.02, '[001]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(ax, x011 + 0.03, y011 - 0.02, '[011]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(ax, x111, y111 + 0.03, '[111]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Format axes
    axis(ax, 'equal');
    axis(ax, 'off');
    title(ax, 'IPF Standard Triangle (Cubic)', 'FontSize', 14);

    hold(ax, 'off');
end
