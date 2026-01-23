function h = plotOrientationCube(orientationInput, varargin)
% PLOTORIENTATIONCUBE Plot a cube representing crystallographic orientation in 3D
%
%   H = PLOTORIENTATIONCUBE(ORIENTATIONINPUT) plots an orientation cube
%   at the origin with default size.
%
%   H = PLOTORIENTATIONCUBE(ORIENTATIONINPUT, 'Name', Value, ...) specifies
%   additional options using name-value pairs.
%
%   Input:
%       orientationInput - Either:
%                          - 3x3 rotation matrix describing the orientation
%                            of the crystal coordinate system relative to the
%                            world coordinate system (crystal to world transform)
%                          - crystalOrientation object
%
%   Optional Name-Value Pairs:
%       'position'   - [x, y, z] position for cube center (default: [0, 0, 0])
%                      If a crystalOrientation object is passed and has a non-zero
%                      position, that position is used unless overridden.
%       'mesh'       - Mesh structure with 'vertices' field. If provided,
%                      cube is placed at center of gravity (mean of vertices).
%                      Overrides 'position' if both are specified.
%       'size'       - Cube size (default: 1)
%       'axes'       - Axes handle to plot into (default: current axes or new figure)
%       'showAxes'   - Logical, show crystal axis arrows and labels (default: true)
%       'faceAlpha'  - Face transparency 0-1 (default: 0.8)
%       'edgeColor'  - Edge color (default: 'k')
%       'lineWidth'  - Edge line width (default: 1.5)
%
%   Output:
%       h - Structure containing handles to graphical objects:
%           h.patch  - Handle to cube patch object
%           h.arrows - Handles to axis arrows (quiver3)
%           h.labels - Handles to axis labels (text)
%
%   Example:
%       % Plot cube at origin with identity orientation
%       R = eye(3);
%       plotOrientationCube(R);
%
%       % Using crystalOrientation object
%       ori = crystalOrientation(eye(3), 'cubic', 'position', [10, 20, 30]);
%       plotOrientationCube(ori);  % Uses position from object
%
%       % Plot cube at specific position with 45° rotation around Z
%       theta = pi/4;
%       R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
%       plotOrientationCube(R, 'position', [10, 20, 30], 'size', 5);
%
%       % Plot cube at center of a mesh
%       mesh.vertices = rand(100, 3) * 50;  % Random vertices
%       plotOrientationCube(R, 'mesh', mesh, 'size', 3);
%
%   See also: crystalOrientation, stereoProj

    % Determine input type and extract rotation matrix
    if isa(orientationInput, 'crystalOrientation')
        R = orientationInput.getRotationMatrix();
        defaultPosition = orientationInput.position;
    elseif isnumeric(orientationInput) && all(size(orientationInput) == [3, 3])
        R = orientationInput;
        defaultPosition = [0, 0, 0];
    else
        error('Input must be a 3x3 rotation matrix or a crystalOrientation object');
    end

    % Parse inputs
    p = inputParser;
    addRequired(p, 'orientationInput');
    addParameter(p, 'position', defaultPosition, @(x) isnumeric(x) && numel(x)==3);
    addParameter(p, 'mesh', [], @(x) isempty(x) || (isstruct(x) && isfield(x, 'vertices')));
    addParameter(p, 'size', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'axes', [], @(x) isempty(x) || ishandle(x));
    addParameter(p, 'showAxes', true, @islogical);
    addParameter(p, 'faceAlpha', 0.8, @(x) isnumeric(x) && x >= 0 && x <= 1);
    addParameter(p, 'edgeColor', 'k', @(x) ischar(x) || (isnumeric(x) && numel(x)==3));
    addParameter(p, 'lineWidth', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(p, orientationInput, varargin{:});
    position = p.Results.position(:)';
    meshData = p.Results.mesh;
    cubeSize = p.Results.size;
    ax = p.Results.axes;
    showAxes = p.Results.showAxes;
    faceAlpha = p.Results.faceAlpha;
    edgeColor = p.Results.edgeColor;
    lineWidth = p.Results.lineWidth;

    % If mesh is provided, calculate center of gravity from vertices
    if ~isempty(meshData)
        vertices = meshData.vertices;
        position = mean(vertices, 1);
    end

    % Get or create axes
    if isempty(ax)
        ax = gca;
    end

    holdState = ishold(ax);
    hold(ax, 'on');

    % Define unit cube vertices (centered at origin)
    s = 0.5 * cubeSize;  % half-size
    cubeVerts = [
        -1 -1 -1;
         1 -1 -1;
         1  1 -1;
        -1  1 -1;
        -1 -1  1;
         1 -1  1;
         1  1  1;
        -1  1  1
    ] * s;

    % Rotate vertices by crystal orientation
    rotVerts = (R * cubeVerts')';

    % Translate to position
    rotVerts = rotVerts + position;

    % Define faces (6 faces of cube)
    faces = [
        1 2 3 4;  % bottom (-Z)
        5 6 7 8;  % top (+Z)
        1 2 6 5;  % front (-Y)
        3 4 8 7;  % back (+Y)
        1 4 8 5;  % left (-X)
        2 3 7 6   % right (+X)
    ];

    % Face colors based on crystal axis directions
    % Blue for [001] faces, Green for [010] faces, Red for [100] faces
    faceColors = [
        0    0    1;    % bottom (-Z) - Blue
        0    0    1;    % top (+Z) - Blue
        0    0.7  0;    % front (-Y) - Green
        0    0.7  0;    % back (+Y) - Green
        1    0    0;    % left (-X) - Red
        1    0    0     % right (+X) - Red
    ];

    % Draw the cube faces
    hPatch = patch(ax, 'Vertices', rotVerts, 'Faces', faces, ...
          'FaceVertexCData', faceColors, ...
          'FaceColor', 'flat', ...
          'EdgeColor', edgeColor, ...
          'LineWidth', lineWidth, ...
          'FaceAlpha', faceAlpha);

    % Initialize output handles
    h.patch = hPatch;
    h.arrows = [];
    h.labels = [];

    % Add crystal axis arrows and labels if requested
    if showAxes
        axisLength = 0.9 * cubeSize;
        axisColors = [1 0 0; 0 0.7 0; 0 0 1];  % RGB for X, Y, Z
        axisLabels = {'[100]', '[010]', '[001]'};

        hArrows = gobjects(3, 1);
        hLabels = gobjects(3, 1);

        for i = 1:3
            dir = zeros(1, 3);
            dir(i) = 1;
            rotDir = (R * dir')';

            % Draw axis arrow from cube center
            hArrows(i) = quiver3(ax, position(1), position(2), position(3), ...
                    rotDir(1)*axisLength, rotDir(2)*axisLength, rotDir(3)*axisLength, ...
                    'Color', axisColors(i,:), 'LineWidth', 2, 'MaxHeadSize', 0.3, ...
                    'AutoScale', 'off');

            % Add label at arrow tip
            labelPos = position + rotDir * axisLength * 1.15;
            hLabels(i) = text(ax, labelPos(1), labelPos(2), labelPos(3), ...
                 axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold', ...
                 'Color', axisColors(i,:), 'HorizontalAlignment', 'center');
        end

        h.arrows = hArrows;
        h.labels = hLabels;
    end

    % Restore hold state
    if ~holdState
        hold(ax, 'off');
    end

    % Set axis properties for 3D viewing if this is a new plot
    axis(ax, 'equal');
    axis(ax, 'vis3d');

    % Add rotation capability
    rotate3d(ax, 'on');
end
