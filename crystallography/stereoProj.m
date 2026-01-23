function [ori, ax] = stereoProj(orientationInput, varargin)
% STEREOPROJ Creates a stereographic projection of crystallographic poles
%
%   [ORI, AX] = STEREOPROJ(CRYSTALORIENTATION) creates a stereographic projection
%   showing the poles of common crystallographic plane families for a crystal
%   with the given orientation.
%
%   The projection uses the APT (atom probe tomography) convention:
%   - View is looking down the +Z axis (analysis direction)
%   - +Z direction is at the center of the projection
%   - X-Y plane forms the projection plane (+X right, +Y up)
%
%   Input:
%       orientationInput - Either:
%                          - 3x3 rotation matrix describing the orientation
%                            of the crystal coordinate system relative to the
%                            world coordinate system (crystal to world transform)
%                          - crystalOrientation object
%
%   Optional Name-Value Pairs:
%       'crystalSystem' - Crystal system name (default: 'cubic')
%                         Only used if input is a rotation matrix.
%                         Options: 'cubic', 'hexagonal', 'tetragonal',
%                                  'orthorhombic', 'monoclinic', 'triclinic'
%       'poles'       - Cell array of pole families to display
%                       (default: {'001', '011', '111'})
%       'showWulff'   - logical, show Wulff net grid (default: true)
%       'showLabels'  - logical, show pole labels (default: true)
%       'upperOnly'   - logical, show only upper hemisphere (default: true)
%       'markerSize'  - size of pole markers (default: 80)
%       'showCube'    - logical, show orientation cube (default: true)
%
%   Output:
%       ori - crystalOrientation object with current orientation
%       ax  - Axes handle
%
%   Interactive Use:
%       After interactive manipulation, retrieve the current orientation:
%           ori = getCurrentOrientation(figHandle)
%
%   Pole family notation:
%       '001' - {100} family: (100), (010), (001), etc.
%       '011' - {110} family: (110), (101), (011), etc.
%       '111' - {111} family: (111), etc.
%       '012' - {012} family: (012), (021), (102), etc.
%       '112' - {112} family
%       '122' - {122} family
%       '113' - {113} family
%
%   Example:
%       % Identity orientation (crystal aligned with world)
%       R = eye(3);
%       [ori, ax] = stereoProj(R);
%
%       % Using crystalOrientation object
%       ori = crystalOrientation(eye(3), 'cubic');
%       [ori, ax] = stereoProj(ori);
%
%       % Show specific pole families with hexagonal system
%       R = eye(3);
%       stereoProj(R, 'crystalSystem', 'hexagonal', 'poles', {'001', '011', '111'});
%
%       % Get orientation after interactive manipulation
%       [~, ax] = stereoProj(eye(3));
%       % ... user drags poles interactively ...
%       ori = getCurrentOrientation(ancestor(ax, 'figure'));
%
%   See also: crystalOrientation, ipfColor, ipfMesh, getCurrentOrientation

    % Determine input type and extract rotation matrix
    if isa(orientationInput, 'crystalOrientation')
        R = orientationInput.getRotationMatrix();
        inputCrystalSystem = orientationInput.crystalSystem;
        inputOriObject = orientationInput;
    elseif isnumeric(orientationInput) && all(size(orientationInput) == [3, 3])
        R = orientationInput;
        inputCrystalSystem = '';  % Will use parameter or default
        inputOriObject = [];
    else
        error('Input must be a 3x3 rotation matrix or a crystalOrientation object');
    end

    % Parse inputs
    p = inputParser;
    addRequired(p, 'orientationInput');
    addParameter(p, 'crystalSystem', 'cubic', @(x) ischar(x) || isstring(x));
    addParameter(p, 'poles', {'001', '011', '111'}, @iscell);
    addParameter(p, 'showWulff', true, @islogical);
    addParameter(p, 'showLabels', true, @islogical);
    addParameter(p, 'upperOnly', true, @islogical);
    addParameter(p, 'markerSize', 80, @isnumeric);
    addParameter(p, 'showCube', true, @islogical);
    parse(p, orientationInput, varargin{:});

    % Determine crystal system (input object takes precedence over parameter)
    if ~isempty(inputCrystalSystem)
        cSystem = inputCrystalSystem;
    else
        cSystem = char(p.Results.crystalSystem);
    end

    % Create crystalOrientation object if not provided
    if isempty(inputOriObject)
        ori = crystalOrientation(R, cSystem);
    else
        ori = inputOriObject;
    end

    poleFamilies = p.Results.poles;
    showWulff = p.Results.showWulff;
    showLabels = p.Results.showLabels;
    upperOnly = p.Results.upperOnly;
    markerSize = p.Results.markerSize;
    showCube = p.Results.showCube;

    % Create figure
    figure('Color', 'w');

    if showCube
        % Create subplot layout: stereographic projection on left, cube on right
        ax = subplot(1, 5, [1 2 3 4]);
    else
        ax = axes;
    end
    hold(ax, 'on');

    % Draw the Wulff net if requested
    if showWulff
        drawWulffNet(ax);
    else
        % Just draw the outer circle
        theta = linspace(0, 2*pi, 100);
        plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5);
    end

    % Define colors and markers for different pole families
    poleStyles = struct();
    poleStyles.p001 = struct('color', [0 0 1], 'marker', 's');      % Blue squares
    poleStyles.p011 = struct('color', [0 0.7 0], 'marker', 'd');    % Green diamonds
    poleStyles.p111 = struct('color', [1 0 0], 'marker', '^');      % Red triangles
    poleStyles.p012 = struct('color', [1 0.5 0], 'marker', 'o');    % Orange circles
    poleStyles.p112 = struct('color', [0.5 0 0.5], 'marker', 'v');  % Purple down-triangles
    poleStyles.p122 = struct('color', [0 0.5 0.5], 'marker', 'p');  % Teal pentagons
    poleStyles.p113 = struct('color', [0.5 0.5 0], 'marker', 'h');  % Olive hexagons
    poleStyles.p123 = struct('color', [0.3 0.3 0.3], 'marker', 'o'); % Gray circles

    % Plot each pole family
    legendEntries = {};
    legendHandles = [];

    for i = 1:length(poleFamilies)
        family = poleFamilies{i};
        poles = generatePoleFamily(family);

        % Get style for this family
        styleField = ['p', family];
        if isfield(poleStyles, styleField)
            style = poleStyles.(styleField);
        else
            style = struct('color', rand(1,3), 'marker', 'o');
        end

        % Track projected positions to avoid duplicates
        projectedPositions = [];
        projectedLabels = {};

        % Transform poles to world coordinates and project
        for j = 1:size(poles, 1)
            pole = poles(j, :);
            pole = pole / norm(pole);

            % Transform to world coordinates
            worldPole = (R * pole')';

            % Check hemisphere
            if upperOnly && worldPole(3) < 0
                worldPole = -worldPole;  % Project antipodal point
            end

            if worldPole(3) >= 0 || ~upperOnly
                % Stereographic projection
                [x, y] = stereoProject(worldPole);

                % Check if we already have a pole at this position
                isDuplicate = false;
                for k = 1:size(projectedPositions, 1)
                    if norm([x, y] - projectedPositions(k, :)) < 0.01
                        isDuplicate = true;
                        break;
                    end
                end

                % Plot if within bounds and not a duplicate
                if x^2 + y^2 <= 1.01 && ~isDuplicate
                    h = scatter(ax, x, y, markerSize, style.color, 'filled', ...
                               'Marker', style.marker, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

                    % Track this position
                    projectedPositions = [projectedPositions; x, y];

                    % Add label if requested
                    if showLabels
                        labelStr = sprintf('(%d%d%d)', poles(j,1), poles(j,2), poles(j,3));
                        text(ax, x + 0.03, y + 0.03, labelStr, 'FontSize', 8);
                    end
                end
            end
        end

        % Store for legend
        legendHandles(end+1) = scatter(ax, nan, nan, markerSize, style.color, 'filled', ...
                                       'Marker', style.marker, 'MarkerEdgeColor', 'k');
        legendEntries{end+1} = ['{', family(1), family(2), family(3), '}'];
    end

    % Add coordinate markers for APT convention (looking down +z)
    % Center is +z direction, x-y plane is the projection plane
    text(ax, 1.12, 0, '+X', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(ax, -1.12, 0, '-X', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(ax, 0, 1.12, '+Y', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(ax, 0, -1.12, '-Y', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(ax, 0.05, 0.05, '+Z', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    % Format axes
    axis(ax, 'equal');
    xlim(ax, [-1.3, 1.3]);
    ylim(ax, [-1.3, 1.3]);
    axis(ax, 'off');
    title(ax, 'Stereographic Projection', 'FontSize', 14);

    % Add legend
    if ~isempty(legendHandles)
        legend(ax, legendHandles, legendEntries, 'Location', 'eastoutside');
    end

    hold(ax, 'off');

    % Draw orientation cube
    if showCube
        axCube = subplot(1, 5, 5);
        drawOrientationCube(axCube, R);

        % Tag the axes for reliable retrieval in callbacks
        ax.Tag = 'stereoAxis';
        axCube.Tag = 'cubeAxis';

        % Store data for interactive updates
        figHandle = gcf;
        setappdata(figHandle, 'rotationMatrix', R);
        setappdata(figHandle, 'baseRotationMatrix', R);  % Original orientation (doesn't change)
        setappdata(figHandle, 'crystalOrientationObject', ori);
        setappdata(figHandle, 'crystalSystem', cSystem);
        setappdata(figHandle, 'poleFamilies', poleFamilies);
        setappdata(figHandle, 'showWulff', showWulff);
        setappdata(figHandle, 'showLabels', showLabels);
        setappdata(figHandle, 'upperOnly', upperOnly);
        setappdata(figHandle, 'markerSize', markerSize);
        setappdata(figHandle, 'poleStyles', poleStyles);

        % Enable interactive rotation on cube only
        h = rotate3d(figHandle);
        h.Enable = 'on';
        h.ActionPostCallback = @updateAfterRotation;
        setAllowAxesRotate(h, ax, false);      % Disable rotation on stereo axis
        setAllowAxesRotate(h, axCube, true);   % Enable rotation on cube axis

        % Add instruction text for cube
        annotation(figHandle, 'textbox', [0.75, 0.02, 0.25, 0.05], ...
                   'String', 'Drag cube to rotate crystal', ...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                   'FontSize', 9, 'FontAngle', 'italic');

        % Add orientation control panel
        createOrientationPanel(figHandle, R);
    end
end

function createOrientationPanel(figHandle, R)
% Create editable UI controls for orientation

    % Get Euler angles
    [phi1, Phi, phi2] = rotmat2euler(R);

    % Create panel
    hPanel = uipanel(figHandle, 'Title', 'Orientation', ...
                     'Units', 'normalized', ...
                     'Position', [0.01, 0.01, 0.28, 0.22], ...
                     'FontSize', 9);

    % Euler angles section
    uicontrol(hPanel, 'Style', 'text', 'String', 'Euler Angles (Bunge ZXZ):', ...
              'Units', 'normalized', 'Position', [0.02, 0.82, 0.96, 0.15], ...
              'HorizontalAlignment', 'left', 'FontSize', 8);

    % phi1
    uicontrol(hPanel, 'Style', 'text', 'String', 'φ1:', ...
              'Units', 'normalized', 'Position', [0.02, 0.62, 0.15, 0.15], ...
              'HorizontalAlignment', 'right', 'FontSize', 8);
    hPhi1 = uicontrol(hPanel, 'Style', 'edit', 'String', sprintf('%.1f', phi1), ...
              'Units', 'normalized', 'Position', [0.18, 0.64, 0.25, 0.15], ...
              'FontSize', 8, 'Tag', 'phi1Edit', ...
              'Callback', @(src,~) eulerEditCallback(figHandle));

    % Phi
    uicontrol(hPanel, 'Style', 'text', 'String', 'Φ:', ...
              'Units', 'normalized', 'Position', [0.45, 0.62, 0.12, 0.15], ...
              'HorizontalAlignment', 'right', 'FontSize', 8);
    hPhi = uicontrol(hPanel, 'Style', 'edit', 'String', sprintf('%.1f', Phi), ...
              'Units', 'normalized', 'Position', [0.58, 0.64, 0.25, 0.15], ...
              'FontSize', 8, 'Tag', 'PhiEdit', ...
              'Callback', @(src,~) eulerEditCallback(figHandle));

    % phi2
    uicontrol(hPanel, 'Style', 'text', 'String', 'φ2:', ...
              'Units', 'normalized', 'Position', [0.02, 0.42, 0.15, 0.15], ...
              'HorizontalAlignment', 'right', 'FontSize', 8);
    hPhi2 = uicontrol(hPanel, 'Style', 'edit', 'String', sprintf('%.1f', phi2), ...
              'Units', 'normalized', 'Position', [0.18, 0.44, 0.25, 0.15], ...
              'FontSize', 8, 'Tag', 'phi2Edit', ...
              'Callback', @(src,~) eulerEditCallback(figHandle));

    % Degree symbol labels
    uicontrol(hPanel, 'Style', 'text', 'String', '°', ...
              'Units', 'normalized', 'Position', [0.43, 0.62, 0.05, 0.15], ...
              'HorizontalAlignment', 'left', 'FontSize', 8);
    uicontrol(hPanel, 'Style', 'text', 'String', '°', ...
              'Units', 'normalized', 'Position', [0.83, 0.62, 0.05, 0.15], ...
              'HorizontalAlignment', 'left', 'FontSize', 8);
    uicontrol(hPanel, 'Style', 'text', 'String', '°', ...
              'Units', 'normalized', 'Position', [0.43, 0.42, 0.05, 0.15], ...
              'HorizontalAlignment', 'left', 'FontSize', 8);

    % Rotation matrix display (read-only)
    uicontrol(hPanel, 'Style', 'text', 'String', 'Rotation Matrix:', ...
              'Units', 'normalized', 'Position', [0.02, 0.22, 0.96, 0.15], ...
              'HorizontalAlignment', 'left', 'FontSize', 8);

    matStr = sprintf('%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f', ...
                     R(1,1), R(1,2), R(1,3), R(2,1), R(2,2), R(2,3), R(3,1), R(3,2), R(3,3));
    hMat = uicontrol(hPanel, 'Style', 'text', 'String', matStr, ...
              'Units', 'normalized', 'Position', [0.02, 0.02, 0.96, 0.22], ...
              'HorizontalAlignment', 'left', 'FontSize', 7, ...
              'FontName', 'FixedWidth', 'Tag', 'matrixDisplay');

    % Store handles
    setappdata(figHandle, 'hPhi1', hPhi1);
    setappdata(figHandle, 'hPhi', hPhi);
    setappdata(figHandle, 'hPhi2', hPhi2);
    setappdata(figHandle, 'hMatrixDisplay', hMat);
    setappdata(figHandle, 'hOrientationPanel', hPanel);
end

function eulerEditCallback(figHandle)
% Callback when Euler angles are edited

    % Get edit box handles
    hPhi1 = getappdata(figHandle, 'hPhi1');
    hPhi = getappdata(figHandle, 'hPhi');
    hPhi2 = getappdata(figHandle, 'hPhi2');

    % Parse values
    phi1 = str2double(get(hPhi1, 'String'));
    Phi = str2double(get(hPhi, 'String'));
    phi2 = str2double(get(hPhi2, 'String'));

    % Validate
    if isnan(phi1) || isnan(Phi) || isnan(phi2)
        return;  % Invalid input, do nothing
    end

    % Convert Euler angles to rotation matrix
    R_new = euler2rotmat(phi1, Phi, phi2);

    % Update stored orientation
    setappdata(figHandle, 'rotationMatrix', R_new);
    setappdata(figHandle, 'baseRotationMatrix', R_new);
    % Update crystalOrientation object
    cSystem = getappdata(figHandle, 'crystalSystem');
    setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_new, cSystem));

    % Update matrix display
    hMat = getappdata(figHandle, 'hMatrixDisplay');
    matStr = sprintf('%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f', ...
                     R_new(1,1), R_new(1,2), R_new(1,3), ...
                     R_new(2,1), R_new(2,2), R_new(2,3), ...
                     R_new(3,1), R_new(3,2), R_new(3,3));
    set(hMat, 'String', matStr);

    % Reset cube view to default
    axCube = findobj(figHandle, 'Type', 'axes', 'Tag', 'cubeAxis');
    if ~isempty(axCube)
        view(axCube, [-90, 90]);
        cla(axCube);
        drawOrientationCube(axCube, R_new);
    end

    % Update stereographic projection
    updateStereoProjection(figHandle, R_new);
end

function updateStereoProjection(figHandle, R_effective)
% Update the stereographic projection with new orientation

    poleFamilies = getappdata(figHandle, 'poleFamilies');
    showWulff = getappdata(figHandle, 'showWulff');
    showLabels = getappdata(figHandle, 'showLabels');
    upperOnly = getappdata(figHandle, 'upperOnly');
    markerSize = getappdata(figHandle, 'markerSize');

    % Delete old stereo axes and create new one
    axStereo = findobj(figHandle, 'Type', 'axes', 'Tag', 'stereoAxis');
    if ~isempty(axStereo)
        delete(axStereo);
    end

    % Create fresh axes for stereographic projection
    axStereo = subplot(1, 5, [1 2 3 4], 'Parent', figHandle);
    axStereo.Tag = 'stereoAxis';
    hold(axStereo, 'on');

    % Redraw Wulff net
    if showWulff
        drawWulffNetOnAxis(axStereo);
    else
        theta = linspace(0, 2*pi, 100);
        plot(axStereo, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    % Get pole styles
    poleStyles = getappdata(figHandle, 'poleStyles');

    % Legend tracking
    legendHandles = [];
    legendEntries = {};

    for i = 1:length(poleFamilies)
        family = poleFamilies{i};
        poles = generatePoleFamily(family);

        styleField = ['p', family];
        if isfield(poleStyles, styleField)
            style = poleStyles.(styleField);
        else
            style = struct('color', rand(1,3), 'marker', 'o');
        end

        projectedPositions = [];

        for j = 1:size(poles, 1)
            pole = poles(j, :);
            pole = pole / norm(pole);
            worldPole = (R_effective * pole')';

            if upperOnly && worldPole(3) < 0
                worldPole = -worldPole;
            end

            if worldPole(3) >= 0 || ~upperOnly
                [x, y] = stereoProject(worldPole);

                isDuplicate = false;
                for k = 1:size(projectedPositions, 1)
                    if norm([x, y] - projectedPositions(k, :)) < 0.01
                        isDuplicate = true;
                        break;
                    end
                end

                if x^2 + y^2 <= 1.01 && ~isDuplicate
                    scatter(axStereo, x, y, markerSize, style.color, 'filled', ...
                           'Marker', style.marker, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                           'HandleVisibility', 'off');
                    projectedPositions = [projectedPositions; x, y];

                    if showLabels
                        labelStr = sprintf('(%d%d%d)', poles(j,1), poles(j,2), poles(j,3));
                        text(axStereo, x + 0.03, y + 0.03, labelStr, 'FontSize', 8);
                    end
                end
            end
        end

        % Add legend entry
        legendHandles(end+1) = scatter(axStereo, nan, nan, markerSize, style.color, 'filled', ...
                                       'Marker', style.marker, 'MarkerEdgeColor', 'k');
        legendEntries{end+1} = ['{', family(1), family(2), family(3), '}'];
    end

    % Coordinate markers
    text(axStereo, 1.12, 0, '+X', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(axStereo, -1.12, 0, '-X', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(axStereo, 0, 1.12, '+Y', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(axStereo, 0, -1.12, '-Y', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(axStereo, 0.05, 0.05, '+Z', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    % Legend
    if ~isempty(legendHandles)
        legend(axStereo, legendHandles, legendEntries, 'Location', 'eastoutside');
    end

    axis(axStereo, 'equal');
    xlim(axStereo, [-1.3, 1.3]);
    ylim(axStereo, [-1.3, 1.3]);
    axis(axStereo, 'off');
    title(axStereo, 'Stereographic Projection', 'FontSize', 14);

    % Prevent 3D rotation on stereo axes
    h = rotate3d(figHandle);
    setAllowAxesRotate(h, axStereo, false);

    hold(axStereo, 'off');
end

function updateAfterRotation(~, evd)
% Callback function after cube rotation - updates stereographic projection

    axCube = evd.Axes;
    figHandle = ancestor(axCube, 'figure');

    % Get stored data
    R_base = getappdata(figHandle, 'baseRotationMatrix');  % Original crystal orientation

    % Get the view angles from the cube axes
    [az, el] = view(axCube);

    % Use MATLAB's viewmtx to get the transformation matrix
    % viewmtx returns a 4x4 matrix that transforms world to view coordinates
    T = viewmtx(az, el);

    % Extract 3x3 rotation part (upper-left)
    % This transforms world directions to view directions
    % In view coords: +X is right, +Y is up, +Z is out of screen (toward viewer)
    R_view = T(1:3, 1:3);

    % For stereo projection, we need directions in a frame where
    % +Z points toward the viewer (center of projection)
    % The view matrix already does this transformation

    % Effective crystal orientation for stereo projection
    R_effective = R_view * R_base;

    % Store for potential retrieval
    setappdata(figHandle, 'rotationMatrix', R_effective);
    % Update crystalOrientation object
    cSystem = getappdata(figHandle, 'crystalSystem');
    setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_effective, cSystem));

    % Update stereo projection
    updateStereoProjection(figHandle, R_effective);

    % Update orientation UI panel
    updateOrientationUI(figHandle, R_effective);
end

function updateOrientationUI(figHandle, R)
% Update the Euler angle edit boxes and matrix display

    % Compute Euler angles
    [phi1, Phi, phi2] = rotmat2euler(R);

    % Update edit boxes
    hPhi1 = getappdata(figHandle, 'hPhi1');
    hPhi = getappdata(figHandle, 'hPhi');
    hPhi2 = getappdata(figHandle, 'hPhi2');
    hMat = getappdata(figHandle, 'hMatrixDisplay');

    if ~isempty(hPhi1) && isvalid(hPhi1)
        set(hPhi1, 'String', sprintf('%.1f', phi1));
    end
    if ~isempty(hPhi) && isvalid(hPhi)
        set(hPhi, 'String', sprintf('%.1f', Phi));
    end
    if ~isempty(hPhi2) && isvalid(hPhi2)
        set(hPhi2, 'String', sprintf('%.1f', phi2));
    end

    % Update matrix display
    if ~isempty(hMat) && isvalid(hMat)
        matStr = sprintf('%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f', ...
                         R(1,1), R(1,2), R(1,3), ...
                         R(2,1), R(2,2), R(2,3), ...
                         R(3,1), R(3,2), R(3,3));
        set(hMat, 'String', matStr);
    end
end

function drawWulffNetOnAxis(ax)
% Draw Wulff net on specified axis (for callback use)
% All elements have HandleVisibility off to exclude from legend

    theta = linspace(0, 2*pi, 100);
    plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Draw latitude circles with labels
    latitudes = [15, 30, 45, 60, 75];
    for lat = latitudes
        latRad = lat * pi / 180;
        r = tan(latRad / 2);
        if r <= 1
            x = r * cos(theta);
            y = r * sin(theta);
            plot(ax, x, y, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
            % Add label
            text(ax, r + 0.02, 0, sprintf('%d°', lat), 'FontSize', 7, ...
                 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'left', ...
                 'VerticalAlignment', 'middle');
        end
    end

    % Add 90° label at edge
    text(ax, 1.02, 0, '90°', 'FontSize', 7, 'Color', [0.5 0.5 0.5], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

    % Draw longitude lines
    longitudes = 0:15:180-15;
    for lon = longitudes
        lonRad = lon * pi / 180;
        x = [-1, 1] * cos(lonRad);
        y = [-1, 1] * sin(lonRad);
        plot(ax, x, y, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end

    % Add longitude labels at edge (every 30°)
    for lon = 0:30:150
        lonRad = deg2rad(lon);
        x = 1.05 * cos(lonRad);
        y = 1.05 * sin(lonRad);
        text(ax, x, y, sprintf('%d°', lon), 'FontSize', 7, ...
             'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle');
        if lon > 0
            x = 1.05 * cos(lonRad + pi);
            y = 1.05 * sin(lonRad + pi);
            text(ax, x, y, sprintf('%d°', lon + 180), 'FontSize', 7, ...
                 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle');
        end
    end
end

function [x, y] = stereoProject(pole)
% STEREOPROJECT Stereographic projection from upper hemisphere
%   Projects a 3D unit vector onto the equatorial plane
%   Projection is from the south pole (0, 0, -1)

    % For upper hemisphere (z >= 0):
    % x_proj = x / (1 + z)
    % y_proj = y / (1 + z)

    if pole(3) >= 0
        denom = 1 + pole(3);
        if denom > 1e-10
            x = pole(1) / denom;
            y = pole(2) / denom;
        else
            x = pole(1);
            y = pole(2);
        end
    else
        % Lower hemisphere - project from north pole
        denom = 1 - pole(3);
        if denom > 1e-10
            x = pole(1) / denom;
            y = pole(2) / denom;
        else
            x = pole(1);
            y = pole(2);
        end
    end
end

function poles = generatePoleFamily(family)
% GENERATEPOLEFAMILY Generate all poles in a crystallographic family
%   For cubic crystals, generates all symmetrically equivalent poles

    % Parse the family string (e.g., '001', '011', '111', '112')
    h = str2double(family(1));
    k = str2double(family(2));
    l = str2double(family(3));

    % Generate all permutations and sign combinations
    baseIndices = [h, k, l];
    permutations = unique(perms(baseIndices), 'rows');

    poles = [];

    for i = 1:size(permutations, 1)
        p = permutations(i, :);

        % Generate all sign combinations
        for s1 = [1, -1]
            for s2 = [1, -1]
                for s3 = [1, -1]
                    newPole = [s1*p(1), s2*p(2), s3*p(3)];

                    % Skip zero vector
                    if norm(newPole) > 0
                        % Check if this exact pole is already in the list
                        isNew = true;
                        for j = 1:size(poles, 1)
                            if norm(poles(j,:) - newPole) < 1e-10
                                isNew = false;
                                break;
                            end
                        end
                        if isNew
                            poles = [poles; newPole];
                        end
                    end
                end
            end
        end
    end
end

function drawWulffNet(ax)
% DRAWWULFFNET Draw the Wulff net (stereographic grid)
% For projection looking down +z axis (APT convention)

    % Outer circle (equator, z=0 plane)
    theta = linspace(0, 2*pi, 100);
    plot(ax, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    % Draw latitude circles (small circles at constant polar angle from +z)
    latitudes = [15, 30, 45, 60, 75];  % degrees from +z pole
    for lat = latitudes
        drawLatitudeCircle(ax, lat);
        % Add label at the right side of each circle
        r = tan(deg2rad(lat) / 2);
        text(ax, r + 0.02, 0, sprintf('%d°', lat), 'FontSize', 7, ...
             'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'left', ...
             'VerticalAlignment', 'middle');
    end

    % Add 90° label at the edge
    text(ax, 1.02, 0, '90°', 'FontSize', 7, 'Color', [0.5 0.5 0.5], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

    % Draw longitude lines (radial lines from center to edge)
    % These are great circles through the +z pole, which project as straight lines
    longitudes = 0:15:180-15;  % degrees
    for lon = longitudes
        drawLongitudeLine(ax, lon);
    end

    % Add longitude labels at the edge (every 30°)
    for lon = 0:30:150
        lonRad = deg2rad(lon);
        x = 1.05 * cos(lonRad);
        y = 1.05 * sin(lonRad);
        text(ax, x, y, sprintf('%d°', lon), 'FontSize', 7, ...
             'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle');
        % Also label the opposite side (lon + 180)
        if lon > 0
            x = 1.05 * cos(lonRad + pi);
            y = 1.05 * sin(lonRad + pi);
            text(ax, x, y, sprintf('%d°', lon + 180), 'FontSize', 7, ...
                 'Color', [0.5 0.5 0.5], 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle');
        end
    end
end

function drawLatitudeCircle(ax, latDeg)
% Draw a latitude circle at the given polar angle from +z

    latRad = latDeg * pi / 180;

    % In stereographic projection, latitude circles become circles
    % The radius in the projection is tan(lat/2) for upper hemisphere
    r = tan(latRad / 2);

    if r <= 1
        theta = linspace(0, 2*pi, 100);
        x = r * cos(theta);
        y = r * sin(theta);
        plot(ax, x, y, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end

function drawLongitudeLine(ax, lonDeg)
% Draw a longitude line (great circle through the +z pole)
% These project as straight lines through the center

    lonRad = lonDeg * pi / 180;

    % Straight line from edge to edge through center
    x = [-1, 1] * cos(lonRad);
    y = [-1, 1] * sin(lonRad);
    plot(ax, x, y, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

function drawOrientationCube(ax, R)
% DRAWORIENTATIONCUBE Draw a unit cube showing crystal orientation
%   R is the rotation matrix (crystal to world coordinates)

    axes(ax);
    hold(ax, 'on');

    % Define unit cube vertices (centered at origin)
    s = 0.5;  % half-size
    cubeVerts = s * [
        -1 -1 -1;
         1 -1 -1;
         1  1 -1;
        -1  1 -1;
        -1 -1  1;
         1 -1  1;
         1  1  1;
        -1  1  1
    ];

    % Rotate vertices by crystal orientation
    rotVerts = (R * cubeVerts')';

    % Define faces (6 faces of cube)
    faces = [
        1 2 3 4;  % bottom
        5 6 7 8;  % top
        1 2 6 5;  % front
        3 4 8 7;  % back
        1 4 8 5;  % left
        2 3 7 6   % right
    ];

    % Face colors based on crystal axis directions (fixed colors)
    % Blue for [001] faces, Green for [010] faces, Red for [100] faces
    faceColors = [
        0    0    1;    % bottom (-Z) - Blue
        0    0    1;    % top (+Z) - Blue
        0    1    0;    % front (-Y) - Green
        0    1    0;    % back (+Y) - Green
        1    0    0;    % left (-X) - Red
        1    0    0     % right (+X) - Red
    ];

    % Draw the cube faces
    patch(ax, 'Vertices', rotVerts, 'Faces', faces, ...
          'FaceVertexCData', faceColors, ...
          'FaceColor', 'flat', ...
          'EdgeColor', 'k', ...
          'LineWidth', 1.5, ...
          'FaceAlpha', 0.8);

    % Add crystal axis labels at the rotated axis directions
    axisLength = 0.9;
    axisColors = [1 0 0; 0 1 0; 0 0 1];  % RGB for X, Y, Z
    axisLabels = {'[100]', '[010]', '[001]'};

    for i = 1:3
        dir = zeros(1, 3);
        dir(i) = 1;
        rotDir = (R * dir')';

        % Draw axis line
        quiver3(ax, 0, 0, 0, rotDir(1)*axisLength, rotDir(2)*axisLength, rotDir(3)*axisLength, ...
                'Color', axisColors(i,:), 'LineWidth', 2, 'MaxHeadSize', 0.3);

        % Add label
        text(ax, rotDir(1)*axisLength*1.15, rotDir(2)*axisLength*1.15, rotDir(3)*axisLength*1.15, ...
             axisLabels{i}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', axisColors(i,:), ...
             'HorizontalAlignment', 'center');
    end

    % Add world coordinate axes (thin gray lines)
    line(ax, [-0.7 0.7], [0 0], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line(ax, [0 0], [-0.7 0.7], [0 0], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    line(ax, [0 0], [0 0], [-0.7 0.7], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

    % World axis labels
    text(ax, 0.8, 0, 0, 'X', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    text(ax, 0, 0.8, 0, 'Y', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    text(ax, 0, 0, 0.8, 'Z', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);

    % Format axes
    axis(ax, 'equal');
    axis(ax, 'vis3d');
    xlim(ax, [-1.2 1.2]);
    ylim(ax, [-1.2 1.2]);
    zlim(ax, [-1.2 1.2]);

    % View from +Z looking down (matching stereographic projection orientation)
    % azimuth = -90 puts +X to the right, elevation = 90 looks straight down +Z
    view(ax, [-90 90]);
    axis(ax, 'off');
    title(ax, 'Crystal Orientation (view down +Z)', 'FontSize', 12);

    hold(ax, 'off');
end

function [phi1, Phi, phi2] = rotmat2euler(R)
% ROTMAT2EULER Convert rotation matrix to Euler angles (Bunge convention)
%   Bunge convention: ZXZ (phi1 around Z, Phi around X', phi2 around Z'')
%   All angles returned in degrees

    % Handle numerical precision issues
    R = max(-1, min(1, R));

    % Phi is the angle between Z axes (0 to 180 degrees)
    Phi = acosd(R(3,3));

    % Handle gimbal lock cases
    if abs(R(3,3)) > 0.9999
        % Phi ≈ 0 or 180, gimbal lock
        phi1 = atan2d(-R(1,2), R(1,1));
        phi2 = 0;
        if R(3,3) < 0
            Phi = 180;
        else
            Phi = 0;
        end
    else
        % General case
        phi1 = atan2d(R(3,1), -R(3,2));
        phi2 = atan2d(R(1,3), R(2,3));
    end

    % Normalize angles to [0, 360) range
    phi1 = mod(phi1, 360);
    phi2 = mod(phi2, 360);
end

function R = euler2rotmat(phi1, Phi, phi2)
% EULER2ROTMAT Convert Euler angles to rotation matrix (Bunge convention)
%   Bunge convention: ZXZ (phi1 around Z, Phi around X', phi2 around Z'')
%   Input angles in degrees

    % Convert to radians
    phi1 = deg2rad(phi1);
    Phi = deg2rad(Phi);
    phi2 = deg2rad(phi2);

    % Rotation matrices
    c1 = cos(phi1); s1 = sin(phi1);
    c = cos(Phi);   s = sin(Phi);
    c2 = cos(phi2); s2 = sin(phi2);

    % Combined rotation matrix (ZXZ convention)
    R = [c1*c2 - s1*c*s2,  -c1*s2 - s1*c*c2,  s1*s;
         s1*c2 + c1*c*s2,  -s1*s2 + c1*c*c2, -c1*s;
         s*s2,              s*c2,             c];
end

%% =========================================================================
%  DRAG-AND-DROP POLE INTERACTION FUNCTIONS
%  =========================================================================

function dir = stereoInverse(x, y)
% STEREOINVERSE Convert stereographic projection coordinates to 3D direction
%   Given (x, y) on the projection plane, returns the unit vector direction
%   Uses upper hemisphere projection from south pole

    r = sqrt(x^2 + y^2);

    if r < 1e-10
        % At center, direction is +Z
        dir = [0, 0, 1];
        return;
    end

    % Inverse stereographic projection
    % theta = polar angle from +Z axis
    theta = 2 * atan(r);
    phi = atan2(y, x);

    % Convert to Cartesian
    dir = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
end

function R = rodriguesRotation(axis, angle)
% RODRIGUESROTATION Build rotation matrix from axis-angle representation
%   axis  - 3x1 or 1x3 unit vector defining rotation axis
%   angle - rotation angle in radians
%   Returns 3x3 rotation matrix

    % Ensure axis is normalized
    axis = axis(:)' / norm(axis);

    if abs(angle) < 1e-10
        R = eye(3);
        return;
    end

    % Rodrigues' rotation formula
    K = [0, -axis(3), axis(2);
         axis(3), 0, -axis(1);
         -axis(2), axis(1), 0];

    R = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
end

function allowRotate = rotateButtonDownFilter(~, eventData)
% ROTATEBUTTONDOWNFILTER Filter to allow pole clicks to pass through rotate3d
%   Returns false for scatter markers (poles) so their ButtonDownFcn fires
%   Returns true for other objects so rotate3d handles them

    % Get the object that was clicked
    hitObj = eventData.HitObject;

    % Check if it's a scatter marker (pole)
    if isa(hitObj, 'matlab.graphics.chart.primitive.Scatter')
        allowRotate = false;  % Don't let rotate3d handle it - pass to scatter's callback
    else
        allowRotate = true;   % Let rotate3d handle it (e.g., cube rotation)
    end
end

function initDragState(figHandle)
% INITDRAGSTATE Initialize drag state variables for two-pole interaction
%   States: 'idle', 'dragging', 'first_pole_set'

    setappdata(figHandle, 'dragState', 'idle');
    setappdata(figHandle, 'draggedPole', []);           % Crystal direction being dragged
    setappdata(figHandle, 'draggedPoleWorld', []);      % World direction of dragged pole
    setappdata(figHandle, 'firstPoleFixed', false);     % Is first pole locked?
    setappdata(figHandle, 'firstPoleCrystal', []);      % Crystal direction of first pole
    setappdata(figHandle, 'firstPoleWorld', []);        % World direction first pole is fixed to
    setappdata(figHandle, 'ghostMarker', []);           % Handle to ghost marker during drag
    setappdata(figHandle, 'lockedMarker', []);          % Handle to locked pole indicator
    setappdata(figHandle, 'previewMarkers', []);        % Handles to preview markers during drag
    setappdata(figHandle, 'dragStartOrientation', []);  % Orientation at start of drag
    setappdata(figHandle, 'previewOrientation', []);    % Current preview orientation
end

function poleClickCallback(src, ~)
% POLECLICKCALLBACK Handle mouse click on a pole marker

    figHandle = ancestor(src, 'figure');
    dragState = getappdata(figHandle, 'dragState');

    % Get the crystal direction stored in this marker
    crystalDir = get(src, 'UserData');
    if isempty(crystalDir)
        return;
    end

    % Get current orientation and compute world direction
    R = getappdata(figHandle, 'rotationMatrix');
    worldDir = (R * crystalDir(:))';
    worldDir = worldDir / norm(worldDir);

    % Handle based on current state
    firstPoleFixed = getappdata(figHandle, 'firstPoleFixed');

    if firstPoleFixed
        % Check if clicking the same pole that's already fixed
        firstPoleCrystal = getappdata(figHandle, 'firstPoleCrystal');
        if norm(abs(crystalDir(:)) - abs(firstPoleCrystal(:))) < 0.01
            % Clicking the locked pole - ignore or could reset
            return;
        end
    end

    % Store dragged pole info and starting orientation for preview
    setappdata(figHandle, 'draggedPole', crystalDir(:)');
    setappdata(figHandle, 'draggedPoleWorld', worldDir);
    setappdata(figHandle, 'dragStartOrientation', R);
    setappdata(figHandle, 'dragState', 'dragging');

    % Set up motion and release callbacks
    set(figHandle, 'WindowButtonMotionFcn', @figureDragMotion);
    set(figHandle, 'WindowButtonUpFcn', @figureDragRelease);
end

function figureDragMotion(figHandle, ~)
% FIGUREDRAGMOTION Handle mouse motion during pole drag - real-time preview

    dragState = getappdata(figHandle, 'dragState');
    if ~strcmp(dragState, 'dragging')
        return;
    end

    % Get mouse position in stereo axes coordinates
    axStereo = findobj(figHandle, 'Type', 'axes', 'Tag', 'stereoAxis');
    if isempty(axStereo)
        return;
    end

    % Get current point in axes coordinates
    cp = get(axStereo, 'CurrentPoint');
    x = cp(1, 1);
    y = cp(1, 2);

    % Clamp to unit circle
    r = sqrt(x^2 + y^2);
    if r > 1
        x = x / r;
        y = y / r;
    end

    % Get target direction from stereo coordinates
    targetDir = stereoInverse(x, y);

    % Get dragged pole info and current orientation
    crystalDir = getappdata(figHandle, 'draggedPole');
    firstPoleFixed = getappdata(figHandle, 'firstPoleFixed');
    R_base = getappdata(figHandle, 'dragStartOrientation');

    if isempty(R_base)
        R_base = getappdata(figHandle, 'rotationMatrix');
    end

    % Calculate preview rotation
    if ~firstPoleFixed
        % First pole - minimum rotation preview
        R_preview = computePreviewRotation(R_base, crystalDir, targetDir, ...
                                           getappdata(figHandle, 'upperOnly'));
    else
        % Second pole - constrained rotation around first pole
        firstPoleWorld = getappdata(figHandle, 'firstPoleWorld');
        R_preview = computeConstrainedPreviewRotation(R_base, crystalDir, targetDir, ...
                                                       firstPoleWorld, ...
                                                       getappdata(figHandle, 'upperOnly'));
    end

    % Update display with preview orientation
    updatePreviewDisplay(figHandle, axStereo, R_preview);
end

function figureDragRelease(figHandle, ~)
% FIGUREDRAGRELEASE Handle mouse release after pole drag

    dragState = getappdata(figHandle, 'dragState');
    if ~strcmp(dragState, 'dragging')
        return;
    end

    % Clear motion/release callbacks
    set(figHandle, 'WindowButtonMotionFcn', '');
    set(figHandle, 'WindowButtonUpFcn', '');

    % Delete preview markers
    hPreview = getappdata(figHandle, 'previewMarkers');
    if ~isempty(hPreview)
        for i = 1:length(hPreview)
            if isvalid(hPreview(i))
                delete(hPreview(i));
            end
        end
    end
    setappdata(figHandle, 'previewMarkers', []);

    % Get the preview orientation that was computed during drag
    R_preview = getappdata(figHandle, 'previewOrientation');

    % Get dragged pole info
    crystalDir = getappdata(figHandle, 'draggedPole');
    firstPoleFixed = getappdata(figHandle, 'firstPoleFixed');

    if isempty(R_preview)
        % No preview was computed, just reset
        setappdata(figHandle, 'dragState', 'idle');
        return;
    end

    if ~firstPoleFixed
        % First pole drag - use the preview rotation and lock the pole
        % Get target direction from final mouse position
        axStereo = findobj(figHandle, 'Type', 'axes', 'Tag', 'stereoAxis');
        cp = get(axStereo, 'CurrentPoint');
        x = cp(1, 1);
        y = cp(1, 2);
        r = sqrt(x^2 + y^2);
        if r > 1
            x = x / r;
            y = y / r;
        end
        targetDir = stereoInverse(x, y);

        % Store the fixed first pole info
        setappdata(figHandle, 'firstPoleFixed', true);
        setappdata(figHandle, 'firstPoleCrystal', crystalDir(:)');
        setappdata(figHandle, 'firstPoleWorld', targetDir(:)');
        setappdata(figHandle, 'rotationMatrix', R_preview);
        setappdata(figHandle, 'baseRotationMatrix', R_preview);
        % Update crystalOrientation object
        cSystem = getappdata(figHandle, 'crystalSystem');
        setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_preview, cSystem));

        % Update displays
        updateStereoProjectionWithDrag(figHandle, R_preview);
        updateCubeDisplay(figHandle, R_preview);
        updateOrientationUI(figHandle, R_preview);
        updateStatusText(figHandle, 'first_set');
    else
        % Second pole drag - finalize orientation
        setappdata(figHandle, 'rotationMatrix', R_preview);
        setappdata(figHandle, 'baseRotationMatrix', R_preview);
        % Update crystalOrientation object
        cSystem = getappdata(figHandle, 'crystalSystem');
        setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_preview, cSystem));

        % Reset drag state (orientation is now fully determined)
        resetDragState(figHandle);

        % Update displays
        updateStereoProjectionWithDrag(figHandle, R_preview);
        updateCubeDisplay(figHandle, R_preview);
        updateOrientationUI(figHandle, R_preview);
        updateStatusText(figHandle, 'idle');
    end

    setappdata(figHandle, 'dragState', 'idle');
    setappdata(figHandle, 'dragStartOrientation', []);
    setappdata(figHandle, 'previewOrientation', []);
end

function applyFirstPoleRotation(figHandle, crystalDir, targetDir)
% APPLYFIRSTPOLEROTATION Apply minimum rotation to move pole to target

    R_current = getappdata(figHandle, 'rotationMatrix');

    % Current world direction of this pole
    currentWorldDir = (R_current * crystalDir(:))';
    currentWorldDir = currentWorldDir / norm(currentWorldDir);

    % Handle hemisphere (if pole was in lower hemisphere, it was flipped)
    upperOnly = getappdata(figHandle, 'upperOnly');
    if upperOnly && currentWorldDir(3) < 0
        currentWorldDir = -currentWorldDir;
    end

    % Compute rotation from current to target
    targetDir = targetDir(:)';
    dotProd = dot(currentWorldDir, targetDir);
    dotProd = max(-1, min(1, dotProd));  % Clamp for numerical stability

    if abs(dotProd - 1) < 1e-10
        % Already at target, no rotation needed
        R_delta = eye(3);
    elseif abs(dotProd + 1) < 1e-10
        % Opposite direction - rotate 180 around any perpendicular axis
        if abs(currentWorldDir(1)) < 0.9
            perpAxis = cross(currentWorldDir, [1, 0, 0]);
        else
            perpAxis = cross(currentWorldDir, [0, 1, 0]);
        end
        perpAxis = perpAxis / norm(perpAxis);
        R_delta = rodriguesRotation(perpAxis, pi);
    else
        % General case - minimum rotation
        rotAxis = cross(currentWorldDir, targetDir);
        rotAxis = rotAxis / norm(rotAxis);
        rotAngle = acos(dotProd);
        R_delta = rodriguesRotation(rotAxis, rotAngle);
    end

    % Apply rotation to current orientation
    % R_new transforms crystal to new world coordinates
    R_new = R_delta * R_current;

    % Store the fixed first pole info
    setappdata(figHandle, 'firstPoleFixed', true);
    setappdata(figHandle, 'firstPoleCrystal', crystalDir(:)');
    setappdata(figHandle, 'firstPoleWorld', targetDir);
    setappdata(figHandle, 'rotationMatrix', R_new);
    setappdata(figHandle, 'baseRotationMatrix', R_new);
    % Update crystalOrientation object
    cSystem = getappdata(figHandle, 'crystalSystem');
    setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_new, cSystem));

    % Update displays
    updateStereoProjectionWithDrag(figHandle, R_new);
    updateCubeDisplay(figHandle, R_new);
    updateOrientationUI(figHandle, R_new);
    updateStatusText(figHandle, 'first_set');
end

function applySecondPoleRotation(figHandle, crystalDir, targetDir)
% APPLYSECONDPOLEROTATION Apply constrained rotation around first pole

    R_current = getappdata(figHandle, 'rotationMatrix');
    firstPoleWorld = getappdata(figHandle, 'firstPoleWorld');

    % Current world direction of second pole
    currentWorldDir = (R_current * crystalDir(:))';
    currentWorldDir = currentWorldDir / norm(currentWorldDir);

    % Handle hemisphere
    upperOnly = getappdata(figHandle, 'upperOnly');
    if upperOnly && currentWorldDir(3) < 0
        currentWorldDir = -currentWorldDir;
    end

    targetDir = targetDir(:)';
    firstPoleWorld = firstPoleWorld(:)';

    % Project both current and target onto plane perpendicular to first pole
    % This gives us the angle to rotate around the first pole axis

    % Remove component along first pole axis
    currentProj = currentWorldDir - dot(currentWorldDir, firstPoleWorld) * firstPoleWorld;
    targetProj = targetDir - dot(targetDir, firstPoleWorld) * firstPoleWorld;

    currentProjNorm = norm(currentProj);
    targetProjNorm = norm(targetProj);

    if currentProjNorm < 1e-10 || targetProjNorm < 1e-10
        % Second pole is along first pole axis - can't determine rotation
        % Reset state and return
        resetDragState(figHandle);
        return;
    end

    currentProj = currentProj / currentProjNorm;
    targetProj = targetProj / targetProjNorm;

    % Calculate rotation angle around first pole axis
    dotProd = dot(currentProj, targetProj);
    dotProd = max(-1, min(1, dotProd));

    % Determine sign of rotation using cross product
    crossProd = cross(currentProj, targetProj);
    sign = 1;
    if dot(crossProd, firstPoleWorld) < 0
        sign = -1;
    end

    rotAngle = sign * acos(dotProd);

    % Create rotation around first pole axis
    R_delta = rodriguesRotation(firstPoleWorld, rotAngle);

    % Apply rotation
    R_new = R_delta * R_current;

    % Update orientation
    setappdata(figHandle, 'rotationMatrix', R_new);
    setappdata(figHandle, 'baseRotationMatrix', R_new);
    % Update crystalOrientation object
    cSystem = getappdata(figHandle, 'crystalSystem');
    setappdata(figHandle, 'crystalOrientationObject', crystalOrientation(R_new, cSystem));

    % Reset drag state (orientation is now fully determined)
    resetDragState(figHandle);

    % Update displays
    updateStereoProjectionWithDrag(figHandle, R_new);
    updateCubeDisplay(figHandle, R_new);
    updateOrientationUI(figHandle, R_new);
    updateStatusText(figHandle, 'idle');
end

function resetDragState(figHandle)
% RESETDRAGSTATE Reset the drag state to idle

    setappdata(figHandle, 'dragState', 'idle');
    setappdata(figHandle, 'firstPoleFixed', false);
    setappdata(figHandle, 'firstPoleCrystal', []);
    setappdata(figHandle, 'firstPoleWorld', []);
    setappdata(figHandle, 'draggedPole', []);
    setappdata(figHandle, 'draggedPoleWorld', []);

    % Delete locked marker if exists
    hLocked = getappdata(figHandle, 'lockedMarker');
    if ~isempty(hLocked) && isvalid(hLocked)
        delete(hLocked);
    end
    setappdata(figHandle, 'lockedMarker', []);
end

function updateCubeDisplay(figHandle, R)
% UPDATECUBEDISPLAY Update the orientation cube with new rotation

    axCube = findobj(figHandle, 'Type', 'axes', 'Tag', 'cubeAxis');
    if isempty(axCube)
        return;
    end

    % Clear and redraw cube
    cla(axCube);
    drawOrientationCube(axCube, R);

    % Reset view to look down +Z
    view(axCube, [-90, 90]);
end

function updateStatusText(figHandle, state)
% UPDATESTATUSTEXT Update the status text based on current drag state

    hStatus = getappdata(figHandle, 'hStatusText');

    if isempty(hStatus) || ~isvalid(hStatus)
        return;
    end

    switch state
        case 'idle'
            set(hStatus, 'String', 'Drag a pole to set orientation (two-pole mode)');
        case 'first_set'
            set(hStatus, 'String', 'First pole locked. Drag second pole to complete.');
        case 'dragging'
            set(hStatus, 'String', 'Dragging...');
    end
end

function updateStereoProjectionWithDrag(figHandle, R_effective)
% UPDATESTEROPROJECTIONWITHDRAG Update stereo projection with draggable poles
%   Similar to updateStereoProjection but makes poles clickable

    poleFamilies = getappdata(figHandle, 'poleFamilies');
    showWulff = getappdata(figHandle, 'showWulff');
    showLabels = getappdata(figHandle, 'showLabels');
    upperOnly = getappdata(figHandle, 'upperOnly');
    markerSize = getappdata(figHandle, 'markerSize');
    firstPoleFixed = getappdata(figHandle, 'firstPoleFixed');
    firstPoleCrystal = getappdata(figHandle, 'firstPoleCrystal');

    % Delete old stereo axes and create new one
    axStereo = findobj(figHandle, 'Type', 'axes', 'Tag', 'stereoAxis');
    if ~isempty(axStereo)
        delete(axStereo);
    end

    % Create fresh axes for stereographic projection
    axStereo = subplot(1, 5, [1 2 3 4], 'Parent', figHandle);
    axStereo.Tag = 'stereoAxis';
    hold(axStereo, 'on');

    % Redraw Wulff net
    if showWulff
        drawWulffNetOnAxis(axStereo);
    else
        theta = linspace(0, 2*pi, 100);
        plot(axStereo, cos(theta), sin(theta), 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    % Get pole styles
    poleStyles = getappdata(figHandle, 'poleStyles');

    % Legend tracking
    legendHandles = [];
    legendEntries = {};

    for i = 1:length(poleFamilies)
        family = poleFamilies{i};
        poles = generatePoleFamily(family);

        styleField = ['p', family];
        if isfield(poleStyles, styleField)
            style = poleStyles.(styleField);
        else
            style = struct('color', rand(1,3), 'marker', 'o');
        end

        projectedPositions = [];

        for j = 1:size(poles, 1)
            pole = poles(j, :);
            pole = pole / norm(pole);
            worldPole = (R_effective * pole')';

            originalPole = pole;  % Store original crystal direction

            if upperOnly && worldPole(3) < 0
                worldPole = -worldPole;
                originalPole = -pole;  % Flip crystal direction too
            end

            if worldPole(3) >= 0 || ~upperOnly
                [x, y] = stereoProject(worldPole);

                isDuplicate = false;
                for k = 1:size(projectedPositions, 1)
                    if norm([x, y] - projectedPositions(k, :)) < 0.01
                        isDuplicate = true;
                        break;
                    end
                end

                if x^2 + y^2 <= 1.01 && ~isDuplicate
                    % Check if this is the locked first pole
                    isLockedPole = false;
                    if firstPoleFixed && ~isempty(firstPoleCrystal)
                        if norm(abs(originalPole(:)) - abs(firstPoleCrystal(:))) < 0.01
                            isLockedPole = true;
                        end
                    end

                    % Create clickable scatter marker
                    if isLockedPole
                        % Locked pole - different appearance
                        h = scatter(axStereo, x, y, markerSize * 1.5, style.color, 'filled', ...
                                   'Marker', style.marker, 'MarkerEdgeColor', 'k', 'LineWidth', 2, ...
                                   'HandleVisibility', 'off');
                        % Add ring around locked pole
                        hRing = plot(axStereo, x, y, 'o', 'MarkerSize', 18, ...
                                    'MarkerEdgeColor', [1 0.5 0], 'LineWidth', 3, ...
                                    'HandleVisibility', 'off');
                        setappdata(figHandle, 'lockedMarker', hRing);
                    else
                        % Normal draggable pole
                        h = scatter(axStereo, x, y, markerSize, style.color, 'filled', ...
                                   'Marker', style.marker, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5, ...
                                   'HandleVisibility', 'off');
                    end

                    % Store crystal direction and add click callback
                    set(h, 'UserData', originalPole(:)');
                    set(h, 'ButtonDownFcn', @poleClickCallback);

                    projectedPositions = [projectedPositions; x, y];

                    if showLabels
                        labelStr = sprintf('(%d%d%d)', poles(j,1), poles(j,2), poles(j,3));
                        text(axStereo, x + 0.03, y + 0.03, labelStr, 'FontSize', 8);
                    end
                end
            end
        end

        % Add legend entry
        legendHandles(end+1) = scatter(axStereo, nan, nan, markerSize, style.color, 'filled', ...
                                       'Marker', style.marker, 'MarkerEdgeColor', 'k');
        legendEntries{end+1} = ['{', family(1), family(2), family(3), '}'];
    end

    % Coordinate markers
    text(axStereo, 1.12, 0, '+X', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(axStereo, -1.12, 0, '-X', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(axStereo, 0, 1.12, '+Y', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(axStereo, 0, -1.12, '-Y', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(axStereo, 0.05, 0.05, '+Z', 'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

    % Legend
    if ~isempty(legendHandles)
        legend(axStereo, legendHandles, legendEntries, 'Location', 'eastoutside');
    end

    axis(axStereo, 'equal');
    xlim(axStereo, [-1.3, 1.3]);
    ylim(axStereo, [-1.3, 1.3]);
    axis(axStereo, 'off');
    title(axStereo, 'Stereographic Projection', 'FontSize', 14);

    % Prevent 3D rotation on stereo axes
    h = rotate3d(figHandle);
    setAllowAxesRotate(h, axStereo, false);

    hold(axStereo, 'off');
end

function R_preview = computePreviewRotation(R_base, crystalDir, targetDir, upperOnly)
% COMPUTEPREVIEWROTATION Compute preview rotation for first pole drag

    % Current world direction of this pole
    currentWorldDir = (R_base * crystalDir(:))';
    currentWorldDir = currentWorldDir / norm(currentWorldDir);

    % Handle hemisphere
    if upperOnly && currentWorldDir(3) < 0
        currentWorldDir = -currentWorldDir;
    end

    % Compute rotation from current to target
    targetDir = targetDir(:)';
    dotProd = dot(currentWorldDir, targetDir);
    dotProd = max(-1, min(1, dotProd));

    if abs(dotProd - 1) < 1e-10
        R_delta = eye(3);
    elseif abs(dotProd + 1) < 1e-10
        if abs(currentWorldDir(1)) < 0.9
            perpAxis = cross(currentWorldDir, [1, 0, 0]);
        else
            perpAxis = cross(currentWorldDir, [0, 1, 0]);
        end
        perpAxis = perpAxis / norm(perpAxis);
        R_delta = rodriguesRotation(perpAxis, pi);
    else
        rotAxis = cross(currentWorldDir, targetDir);
        rotAxis = rotAxis / norm(rotAxis);
        rotAngle = acos(dotProd);
        R_delta = rodriguesRotation(rotAxis, rotAngle);
    end

    R_preview = R_delta * R_base;
end

function R_preview = computeConstrainedPreviewRotation(R_base, crystalDir, targetDir, firstPoleWorld, upperOnly)
% COMPUTECONSTRAINEDPREVIEWROTATION Compute preview rotation for second pole (constrained)

    % Current world direction of second pole
    currentWorldDir = (R_base * crystalDir(:))';
    currentWorldDir = currentWorldDir / norm(currentWorldDir);

    if upperOnly && currentWorldDir(3) < 0
        currentWorldDir = -currentWorldDir;
    end

    targetDir = targetDir(:)';
    firstPoleWorld = firstPoleWorld(:)';

    % Project both onto plane perpendicular to first pole
    currentProj = currentWorldDir - dot(currentWorldDir, firstPoleWorld) * firstPoleWorld;
    targetProj = targetDir - dot(targetDir, firstPoleWorld) * firstPoleWorld;

    currentProjNorm = norm(currentProj);
    targetProjNorm = norm(targetProj);

    if currentProjNorm < 1e-10 || targetProjNorm < 1e-10
        R_preview = R_base;
        return;
    end

    currentProj = currentProj / currentProjNorm;
    targetProj = targetProj / targetProjNorm;

    dotProd = dot(currentProj, targetProj);
    dotProd = max(-1, min(1, dotProd));

    crossProd = cross(currentProj, targetProj);
    rotSign = 1;
    if dot(crossProd, firstPoleWorld) < 0
        rotSign = -1;
    end

    rotAngle = rotSign * acos(dotProd);
    R_delta = rodriguesRotation(firstPoleWorld, rotAngle);
    R_preview = R_delta * R_base;
end

function updatePreviewDisplay(figHandle, axStereo, R_preview)
% UPDATEPREVIEWDISPLAY Update pole positions in real-time during drag

    poleFamilies = getappdata(figHandle, 'poleFamilies');
    upperOnly = getappdata(figHandle, 'upperOnly');
    markerSize = getappdata(figHandle, 'markerSize');
    poleStyles = getappdata(figHandle, 'poleStyles');
    firstPoleFixed = getappdata(figHandle, 'firstPoleFixed');
    firstPoleCrystal = getappdata(figHandle, 'firstPoleCrystal');

    % Get or create preview markers
    hPreview = getappdata(figHandle, 'previewMarkers');

    % Delete old preview markers
    if ~isempty(hPreview)
        for i = 1:length(hPreview)
            if isvalid(hPreview(i))
                delete(hPreview(i));
            end
        end
    end

    % Create new preview markers
    hPreview = gobjects(0);
    hold(axStereo, 'on');

    for i = 1:length(poleFamilies)
        family = poleFamilies{i};
        poles = generatePoleFamily(family);

        styleField = ['p', family];
        if isfield(poleStyles, styleField)
            style = poleStyles.(styleField);
        else
            style = struct('color', rand(1,3), 'marker', 'o');
        end

        projectedPositions = [];

        for j = 1:size(poles, 1)
            pole = poles(j, :);
            pole = pole / norm(pole);
            worldPole = (R_preview * pole')';
            originalPole = pole;

            if upperOnly && worldPole(3) < 0
                worldPole = -worldPole;
                originalPole = -pole;
            end

            if worldPole(3) >= 0 || ~upperOnly
                [x, y] = stereoProject(worldPole);

                isDuplicate = false;
                for k = 1:size(projectedPositions, 1)
                    if norm([x, y] - projectedPositions(k, :)) < 0.01
                        isDuplicate = true;
                        break;
                    end
                end

                if x^2 + y^2 <= 1.01 && ~isDuplicate
                    % Check if locked pole
                    isLockedPole = false;
                    if firstPoleFixed && ~isempty(firstPoleCrystal)
                        if norm(abs(originalPole(:)) - abs(firstPoleCrystal(:))) < 0.01
                            isLockedPole = true;
                        end
                    end

                    % Create preview marker (slightly transparent)
                    if isLockedPole
                        h = scatter(axStereo, x, y, markerSize * 1.5, style.color, 'filled', ...
                                   'Marker', style.marker, 'MarkerEdgeColor', [1 0.5 0], ...
                                   'LineWidth', 2, 'MarkerFaceAlpha', 0.8, ...
                                   'HandleVisibility', 'off');
                    else
                        h = scatter(axStereo, x, y, markerSize, style.color, 'filled', ...
                                   'Marker', style.marker, 'MarkerEdgeColor', 'k', ...
                                   'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8, ...
                                   'HandleVisibility', 'off');
                    end

                    % Store for cleanup but don't add click callback during preview
                    hPreview(end+1) = h;
                    projectedPositions = [projectedPositions; x, y];
                end
            end
        end
    end

    hold(axStereo, 'off');
    setappdata(figHandle, 'previewMarkers', hPreview);

    % Store preview orientation for final application
    setappdata(figHandle, 'previewOrientation', R_preview);
end

function ori = getCurrentOrientation(figHandle)
% GETCURRENTORIENTATION Retrieve the current crystalOrientation from stereoProj figure
%
%   ORI = GETCURRENTORIENTATION(FIGHANDLE) returns the current crystalOrientation
%   object from an interactive stereoProj figure. Use this after interactive
%   manipulation to get the final orientation.
%
%   Input:
%       figHandle - Handle to the stereoProj figure
%
%   Output:
%       ori - crystalOrientation object with the current orientation
%
%   Example:
%       [~, ax] = stereoProj(eye(3));
%       % ... user interacts with the widget ...
%       ori = getCurrentOrientation(ancestor(ax, 'figure'));
%       disp(ori);
%
%   See also: stereoProj, crystalOrientation

    if ~ishandle(figHandle) || ~strcmp(get(figHandle, 'Type'), 'figure')
        error('Input must be a valid figure handle');
    end

    ori = getappdata(figHandle, 'crystalOrientationObject');

    if isempty(ori)
        error('No crystalOrientation found. Ensure the figure was created by stereoProj.');
    end
end
