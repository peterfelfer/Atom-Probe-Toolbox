classdef crystalOrientation
% CRYSTALORIENTATION Represents crystallographic orientation of a location/volume
%
%   The crystalOrientation class stores and manipulates crystal orientations
%   with support for multiple representation formats and all crystal systems.
%   Can also store full crystal structure data from CIF files.
%
%   Constructor Signatures:
%       ori = crystalOrientation(R, crystalSystem)
%       ori = crystalOrientation(R, crystalSystem, 'Name', Value, ...)
%       ori = crystalOrientation('euler', [phi1, Phi, phi2], crystalSystem, ...)
%       ori = crystalOrientation('quaternion', [w, x, y, z], crystalSystem, ...)
%
%   Static Factory Methods:
%       ori = crystalOrientation.fromCIF(filename)
%       ori = crystalOrientation.fromCIF(filename, 'Name', Value, ...)
%
%   Properties:
%       rotationMatrix  - 3x3 rotation matrix (crystal to world)
%       eulerAngles     - [phi1, Phi, phi2] in degrees (Bunge ZXZ)
%       quaternion      - [w, x, y, z] unit quaternion
%       position        - [x, y, z] location in world coordinates
%       crystalSystem   - 'cubic', 'hexagonal', 'tetragonal', etc.
%       latticeParams   - struct with a, b, c, alpha, beta, gamma
%       source          - origin of data ('EBSD', 'APT', 'manual', 'CIF')
%       spaceGroup      - Space group name (Hermann-Mauguin notation)
%       atomPositions   - Nx3 array of fractional coordinates
%       atomTypes       - Nx1 cell array of element symbols
%       atomLabels      - Nx1 cell array of atom site labels
%
%   Example:
%       % Create cubic orientation from rotation matrix
%       R = eye(3);
%       ori = crystalOrientation(R, 'cubic');
%
%       % Create from Euler angles with position
%       ori = crystalOrientation('euler', [45, 30, 60], 'cubic', ...
%                                'position', [10, 20, 30]);
%
%       % Create hexagonal with lattice parameters
%       ori = crystalOrientation(R, 'hexagonal', ...
%                                'latticeParams', struct('a', 3.2, 'c', 5.2));
%
%       % Load from CIF file
%       ori = crystalOrientation.fromCIF('structure.cif');
%
%   See also: stereoProj, plotOrientationCube

    properties (SetAccess = private)
        rotationMatrix  % 3x3 rotation matrix (crystal to world)
        eulerAngles     % [phi1, Phi, phi2] in degrees (Bunge ZXZ convention)
        quaternion      % [w, x, y, z] unit quaternion
        position        % [x, y, z] location in world coordinates
        crystalSystem   % Crystal system name
        latticeParams   % Lattice parameters struct
        source          % Data source identifier

        % Crystal structure properties (from CIF or manual input)
        spaceGroup      % Space group name (Hermann-Mauguin notation)
        atomPositions   % Nx3 array of fractional coordinates
        atomTypes       % Nx1 cell array of element symbols
        atomLabels      % Nx1 cell array of atom site labels
    end

    methods
        function obj = crystalOrientation(varargin)
            % CRYSTALORIENTATION Constructor
            %
            %   ori = crystalOrientation(R, crystalSystem)
            %   ori = crystalOrientation('euler', angles, crystalSystem)
            %   ori = crystalOrientation('quaternion', q, crystalSystem)

            if nargin == 0
                error('crystalOrientation requires at least orientation and crystal system inputs');
            end

            % Parse first argument to determine input type
            firstArg = varargin{1};

            if ischar(firstArg) || isstring(firstArg)
                % Input type specified: 'euler' or 'quaternion'
                inputType = lower(char(firstArg));
                if nargin < 3
                    error('Usage: crystalOrientation(''%s'', values, crystalSystem)', inputType);
                end
                orientationData = varargin{2};
                crystalSys = varargin{3};
                extraArgs = varargin(4:end);
            elseif isnumeric(firstArg) && all(size(firstArg) == [3, 3])
                % Rotation matrix input
                inputType = 'rotmat';
                orientationData = firstArg;
                if nargin < 2
                    error('Usage: crystalOrientation(R, crystalSystem)');
                end
                crystalSys = varargin{2};
                extraArgs = varargin(3:end);
            else
                error('First argument must be a 3x3 rotation matrix or ''euler''/''quaternion''');
            end

            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'position', [0, 0, 0], @(x) isnumeric(x) && numel(x)==3);
            addParameter(p, 'latticeParams', [], @(x) isempty(x) || isstruct(x));
            addParameter(p, 'source', 'manual', @(x) ischar(x) || isstring(x));
            addParameter(p, 'spaceGroup', '', @(x) ischar(x) || isstring(x));
            addParameter(p, 'atomPositions', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
            addParameter(p, 'atomTypes', {}, @iscell);
            addParameter(p, 'atomLabels', {}, @iscell);
            parse(p, extraArgs{:});

            obj.position = p.Results.position(:)';
            obj.source = char(p.Results.source);
            obj.spaceGroup = char(p.Results.spaceGroup);
            obj.atomPositions = p.Results.atomPositions;
            obj.atomTypes = p.Results.atomTypes;
            obj.atomLabels = p.Results.atomLabels;

            % Validate and set crystal system
            obj.crystalSystem = obj.validateCrystalSystem(crystalSys);

            % Set lattice parameters (defaults based on crystal system)
            obj.latticeParams = obj.setLatticeParams(p.Results.latticeParams);

            % Convert input to rotation matrix and compute all representations
            switch inputType
                case 'rotmat'
                    R = orientationData;
                    obj = obj.validateAndSetRotationMatrix(R);
                case 'euler'
                    if numel(orientationData) ~= 3
                        error('Euler angles must be [phi1, Phi, phi2]');
                    end
                    R = obj.euler2rotmat(orientationData(1), orientationData(2), orientationData(3));
                    obj = obj.validateAndSetRotationMatrix(R);
                case 'quaternion'
                    if numel(orientationData) ~= 4
                        error('Quaternion must be [w, x, y, z]');
                    end
                    R = obj.quat2rotmat(orientationData);
                    obj = obj.validateAndSetRotationMatrix(R);
                otherwise
                    error('Unknown input type: %s', inputType);
            end
        end

        function R = getRotationMatrix(obj)
            % GETROTATIONMATRIX Return the 3x3 rotation matrix
            R = obj.rotationMatrix;
        end

        function [phi1, Phi, phi2] = getEulerAngles(obj)
            % GETEULERANGLES Return Euler angles in degrees (Bunge ZXZ)
            phi1 = obj.eulerAngles(1);
            Phi = obj.eulerAngles(2);
            phi2 = obj.eulerAngles(3);
        end

        function q = getQuaternion(obj)
            % GETQUATERNION Return unit quaternion [w, x, y, z]
            q = obj.quaternion;
        end

        function obj = setPosition(obj, pos)
            % SETPOSITION Set the position
            if ~isnumeric(pos) || numel(pos) ~= 3
                error('Position must be [x, y, z]');
            end
            obj.position = pos(:)';
        end

        function symOps = getSymmetryOperators(obj)
            % GETSYMMETRYOPERATORS Return symmetry operators for the crystal system
            %   Returns cell array of 3x3 rotation matrices
            symOps = obj.computeSymmetryOperators(obj.crystalSystem);
        end

        function equivOri = getEquivalentOrientations(obj)
            % GETEQUIVALENTORIENTATIONS Return all symmetrically equivalent orientations
            %   Returns array of crystalOrientation objects
            symOps = obj.getSymmetryOperators();
            nOps = length(symOps);
            equivOri = repmat(obj, nOps, 1);

            for i = 1:nOps
                R_equiv = obj.rotationMatrix * symOps{i};
                equivOri(i) = crystalOrientation(R_equiv, obj.crystalSystem, ...
                    'position', obj.position, ...
                    'latticeParams', obj.latticeParams, ...
                    'source', obj.source);
            end
        end

        function angle = misorientationAngle(obj, other)
            % MISORIENTATIONANGLE Calculate minimum misorientation angle (degrees)
            %   angle = misorientationAngle(obj, other)
            %   Accounts for crystal symmetry

            if ~isa(other, 'crystalOrientation')
                error('Input must be a crystalOrientation object');
            end

            % Get symmetry operators
            symOps = obj.getSymmetryOperators();

            % Misorientation: R_mis = R2^T * R1
            R1 = obj.rotationMatrix;
            R2 = other.rotationMatrix;

            minAngle = 180;  % Start with maximum possible

            % Check all symmetrically equivalent misorientations
            for i = 1:length(symOps)
                R_sym = symOps{i};
                R_mis = R2' * R1 * R_sym;

                % Calculate rotation angle from trace
                traceR = trace(R_mis);
                traceR = max(-1, min(3, traceR));  % Clamp for numerical stability
                angleRad = acos((traceR - 1) / 2);
                angleDeg = rad2deg(angleRad);

                if angleDeg < minAngle
                    minAngle = angleDeg;
                end
            end

            angle = minAngle;
        end

        function [axis, angle] = misorientationAxisAngle(obj, other)
            % MISORIENTATIONAXISANGLE Calculate misorientation axis and angle
            %   [axis, angle] = misorientationAxisAngle(obj, other)

            if ~isa(other, 'crystalOrientation')
                error('Input must be a crystalOrientation object');
            end

            % Get symmetry operators
            symOps = obj.getSymmetryOperators();

            R1 = obj.rotationMatrix;
            R2 = other.rotationMatrix;

            minAngle = 180;
            minAxis = [0, 0, 1];

            % Check all symmetrically equivalent misorientations
            for i = 1:length(symOps)
                R_sym = symOps{i};
                R_mis = R2' * R1 * R_sym;

                % Calculate rotation angle from trace
                traceR = trace(R_mis);
                traceR = max(-1, min(3, traceR));
                angleRad = acos((traceR - 1) / 2);
                angleDeg = rad2deg(angleRad);

                if angleDeg < minAngle
                    minAngle = angleDeg;

                    % Extract axis from rotation matrix
                    if abs(angleRad) < 1e-10
                        minAxis = [0, 0, 1];  % Arbitrary axis for zero rotation
                    else
                        minAxis = [R_mis(3,2) - R_mis(2,3);
                                   R_mis(1,3) - R_mis(3,1);
                                   R_mis(2,1) - R_mis(1,2)] / (2 * sin(angleRad));
                        minAxis = minAxis' / norm(minAxis);
                    end
                end
            end

            axis = minAxis;
            angle = minAngle;
        end

        function h = plotStereographic(obj, varargin)
            % PLOTSTEREOGRAPHIC Plot stereographic projection
            %   h = plotStereographic(obj)
            %   h = plotStereographic(obj, 'Name', Value, ...)
            %
            %   Passes all arguments to stereoProj function

            h = stereoProj(obj.rotationMatrix, varargin{:});
        end

        function h = plotCube(obj, varargin)
            % PLOTCUBE Plot orientation cube at position
            %   h = plotCube(obj)
            %   h = plotCube(obj, 'Name', Value, ...)
            %
            %   Passes all arguments to plotOrientationCube function

            h = plotOrientationCube(obj.rotationMatrix, ...
                'position', obj.position, varargin{:});
        end

        function disp(obj)
            % DISP Display orientation information
            fprintf('crystalOrientation:\n');
            fprintf('  Crystal System: %s\n', obj.crystalSystem);
            if ~isempty(obj.spaceGroup)
                fprintf('  Space Group: %s\n', obj.spaceGroup);
            end
            fprintf('  Position: [%.3f, %.3f, %.3f]\n', obj.position);
            fprintf('  Source: %s\n', obj.source);
            fprintf('  Euler Angles (Bunge ZXZ): [%.2f°, %.2f°, %.2f°]\n', ...
                obj.eulerAngles(1), obj.eulerAngles(2), obj.eulerAngles(3));
            fprintf('  Quaternion: [%.4f, %.4f, %.4f, %.4f]\n', obj.quaternion);
            fprintf('  Rotation Matrix:\n');
            fprintf('    [%7.4f %7.4f %7.4f]\n', obj.rotationMatrix(1,:));
            fprintf('    [%7.4f %7.4f %7.4f]\n', obj.rotationMatrix(2,:));
            fprintf('    [%7.4f %7.4f %7.4f]\n', obj.rotationMatrix(3,:));
            if ~isempty(obj.latticeParams)
                fprintf('  Lattice Parameters:\n');
                fprintf('    a=%.4f, b=%.4f, c=%.4f\n', ...
                    obj.latticeParams.a, obj.latticeParams.b, obj.latticeParams.c);
                fprintf('    alpha=%.2f°, beta=%.2f°, gamma=%.2f°\n', ...
                    obj.latticeParams.alpha, obj.latticeParams.beta, obj.latticeParams.gamma);
            end
            if ~isempty(obj.atomPositions)
                nAtoms = size(obj.atomPositions, 1);
                fprintf('  Atom Sites: %d\n', nAtoms);
                maxShow = min(5, nAtoms);
                for i = 1:maxShow
                    if ~isempty(obj.atomLabels) && length(obj.atomLabels) >= i
                        label = obj.atomLabels{i};
                    else
                        label = sprintf('Site%d', i);
                    end
                    if ~isempty(obj.atomTypes) && length(obj.atomTypes) >= i
                        elem = obj.atomTypes{i};
                    else
                        elem = '?';
                    end
                    fprintf('    %s (%s): [%.4f, %.4f, %.4f]\n', ...
                        label, elem, obj.atomPositions(i,:));
                end
                if nAtoms > maxShow
                    fprintf('    ... and %d more atoms\n', nAtoms - maxShow);
                end
            end
        end

        function atoms = getAtomPositions(obj, varargin)
            % GETATOMPOSITIONS Get atom positions in world or fractional coordinates
            %   atoms = getAtomPositions(obj) returns fractional coordinates
            %   atoms = getAtomPositions(obj, 'cartesian', true) returns Cartesian
            %           coordinates in the crystal frame
            %   atoms = getAtomPositions(obj, 'world', true) returns Cartesian
            %           coordinates in the world frame (applying rotation)

            p = inputParser;
            addParameter(p, 'cartesian', false, @islogical);
            addParameter(p, 'world', false, @islogical);
            parse(p, varargin{:});

            if isempty(obj.atomPositions)
                atoms = [];
                return;
            end

            if ~p.Results.cartesian && ~p.Results.world
                % Return fractional coordinates
                atoms = obj.atomPositions;
            else
                % Convert to Cartesian (crystal frame)
                atoms = obj.fractionalToCartesian(obj.atomPositions);

                if p.Results.world
                    % Transform to world coordinates
                    atoms = (obj.rotationMatrix * atoms')';
                end
            end
        end

        function cartesian = fractionalToCartesian(obj, fractional)
            % FRACTIONALTOCARTESIAN Convert fractional to Cartesian coordinates
            %   Uses lattice parameters to build transformation matrix

            a = obj.latticeParams.a;
            b = obj.latticeParams.b;
            c = obj.latticeParams.c;
            alpha = deg2rad(obj.latticeParams.alpha);
            beta = deg2rad(obj.latticeParams.beta);
            gamma = deg2rad(obj.latticeParams.gamma);

            % Build transformation matrix (fractional to Cartesian)
            % Convention: a along x, b in xy plane
            cosAlpha = cos(alpha);
            cosBeta = cos(beta);
            cosGamma = cos(gamma);
            sinGamma = sin(gamma);

            omega = a * b * c * sqrt(1 - cosAlpha^2 - cosBeta^2 - cosGamma^2 + 2*cosAlpha*cosBeta*cosGamma);

            M = [a,  b*cosGamma,  c*cosBeta;
                 0,  b*sinGamma,  c*(cosAlpha - cosBeta*cosGamma)/sinGamma;
                 0,  0,           omega/(a*b*sinGamma)];

            cartesian = (M * fractional')';
        end
    end

    methods (Static)
        function obj = fromCIF(filename, varargin)
            % FROMCIF Create crystalOrientation from a CIF file
            %
            %   ORI = crystalOrientation.fromCIF(FILENAME) reads crystallographic
            %   data from a CIF file and creates a crystalOrientation object with
            %   identity orientation.
            %
            %   ORI = crystalOrientation.fromCIF(FILENAME, 'Name', Value, ...)
            %   allows specifying additional options.
            %
            %   Optional Name-Value Pairs:
            %       'orientation' - Initial orientation as 3x3 rotation matrix
            %                       or 'euler' followed by [phi1, Phi, phi2]
            %                       (default: eye(3))
            %       'position'    - [x, y, z] position (default: [0,0,0])
            %
            %   Example:
            %       ori = crystalOrientation.fromCIF('NaCl.cif');
            %       ori = crystalOrientation.fromCIF('Fe.cif', 'position', [10,20,30]);

            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
            addParameter(p, 'orientation', eye(3), @(x) isnumeric(x) && all(size(x)==[3,3]));
            addParameter(p, 'position', [0, 0, 0], @(x) isnumeric(x) && numel(x)==3);
            parse(p, filename, varargin{:});

            filename = char(p.Results.filename);
            R = p.Results.orientation;
            position = p.Results.position;

            % Parse the CIF file
            cifData = parseCIFFile(filename);

            % Create the crystalOrientation object
            obj = crystalOrientation(R, cifData.crystalSystem, ...
                'position', position, ...
                'source', 'CIF', ...
                'latticeParams', cifData.latticeParams, ...
                'spaceGroup', cifData.spaceGroup, ...
                'atomPositions', cifData.atomPositions, ...
                'atomTypes', cifData.atomTypes, ...
                'atomLabels', cifData.atomLabels);
        end
    end

    methods (Access = private)
        function obj = validateAndSetRotationMatrix(obj, R)
            % Validate rotation matrix and compute all representations

            % Check orthogonality
            RtR = R' * R;
            if max(abs(RtR - eye(3)), [], 'all') > 1e-6
                warning('Rotation matrix not perfectly orthogonal, orthogonalizing...');
                [U, ~, V] = svd(R);
                R = U * V';
            end

            % Check determinant
            if det(R) < 0
                error('Rotation matrix has negative determinant (reflection)');
            end

            obj.rotationMatrix = R;
            obj.eulerAngles = obj.rotmat2euler(R);
            obj.quaternion = obj.rotmat2quat(R);
        end

        function crystalSys = validateCrystalSystem(~, input)
            % Validate and normalize crystal system name
            validSystems = {'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', ...
                           'monoclinic', 'triclinic'};
            input = lower(char(input));

            if ~ismember(input, validSystems)
                error('Invalid crystal system: %s. Valid options: %s', ...
                    input, strjoin(validSystems, ', '));
            end

            crystalSys = input;
        end

        function params = setLatticeParams(obj, inputParams)
            % Set lattice parameters with defaults based on crystal system

            % Default parameters by crystal system
            switch obj.crystalSystem
                case 'cubic'
                    defaults = struct('a', 1, 'b', 1, 'c', 1, ...
                                     'alpha', 90, 'beta', 90, 'gamma', 90);
                case 'hexagonal'
                    defaults = struct('a', 1, 'b', 1, 'c', 1.633, ...
                                     'alpha', 90, 'beta', 90, 'gamma', 120);
                case 'tetragonal'
                    defaults = struct('a', 1, 'b', 1, 'c', 1.5, ...
                                     'alpha', 90, 'beta', 90, 'gamma', 90);
                case 'orthorhombic'
                    defaults = struct('a', 1, 'b', 1.2, 'c', 1.5, ...
                                     'alpha', 90, 'beta', 90, 'gamma', 90);
                case 'monoclinic'
                    defaults = struct('a', 1, 'b', 1.2, 'c', 1.5, ...
                                     'alpha', 90, 'beta', 100, 'gamma', 90);
                case 'triclinic'
                    defaults = struct('a', 1, 'b', 1.2, 'c', 1.5, ...
                                     'alpha', 80, 'beta', 85, 'gamma', 70);
            end

            if isempty(inputParams)
                params = defaults;
            else
                % Merge input with defaults
                params = defaults;
                fields = fieldnames(inputParams);
                for i = 1:length(fields)
                    params.(fields{i}) = inputParams.(fields{i});
                end
            end
        end

        function symOps = computeSymmetryOperators(~, crystalSys)
            % Compute symmetry operators for the given crystal system
            % Returns cell array of 3x3 rotation matrices

            switch crystalSys
                case 'cubic'
                    % 24 proper rotations of cubic symmetry (m-3m point group)
                    symOps = cell(24, 1);
                    idx = 1;

                    % Identity
                    symOps{idx} = eye(3); idx = idx + 1;

                    % 90° rotations about <100> (6 total)
                    axes = {[1,0,0], [0,1,0], [0,0,1]};
                    for a = 1:3
                        ax = axes{a};
                        for angle = [90, 180, 270]
                            symOps{idx} = rotationMatrix3D(ax, angle);
                            idx = idx + 1;
                        end
                    end

                    % 120° and 240° rotations about <111> (8 total)
                    diag_axes = {[1,1,1], [1,1,-1], [1,-1,1], [-1,1,1]};
                    for a = 1:4
                        ax = diag_axes{a} / norm(diag_axes{a});
                        for angle = [120, 240]
                            symOps{idx} = rotationMatrix3D(ax, angle);
                            idx = idx + 1;
                        end
                    end

                    % 180° rotations about <110> (6 total)
                    edge_axes = {[1,1,0], [1,-1,0], [1,0,1], [1,0,-1], [0,1,1], [0,1,-1]};
                    for a = 1:6
                        ax = edge_axes{a} / norm(edge_axes{a});
                        symOps{idx} = rotationMatrix3D(ax, 180);
                        idx = idx + 1;
                    end

                case 'hexagonal'
                    % 12 proper rotations of hexagonal symmetry (6/mmm)
                    symOps = cell(12, 1);

                    % Rotations about c-axis (6-fold)
                    for i = 0:5
                        angle = i * 60;
                        symOps{i+1} = rotationMatrix3D([0,0,1], angle);
                    end

                    % 180° rotations about a-axes
                    for i = 0:5
                        angle_rad = i * pi / 3;
                        ax = [cos(angle_rad), sin(angle_rad), 0];
                        symOps{7+i} = rotationMatrix3D(ax, 180);
                    end

                case 'tetragonal'
                    % 8 proper rotations (4/mmm)
                    symOps = cell(8, 1);

                    % Rotations about c-axis (4-fold)
                    for i = 0:3
                        angle = i * 90;
                        symOps{i+1} = rotationMatrix3D([0,0,1], angle);
                    end

                    % 180° rotations about <100> and <110>
                    symOps{5} = rotationMatrix3D([1,0,0], 180);
                    symOps{6} = rotationMatrix3D([0,1,0], 180);
                    symOps{7} = rotationMatrix3D([1,1,0]/sqrt(2), 180);
                    symOps{8} = rotationMatrix3D([1,-1,0]/sqrt(2), 180);

                case 'orthorhombic'
                    % 4 proper rotations (mmm)
                    symOps = cell(4, 1);
                    symOps{1} = eye(3);
                    symOps{2} = rotationMatrix3D([1,0,0], 180);
                    symOps{3} = rotationMatrix3D([0,1,0], 180);
                    symOps{4} = rotationMatrix3D([0,0,1], 180);

                case 'monoclinic'
                    % 2 proper rotations (2/m)
                    symOps = cell(2, 1);
                    symOps{1} = eye(3);
                    symOps{2} = rotationMatrix3D([0,1,0], 180);  % 2-fold about b

                case 'triclinic'
                    % Only identity
                    symOps = {eye(3)};

                otherwise
                    error('Unknown crystal system: %s', crystalSys);
            end
        end
    end

    methods (Static, Access = private)
        function euler = rotmat2euler(R)
            % Convert rotation matrix to Euler angles (Bunge ZXZ convention)
            % Returns [phi1, Phi, phi2] in degrees

            R = max(-1, min(1, R));  % Clamp for numerical stability

            Phi = acosd(R(3,3));

            if abs(R(3,3)) > 0.9999
                % Gimbal lock case
                phi1 = atan2d(-R(1,2), R(1,1));
                phi2 = 0;
                if R(3,3) < 0
                    Phi = 180;
                else
                    Phi = 0;
                end
            else
                phi1 = atan2d(R(3,1), -R(3,2));
                phi2 = atan2d(R(1,3), R(2,3));
            end

            % Normalize to [0, 360)
            phi1 = mod(phi1, 360);
            phi2 = mod(phi2, 360);

            euler = [phi1, Phi, phi2];
        end

        function R = euler2rotmat(phi1, Phi, phi2)
            % Convert Euler angles to rotation matrix (Bunge ZXZ convention)
            % Input angles in degrees

            phi1 = deg2rad(phi1);
            Phi = deg2rad(Phi);
            phi2 = deg2rad(phi2);

            c1 = cos(phi1); s1 = sin(phi1);
            c = cos(Phi);   s = sin(Phi);
            c2 = cos(phi2); s2 = sin(phi2);

            R = [c1*c2 - s1*c*s2,  -c1*s2 - s1*c*c2,  s1*s;
                 s1*c2 + c1*c*s2,  -s1*s2 + c1*c*c2, -c1*s;
                 s*s2,              s*c2,             c];
        end

        function q = rotmat2quat(R)
            % Convert rotation matrix to quaternion [w, x, y, z]

            tr = trace(R);

            if tr > 0
                s = sqrt(tr + 1) * 2;
                w = s / 4;
                x = (R(3,2) - R(2,3)) / s;
                y = (R(1,3) - R(3,1)) / s;
                z = (R(2,1) - R(1,2)) / s;
            elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
                s = sqrt(1 + R(1,1) - R(2,2) - R(3,3)) * 2;
                w = (R(3,2) - R(2,3)) / s;
                x = s / 4;
                y = (R(1,2) + R(2,1)) / s;
                z = (R(1,3) + R(3,1)) / s;
            elseif R(2,2) > R(3,3)
                s = sqrt(1 + R(2,2) - R(1,1) - R(3,3)) * 2;
                w = (R(1,3) - R(3,1)) / s;
                x = (R(1,2) + R(2,1)) / s;
                y = s / 4;
                z = (R(2,3) + R(3,2)) / s;
            else
                s = sqrt(1 + R(3,3) - R(1,1) - R(2,2)) * 2;
                w = (R(2,1) - R(1,2)) / s;
                x = (R(1,3) + R(3,1)) / s;
                y = (R(2,3) + R(3,2)) / s;
                z = s / 4;
            end

            q = [w, x, y, z];

            % Normalize
            q = q / norm(q);

            % Ensure w >= 0 for unique representation
            if q(1) < 0
                q = -q;
            end
        end

        function R = quat2rotmat(q)
            % Convert quaternion [w, x, y, z] to rotation matrix

            q = q / norm(q);  % Ensure normalized
            w = q(1); x = q(2); y = q(3); z = q(4);

            R = [1 - 2*y^2 - 2*z^2,     2*x*y - 2*z*w,       2*x*z + 2*y*w;
                 2*x*y + 2*z*w,         1 - 2*x^2 - 2*z^2,   2*y*z - 2*x*w;
                 2*x*z - 2*y*w,         2*y*z + 2*x*w,       1 - 2*x^2 - 2*y^2];
        end
    end
end

function R = rotationMatrix3D(axis, angle_deg)
    % Helper function to create rotation matrix from axis-angle
    % axis: unit vector [x, y, z]
    % angle_deg: rotation angle in degrees

    axis = axis(:)' / norm(axis);
    angle = deg2rad(angle_deg);

    c = cos(angle);
    s = sin(angle);
    t = 1 - c;

    x = axis(1);
    y = axis(2);
    z = axis(3);

    R = [t*x*x + c,    t*x*y - s*z,  t*x*z + s*y;
         t*x*y + s*z,  t*y*y + c,    t*y*z - s*x;
         t*x*z - s*y,  t*y*z + s*x,  t*z*z + c];
end

function cifData = parseCIFFile(filename)
% PARSECIFFILE Parse a CIF (Crystallographic Information File)
%   Returns a structure with crystal structure data
%
%   cifData = parseCIFFile(filename)
%
%   Output fields:
%       crystalSystem  - Detected crystal system
%       spaceGroup     - Space group name (Hermann-Mauguin)
%       latticeParams  - struct with a, b, c, alpha, beta, gamma
%       atomPositions  - Nx3 fractional coordinates
%       atomTypes      - Nx1 cell array of element symbols
%       atomLabels     - Nx1 cell array of atom site labels

    if ~exist(filename, 'file')
        error('CIF file not found: %s', filename);
    end

    % Read the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open CIF file: %s', filename);
    end
    content = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = content{1};

    % Initialize output
    cifData = struct();
    cifData.spaceGroup = '';
    cifData.latticeParams = struct('a', 1, 'b', 1, 'c', 1, 'alpha', 90, 'beta', 90, 'gamma', 90);
    cifData.atomPositions = [];
    cifData.atomTypes = {};
    cifData.atomLabels = {};

    % Parse lattice parameters and space group
    for i = 1:length(lines)
        line = strtrim(lines{i});

        % Skip empty lines and comments
        if isempty(line) || line(1) == '#'
            continue;
        end

        % Lattice parameters
        if startsWith(line, '_cell_length_a')
            cifData.latticeParams.a = parseNumericValue(line);
        elseif startsWith(line, '_cell_length_b')
            cifData.latticeParams.b = parseNumericValue(line);
        elseif startsWith(line, '_cell_length_c')
            cifData.latticeParams.c = parseNumericValue(line);
        elseif startsWith(line, '_cell_angle_alpha')
            cifData.latticeParams.alpha = parseNumericValue(line);
        elseif startsWith(line, '_cell_angle_beta')
            cifData.latticeParams.beta = parseNumericValue(line);
        elseif startsWith(line, '_cell_angle_gamma')
            cifData.latticeParams.gamma = parseNumericValue(line);

        % Space group (various CIF tags)
        elseif startsWith(line, '_symmetry_space_group_name_H-M') || ...
               startsWith(line, '_space_group_name_H-M_alt') || ...
               startsWith(line, '_space_group.name_H-M_alt')
            cifData.spaceGroup = parseStringValue(line);

        elseif startsWith(line, '_symmetry_Int_Tables_number') || ...
               startsWith(line, '_space_group_IT_number')
            % Store space group number if name not found
            if isempty(cifData.spaceGroup)
                cifData.spaceGroup = sprintf('SG#%d', round(parseNumericValue(line)));
            end
        end
    end

    % Parse atom sites (look for loop_ containing atom site data)
    [cifData.atomPositions, cifData.atomTypes, cifData.atomLabels] = parseAtomSites(lines);

    % Determine crystal system from space group or lattice parameters
    cifData.crystalSystem = determineCrystalSystem(cifData.spaceGroup, cifData.latticeParams);
end

function val = parseNumericValue(line)
% Parse numeric value from a CIF line, handling uncertainty notation
    tokens = strsplit(strtrim(line));
    if length(tokens) >= 2
        valStr = tokens{2};
        % Remove uncertainty in parentheses, e.g., "5.431(1)" -> "5.431"
        valStr = regexprep(valStr, '\([0-9]+\)', '');
        val = str2double(valStr);
        if isnan(val)
            val = 0;
        end
    else
        val = 0;
    end
end

function str = parseStringValue(line)
% Parse string value from a CIF line
    tokens = strsplit(strtrim(line));
    if length(tokens) >= 2
        % Rejoin tokens after the tag, removing quotes
        str = strjoin(tokens(2:end), ' ');
        str = strrep(str, '''', '');
        str = strrep(str, '"', '');
        str = strtrim(str);
    else
        str = '';
    end
end

function [positions, types, labels] = parseAtomSites(lines)
% Parse atom site loop from CIF
    positions = [];
    types = {};
    labels = {};

    % Find loop_ containing atom sites
    inLoop = false;
    loopColumns = {};
    xCol = 0; yCol = 0; zCol = 0;
    typeCol = 0; labelCol = 0;

    i = 1;
    while i <= length(lines)
        line = strtrim(lines{i});

        if strcmp(line, 'loop_')
            inLoop = true;
            loopColumns = {};
            xCol = 0; yCol = 0; zCol = 0; typeCol = 0; labelCol = 0;
            i = i + 1;
            continue;
        end

        if inLoop
            % Check if line contains atom_site column header (with or without leading space)
            if contains(line, '_atom_site')
                % Column header - extract the tag name
                loopColumns{end+1} = strtrim(line);

                % Identify which column is which
                colIdx = length(loopColumns);
                if contains(line, '_atom_site_fract_x')
                    xCol = colIdx;
                elseif contains(line, '_atom_site_fract_y')
                    yCol = colIdx;
                elseif contains(line, '_atom_site_fract_z')
                    zCol = colIdx;
                elseif contains(line, '_atom_site_type_symbol')
                    typeCol = colIdx;
                elseif contains(line, '_atom_site_label')
                    labelCol = colIdx;
                end
            elseif ~isempty(line) && ~startsWith(line, '_') && line(1) ~= '#' && ~strcmp(line, 'loop_')
                % Data line - only process if we have the required columns
                if xCol > 0 && yCol > 0 && zCol > 0
                    % Split by whitespace, handling multiple spaces
                    tokens = strsplit(strtrim(line));
                    tokens = tokens(~cellfun('isempty', tokens));  % Remove empty tokens

                    maxCol = max([xCol, yCol, zCol, typeCol, labelCol]);
                    if length(tokens) >= maxCol || length(tokens) >= max([xCol, yCol, zCol])
                        x = parseCoordinate(tokens{xCol});
                        y = parseCoordinate(tokens{yCol});
                        z = parseCoordinate(tokens{zCol});
                        positions = [positions; x, y, z];

                        if typeCol > 0 && length(tokens) >= typeCol
                            types{end+1} = tokens{typeCol};
                        else
                            types{end+1} = '';
                        end

                        if labelCol > 0 && length(tokens) >= labelCol
                            labels{end+1} = tokens{labelCol};
                        else
                            labels{end+1} = '';
                        end
                    end
                end
            elseif strcmp(line, 'loop_') || (startsWith(line, '_') && ~contains(line, '_atom_site'))
                % End of this loop, check if we found atom sites
                if ~isempty(positions)
                    break;
                end
                % Reset for next loop
                inLoop = strcmp(line, 'loop_');
                if inLoop
                    loopColumns = {};
                    xCol = 0; yCol = 0; zCol = 0; typeCol = 0; labelCol = 0;
                end
            end
        end
        i = i + 1;
    end

    types = types(:);
    labels = labels(:);
end

function val = parseCoordinate(str)
% Parse coordinate value, removing uncertainty
    str = regexprep(str, '\([0-9]+\)', '');
    val = str2double(str);
    if isnan(val)
        val = 0;
    end
end

function crystalSystem = determineCrystalSystem(spaceGroup, latticeParams)
% Determine crystal system from space group or lattice parameters

    % Try to determine from space group name
    if ~isempty(spaceGroup)
        sg = upper(spaceGroup);

        % Cubic: contains 'M3' or 'm-3' or specific cubic groups
        if contains(sg, 'M3') || contains(sg, 'M-3') || ...
           contains(sg, '23') || contains(sg, '432') || contains(sg, '-43')
            crystalSystem = 'cubic';
            return;
        end

        % Hexagonal: starts with P6, R3, etc.
        if startsWith(sg, 'P6') || startsWith(sg, 'P-6') || contains(sg, '/M')
            if latticeParams.gamma == 120 || abs(latticeParams.gamma - 120) < 0.1
                crystalSystem = 'hexagonal';
                return;
            end
        end

        % Trigonal/Rhombohedral
        if startsWith(sg, 'R3') || startsWith(sg, 'R-3') || startsWith(sg, 'P3')
            crystalSystem = 'hexagonal';  % Often grouped with hexagonal
            return;
        end

        % Tetragonal: P4, I4
        if startsWith(sg, 'P4') || startsWith(sg, 'I4') || startsWith(sg, 'P-4') || startsWith(sg, 'I-4')
            crystalSystem = 'tetragonal';
            return;
        end

        % Orthorhombic: Pnma, Cmcm, etc.
        if startsWith(sg, 'P') || startsWith(sg, 'C') || startsWith(sg, 'I') || startsWith(sg, 'F')
            if contains(sg, 'MM') || contains(sg, '222') || contains(sg, 'MC') || contains(sg, 'MA')
                crystalSystem = 'orthorhombic';
                return;
            end
        end

        % Monoclinic: P2, C2, etc.
        if startsWith(sg, 'P2') || startsWith(sg, 'C2') || contains(sg, '/M') || contains(sg, '2/M')
            if ~contains(sg, 'MM')
                crystalSystem = 'monoclinic';
                return;
            end
        end
    end

    % Fall back to lattice parameter analysis
    a = latticeParams.a;
    b = latticeParams.b;
    c = latticeParams.c;
    alpha = latticeParams.alpha;
    beta = latticeParams.beta;
    gamma = latticeParams.gamma;

    tol = 0.01;  % Tolerance for comparisons

    % Cubic: a = b = c, all angles 90
    if abs(a-b)/a < tol && abs(b-c)/b < tol && ...
       abs(alpha-90) < 0.1 && abs(beta-90) < 0.1 && abs(gamma-90) < 0.1
        crystalSystem = 'cubic';
        return;
    end

    % Hexagonal: a = b != c, alpha = beta = 90, gamma = 120
    if abs(a-b)/a < tol && abs(alpha-90) < 0.1 && abs(beta-90) < 0.1 && abs(gamma-120) < 0.1
        crystalSystem = 'hexagonal';
        return;
    end

    % Tetragonal: a = b != c, all angles 90
    if abs(a-b)/a < tol && abs(alpha-90) < 0.1 && abs(beta-90) < 0.1 && abs(gamma-90) < 0.1
        crystalSystem = 'tetragonal';
        return;
    end

    % Orthorhombic: a != b != c, all angles 90
    if abs(alpha-90) < 0.1 && abs(beta-90) < 0.1 && abs(gamma-90) < 0.1
        crystalSystem = 'orthorhombic';
        return;
    end

    % Monoclinic: alpha = gamma = 90, beta != 90
    if abs(alpha-90) < 0.1 && abs(gamma-90) < 0.1 && abs(beta-90) >= 0.1
        crystalSystem = 'monoclinic';
        return;
    end

    % Default to triclinic
    crystalSystem = 'triclinic';
end
