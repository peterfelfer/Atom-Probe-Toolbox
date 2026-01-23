function [angle, axis, details] = misorientation(orientation1, orientation2, crystalSystem)
% MISORIENTATION Calculate misorientation between two crystal orientations
%
% [angle, axis, details] = misorientation(orientation1, orientation2)
% [angle, axis, details] = misorientation(orientation1, orientation2, 'cubic')
%
% Calculates the misorientation angle and axis between two crystal
% orientations, accounting for crystal symmetry to find the minimum
% (disorientation) angle.
%
% INPUT:
%   orientation1 - First orientation as:
%                  - 3x3 rotation matrix
%                  - 4-element quaternion [w, x, y, z]
%                  - 3-element Euler angles [phi1, Phi, phi2] in degrees (Bunge convention)
%                  - crystalOrientation object
%   orientation2 - Second orientation (same format as orientation1)
%   crystalSystem - Crystal symmetry (default: 'cubic')
%                   Options: 'cubic', 'hexagonal', 'tetragonal', 'orthorhombic',
%                           'monoclinic', 'triclinic'
%
% OUTPUT:
%   angle   - Misorientation angle in degrees (disorientation - minimum angle)
%   axis    - Rotation axis [x, y, z] (unit vector in crystal frame)
%   details - Structure with additional information:
%       .fullAngle      - Angle before symmetry reduction
%       .quaternion     - Misorientation quaternion
%       .rotationMatrix - Misorientation rotation matrix
%       .eulerAngles    - Euler angles [phi1, Phi, phi2]
%       .symmetryOp     - Symmetry operator used
%       .cslType        - CSL boundary type if applicable (Sigma value)
%
% THEORY:
%   Misorientation: g_mis = g2 * g1^(-1)
%   Disorientation: minimum angle considering all crystal symmetry equivalents
%
%   For cubic crystals, the maximum disorientation angle is 62.8 degrees.
%
% EXAMPLES:
%   % From rotation matrices
%   R1 = eye(3);
%   R2 = rotationMatrix([0 0 1], 45);  % 45 deg rotation about z
%   [angle, axis] = misorientation(R1, R2);
%
%   % From Euler angles (Bunge convention, degrees)
%   euler1 = [0, 0, 0];
%   euler2 = [45, 30, 60];
%   [angle, axis] = misorientation(euler1, euler2, 'cubic');
%
%   % From quaternions
%   q1 = [1, 0, 0, 0];  % Identity
%   q2 = [0.9659, 0, 0, 0.2588];  % 30 deg about z
%   [angle, axis] = misorientation(q1, q2);
%
% REFERENCES:
%   Randle & Engler (2000) "Texture Analysis"
%   Morawiec (2004) "Orientations and Rotations"
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    orientation1
    orientation2
    crystalSystem (1,:) char {mustBeMember(crystalSystem, ...
        {'cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'triclinic'})} = 'cubic'
end

% Convert inputs to rotation matrices
R1 = toRotationMatrix(orientation1);
R2 = toRotationMatrix(orientation2);

% Calculate misorientation
R_mis = R2 * R1';

% Get crystal symmetry operators
symOps = getSymmetryOperators(crystalSystem);
nSym = size(symOps, 3);

% Find minimum angle (disorientation) considering symmetry
minAngle = inf;
minAxis = [0; 0; 1];
minSymOp = 1;
minQuaternion = [1; 0; 0; 0];

for i = 1:nSym
    for j = 1:nSym
        % Apply symmetry: S_i * R_mis * S_j
        R_test = symOps(:,:,i) * R_mis * symOps(:,:,j);

        % Convert to angle-axis
        [testAngle, testAxis] = rotationMatrixToAngleAxis(R_test);

        % Check if this is the minimum angle
        if testAngle < minAngle
            minAngle = testAngle;
            minAxis = testAxis;
            minSymOp = i * 100 + j;  % Encode both operators
            minQuaternion = rotationMatrixToQuaternion(R_test);
        end
    end
end

angle = minAngle;
axis = minAxis;

% Build details structure
details = struct();
details.fullAngle = rad2deg(acos((trace(R_mis) - 1) / 2));
details.quaternion = minQuaternion;
details.rotationMatrix = quaternionToRotationMatrix(minQuaternion);
details.eulerAngles = rotationMatrixToEuler(details.rotationMatrix);
details.symmetryOp = minSymOp;
details.crystalSystem = crystalSystem;

% Check for CSL boundaries
details.cslType = checkCSL(angle, axis, crystalSystem);

end

%% Input Conversion Functions
function R = toRotationMatrix(orientation)
    % Convert various orientation representations to rotation matrix

    if isa(orientation, 'crystalOrientation')
        % crystalOrientation object
        R = orientation.rotationMatrix;

    elseif isnumeric(orientation)
        if isequal(size(orientation), [3, 3])
            % Already a rotation matrix
            R = orientation;

        elseif numel(orientation) == 4
            % Quaternion [w, x, y, z]
            R = quaternionToRotationMatrix(orientation(:));

        elseif numel(orientation) == 3
            % Euler angles [phi1, Phi, phi2] in degrees
            R = eulerToRotationMatrix(orientation(:)');

        else
            error('misorientation:invalidInput', ...
                'Orientation must be 3x3 matrix, 4-element quaternion, or 3-element Euler angles');
        end
    else
        error('misorientation:invalidInput', 'Unsupported orientation type');
    end
end

function R = quaternionToRotationMatrix(q)
    % Convert quaternion [w, x, y, z] to rotation matrix
    q = q(:) / norm(q);  % Normalize
    w = q(1); x = q(2); y = q(3); z = q(4);

    R = [1-2*(y^2+z^2),   2*(x*y-z*w),   2*(x*z+y*w);
         2*(x*y+z*w),   1-2*(x^2+z^2),   2*(y*z-x*w);
         2*(x*z-y*w),     2*(y*z+x*w), 1-2*(x^2+y^2)];
end

function q = rotationMatrixToQuaternion(R)
    % Convert rotation matrix to quaternion [w, x, y, z]
    t = trace(R);

    if t > 0
        s = 0.5 / sqrt(t + 1);
        w = 0.25 / s;
        x = (R(3,2) - R(2,3)) * s;
        y = (R(1,3) - R(3,1)) * s;
        z = (R(2,1) - R(1,2)) * s;
    elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
        s = 2 * sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        w = (R(3,2) - R(2,3)) / s;
        x = 0.25 * s;
        y = (R(1,2) + R(2,1)) / s;
        z = (R(1,3) + R(3,1)) / s;
    elseif R(2,2) > R(3,3)
        s = 2 * sqrt(1 + R(2,2) - R(1,1) - R(3,3));
        w = (R(1,3) - R(3,1)) / s;
        x = (R(1,2) + R(2,1)) / s;
        y = 0.25 * s;
        z = (R(2,3) + R(3,2)) / s;
    else
        s = 2 * sqrt(1 + R(3,3) - R(1,1) - R(2,2));
        w = (R(2,1) - R(1,2)) / s;
        x = (R(1,3) + R(3,1)) / s;
        y = (R(2,3) + R(3,2)) / s;
        z = 0.25 * s;
    end

    q = [w; x; y; z];
    if w < 0
        q = -q;  % Ensure positive scalar part
    end
end

function R = eulerToRotationMatrix(euler)
    % Convert Euler angles (Bunge convention, degrees) to rotation matrix
    % euler = [phi1, Phi, phi2]

    phi1 = deg2rad(euler(1));
    Phi = deg2rad(euler(2));
    phi2 = deg2rad(euler(3));

    c1 = cos(phi1); s1 = sin(phi1);
    c = cos(Phi);   s = sin(Phi);
    c2 = cos(phi2); s2 = sin(phi2);

    R = [c1*c2 - s1*s2*c,  s1*c2 + c1*s2*c,  s2*s;
        -c1*s2 - s1*c2*c, -s1*s2 + c1*c2*c,  c2*s;
         s1*s,            -c1*s,             c];
end

function euler = rotationMatrixToEuler(R)
    % Convert rotation matrix to Euler angles (Bunge convention, degrees)

    if abs(R(3,3)) < 1 - 1e-8
        Phi = acos(R(3,3));
        phi1 = atan2(R(3,1)/sin(Phi), -R(3,2)/sin(Phi));
        phi2 = atan2(R(1,3)/sin(Phi), R(2,3)/sin(Phi));
    else
        % Gimbal lock
        Phi = 0;
        phi1 = atan2(R(1,2), R(1,1));
        phi2 = 0;
    end

    euler = rad2deg([phi1, Phi, phi2]);

    % Ensure positive angles
    euler = mod(euler, 360);
end

function [angle, axis] = rotationMatrixToAngleAxis(R)
    % Convert rotation matrix to angle (degrees) and axis

    % Angle from trace
    cosAngle = (trace(R) - 1) / 2;
    cosAngle = max(-1, min(1, cosAngle));  % Clamp to [-1, 1]
    angle = rad2deg(acos(cosAngle));

    if angle < 1e-6
        % Identity rotation
        axis = [0; 0; 1];
    elseif abs(angle - 180) < 1e-6
        % 180 degree rotation - find axis from R + I
        [~, idx] = max(diag(R));
        axis = (R(:, idx) + eye(3, 1) * (idx == 1) + eye(3, 1) * (idx == 2) * [0;1;0] + eye(3, 1) * (idx == 3) * [0;0;1]);
        axis = axis / norm(axis);
    else
        % General case
        axis = [R(3,2) - R(2,3);
                R(1,3) - R(3,1);
                R(2,1) - R(1,2)];
        axis = axis / (2 * sin(deg2rad(angle)));
    end

    % Ensure axis points in positive direction (conventional)
    if axis(1) < 0 || (axis(1) == 0 && axis(2) < 0) || (axis(1) == 0 && axis(2) == 0 && axis(3) < 0)
        axis = -axis;
    end
end

%% Symmetry Operators
function symOps = getSymmetryOperators(crystalSystem)
    % Return symmetry operators for different crystal systems

    switch lower(crystalSystem)
        case 'cubic'
            % 24 proper rotations of cubic system
            symOps = getCubicSymmetry();

        case 'hexagonal'
            % 12 proper rotations
            symOps = getHexagonalSymmetry();

        case 'tetragonal'
            % 8 proper rotations
            symOps = getTetragonalSymmetry();

        case 'orthorhombic'
            % 4 proper rotations
            symOps = getOrthorhombicSymmetry();

        case 'monoclinic'
            % 2 proper rotations
            symOps = getMonoclinicSymmetry();

        case 'triclinic'
            % Only identity
            symOps = eye(3);
    end
end

function ops = getCubicSymmetry()
    % 24 proper rotations of cubic symmetry (point group 432)

    ops = zeros(3, 3, 24);
    idx = 1;

    % Identity
    ops(:,:,idx) = eye(3); idx = idx + 1;

    % 90, 180, 270 deg rotations about <100>
    for axis = 1:3
        for angle = [90, 180, 270]
            ops(:,:,idx) = axisAngleToMatrix(axis, angle);
            idx = idx + 1;
        end
    end

    % 120, 240 deg rotations about <111>
    diagonals = [1 1 1; 1 1 -1; 1 -1 1; -1 1 1];
    for d = 1:4
        for angle = [120, 240]
            axis = diagonals(d,:);
            ops(:,:,idx) = rodrigues(axis/norm(axis), deg2rad(angle));
            idx = idx + 1;
        end
    end

    % 180 deg rotations about <110>
    axes110 = [1 1 0; 1 -1 0; 1 0 1; 1 0 -1; 0 1 1; 0 1 -1];
    for a = 1:6
        axis = axes110(a,:);
        ops(:,:,idx) = rodrigues(axis/norm(axis), pi);
        idx = idx + 1;
    end
end

function ops = getHexagonalSymmetry()
    % 12 proper rotations of hexagonal symmetry

    ops = zeros(3, 3, 12);

    % Rotations about c-axis (z)
    for i = 0:5
        angle = i * 60;
        ops(:,:,i+1) = axisAngleToMatrix(3, angle);
    end

    % 180 deg rotations about a-axes
    for i = 0:5
        angle_a = i * 60;
        axis = [cosd(angle_a), sind(angle_a), 0];
        ops(:,:,7+i) = rodrigues(axis, pi);
    end
end

function ops = getTetragonalSymmetry()
    % 8 proper rotations

    ops = zeros(3, 3, 8);
    ops(:,:,1) = eye(3);

    % 90, 180, 270 about z
    ops(:,:,2) = axisAngleToMatrix(3, 90);
    ops(:,:,3) = axisAngleToMatrix(3, 180);
    ops(:,:,4) = axisAngleToMatrix(3, 270);

    % 180 about x, y
    ops(:,:,5) = axisAngleToMatrix(1, 180);
    ops(:,:,6) = axisAngleToMatrix(2, 180);

    % 180 about [110], [1-10]
    ops(:,:,7) = rodrigues([1 1 0]/sqrt(2), pi);
    ops(:,:,8) = rodrigues([1 -1 0]/sqrt(2), pi);
end

function ops = getOrthorhombicSymmetry()
    % 4 proper rotations

    ops = zeros(3, 3, 4);
    ops(:,:,1) = eye(3);
    ops(:,:,2) = axisAngleToMatrix(1, 180);
    ops(:,:,3) = axisAngleToMatrix(2, 180);
    ops(:,:,4) = axisAngleToMatrix(3, 180);
end

function ops = getMonoclinicSymmetry()
    % 2 proper rotations

    ops = zeros(3, 3, 2);
    ops(:,:,1) = eye(3);
    ops(:,:,2) = axisAngleToMatrix(2, 180);
end

function R = axisAngleToMatrix(axisIdx, angleDeg)
    % Rotation matrix for rotation about coordinate axis

    c = cosd(angleDeg);
    s = sind(angleDeg);

    switch axisIdx
        case 1  % x-axis
            R = [1 0 0; 0 c -s; 0 s c];
        case 2  % y-axis
            R = [c 0 s; 0 1 0; -s 0 c];
        case 3  % z-axis
            R = [c -s 0; s c 0; 0 0 1];
    end
end

function R = rodrigues(axis, angle)
    % Rodrigues formula for rotation matrix

    axis = axis(:) / norm(axis);
    K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
    R = eye(3) + sin(angle)*K + (1-cos(angle))*(K*K);
end

%% CSL Boundary Detection
function cslType = checkCSL(angle, axis, crystalSystem)
    % Check if misorientation corresponds to a CSL boundary

    cslType = '';

    if ~strcmpi(crystalSystem, 'cubic')
        return;  % CSL primarily defined for cubic
    end

    % Common CSL boundaries for cubic (angle in degrees, axis)
    cslData = {
        3,   60.00, [1 1 1]
        5,   36.87, [1 0 0]
        7,   38.21, [1 1 1]
        9,   38.94, [1 1 0]
        11,  50.48, [1 1 0]
        13,  22.62, [1 0 0]
        15,  48.19, [2 1 0]
        17,  28.07, [1 0 0]
        19,  26.53, [1 1 0]
        21,  21.79, [1 1 1]
        25,  16.26, [1 0 0]
        27,  31.59, [1 1 0]
        29,  43.60, [1 0 0]
        31,  17.90, [1 1 1]
    };

    tolerance = 2.0;  % degrees

    for i = 1:size(cslData, 1)
        sigma = cslData{i, 1};
        cslAngle = cslData{i, 2};
        cslAxis = cslData{i, 3};
        cslAxis = cslAxis / norm(cslAxis);

        % Check angle
        if abs(angle - cslAngle) < tolerance
            % Check axis (considering equivalent directions)
            axisDot = abs(dot(axis(:), cslAxis(:)));
            if axisDot > cosd(tolerance)
                cslType = sprintf('Sigma%d', sigma);
                return;
            end
        end
    end
end
