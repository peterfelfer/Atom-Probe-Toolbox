function rotatedPoints = rotatePointsAxisAngle(points, axis, angle)
    % rotatePointsAxisAngle rotates a set of 3D points about a given axis by a given angle.
    %
    % INPUTS:
    %   points : An N-by-3 matrix of 3D points to rotate. 
    %            Each row is a point [x, y, z].
    %   axis   : A 1-by-3 vector specifying the axis of rotation [ax, ay, az].
    %   angle  : Rotation angle in radians.
    %
    % OUTPUT:
    %   rotatedPoints : The rotated points as an N-by-3 matrix.
    %
    % EXAMPLE:
    %   % Rotate points around the z-axis by 90 degrees (pi/2 radians)
    %   pts = [1 0 0; 0 1 0; 0 0 1];
    %   axis = [0 0 1];
    %   angle = pi/2;
    %   pts_rot = rotatePointsAxisAngle(pts, axis, angle);

    % Ensure axis is a unit vector
    axis = axis / norm(axis);

    % Extract components of axis
    kx = axis(1);
    ky = axis(2);
    kz = axis(3);

    % Compute sine and cosine of angle
    c = cos(angle);
    s = sin(angle);

    % Rodrigues' rotation formula:
    % R = I*c + (1 - c)*k*k' + [k]_x * s
    %
    % where k is the axis and [k]_x is the skew-symmetric matrix:
    % [k]_x = [  0   -kz   ky
    %           kz     0  -kx
    %          -ky    kx    0 ]

    K = [  0   -kz   ky;
          kz    0   -kx;
         -ky   kx    0 ];

    % Rotation matrix
    R = eye(3)*c + (1 - c)*(axis'*axis) + K*s;

    % Apply rotation
    rotatedPoints = (R * points')';  % Rotate and transpose back

end
