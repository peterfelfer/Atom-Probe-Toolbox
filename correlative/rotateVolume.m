function [volRot, voxelSizeOut, outShift] = rotateVolume(vol, voxelSize, axis, angle, varargin)
% rotateVolume rotates a 3D volume about an axis-angle pair.
%
% INPUTS (all in physical units):
%   vol        : 3D volume matrix.
%   voxelSize  : scalar or [x y z] voxel size.
%   axis       : 1x3 rotation axis.
%   angle      : rotation angle in radians.
%
% NAME-VALUE OPTIONS:
%   'Center'     : 1x3 rotation center (default: volume physical center).
%   'Shift'      : 1x3 shift of input grid origin (default: [0 0 0]).
%   'Method'     : 'linear' (default), 'nearest', or 'quadratic'.
%   'OutOfBounds': 'nan' (default), 'zero', or 'linear'.
%
% OUTPUTS:
%   volRot      : rotated volume sampled on an expanded grid.
%   voxelSizeOut: output voxel size (same as input).
%   outShift    : physical shift/origin of the output grid.

p = inputParser;
p.addParameter('Center', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 3));
p.addParameter('Shift', 0, @(x) isnumeric(x));
p.addParameter('Method', 'linear', @(x) ischar(x) || isstring(x));
p.addParameter('OutOfBounds', 'nan', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});

shift = p.Results.Shift;
center = p.Results.Center;
method = char(p.Results.Method);
outOfBounds = char(p.Results.OutOfBounds);

if ndims(vol) ~= 3
    error('rotateVolume expects a 3D volume.');
end

voxelSize = expandToDims(voxelSize, 3, 'voxelSize');
shift = expandToDims(shift, 3, 'shift');

axis = axis(:).';
if numel(axis) ~= 3 || any(~isfinite(axis))
    error('axis must be a finite 1x3 vector.');
end
axisNorm = norm(axis);
if axisNorm == 0
    error('axis must be non-zero.');
end
axis = axis / axisNorm;

nBinsIn = size(vol);
extent = nBinsIn .* voxelSize;
if isempty(center)
    center = shift + extent / 2;
else
    center = reshape(center, 1, 3);
end

inCenters = cell(1,3);
for d = 1:3
    inCenters{d} = (0:nBinsIn(d)-1) * voxelSize(d) + voxelSize(d)/2 + shift(d);
end

R = axisAngleToMatrix(axis, angle);

% Rotate the volume bounds (edges) to determine output grid extent
edges = [shift; shift + extent];
[xc, yc, zc] = ndgrid(edges(:,1), edges(:,2), edges(:,3));
cornerPts = [xc(:), yc(:), zc(:)];
rotCorners = (R * (cornerPts - center)')' + center;
minRot = min(rotCorners, [], 1);
maxRot = max(rotCorners, [], 1);

outShift = minRot;
outExtent = maxRot - minRot;

voxelSizeOut = voxelSize;

nBinsOut = max(1, ceil(outExtent ./ voxelSizeOut));

outCenters = cell(1,3);
for d = 1:3
    outCenters{d} = outShift(d) + (0:nBinsOut(d)-1) * voxelSizeOut(d) + voxelSizeOut(d)/2;
end

[xx, yy, zz] = ndgrid(outCenters{:});

% Inverse rotation to sample input grid
x = xx - center(1);
y = yy - center(2);
z = zz - center(3);

Rt = R.';
xIn = Rt(1,1)*x + Rt(1,2)*y + Rt(1,3)*z + center(1);
yIn = Rt(2,1)*x + Rt(2,2)*y + Rt(2,3)*z + center(2);
zIn = Rt(3,1)*x + Rt(3,2)*y + Rt(3,3)*z + center(3);

method = lower(strtrim(method));
switch method
    case 'linear'
        interpMethod = 'linear';
    case 'nearest'
        interpMethod = 'nearest';
    case 'quadratic'
        interpMethod = 'cubic';
    otherwise
        error('Unsupported Method. Use ''linear'', ''nearest'', or ''quadratic''.');
end

outOfBounds = lower(strtrim(outOfBounds));
switch outOfBounds
    case {'nan', 'none'}
        extrapMethod = 'none';
        fillValue = NaN;
    case {'zero', '0'}
        extrapMethod = 'none';
        fillValue = 0;
    case {'linear', 'extrap'}
        extrapMethod = 'linear';
        fillValue = NaN;
    otherwise
        error('Unsupported OutOfBounds. Use ''nan'', ''zero'', or ''linear''.');
end

F = griddedInterpolant(inCenters, double(vol), interpMethod, extrapMethod);
volRot = F(xIn, yIn, zIn);

if strcmp(extrapMethod, 'none')
    if ~isnan(fillValue)
        volRot(isnan(volRot)) = fillValue;
    end
end

end

function vec = expandToDims(vec, numberDims, name)
if isscalar(vec)
    vec = repmat(vec, 1, numberDims);
elseif numel(vec) ~= numberDims
    error('%s must be a scalar or %d-element vector.', name, numberDims);
else
    vec = reshape(vec, 1, numberDims);
end
end

function R = axisAngleToMatrix(axis, angle)
% Rodrigues' rotation formula
kx = axis(1);
ky = axis(2);
kz = axis(3);

c = cos(angle);
s = sin(angle);

K = [  0   -kz   ky;
      kz    0   -kx;
     -ky   kx    0 ];

R = eye(3)*c + (1 - c)*(axis'*axis) + K*s;
end
