% Basic tests for rotateVolume

clear; clc;
doVisualize = true;

% Gyroid field for visual inspection (periodic, easy to see binning)
gyroid = @(x,y,z, s) sin(x/s).*cos(y/s) + sin(y/s).*cos(z/s) + sin(z/s).*cos(x/s);

%% Test 1: isotropic rotation about z-axis
volSize = [24 24 24];
voxelSize = 1;
shift = [0 0 0];
axis = [0 0 1];
angle = pi/3;
s = 4.0;

inCenters = cell(1,3);
for d = 1:3
    inCenters{d} = (0:volSize(d)-1) * voxelSize + voxelSize/2 + shift(d);
end
[xi, yi, zi] = ndgrid(inCenters{:});
vol = gyroid(xi, yi, zi, s);

[volRot, voxelSizeOut, outShift] = rotateVolume(vol, voxelSize, axis, angle, ...
    'Shift', shift, 'Method', 'linear', 'OutOfBounds', 'nan');

assert(all(voxelSizeOut == voxelSize), 'Output voxel size should match input.');
assert(all(size(volRot) >= volSize), 'Output size should be >= input size for expansion.');

% Build expected volume on the output grid using the analytical gyroid
outCenters = cell(1,3);
for d = 1:3
    outCenters{d} = outShift(d) + (0:size(volRot,d)-1) * voxelSizeOut(d) + voxelSizeOut(d)/2;
end
[xx, yy, zz] = ndgrid(outCenters{:});

extent = volSize .* voxelSize;
center = shift + extent / 2;
R = axisAngleToMatrix(axis / norm(axis), angle);
Rt = R.';

x = xx - center(1);
y = yy - center(2);
z = zz - center(3);

xIn = Rt(1,1)*x + Rt(1,2)*y + Rt(1,3)*z + center(1);
yIn = Rt(2,1)*x + Rt(2,2)*y + Rt(2,3)*z + center(2);
zIn = Rt(3,1)*x + Rt(3,2)*y + Rt(3,3)*z + center(3);

% Mask points outside the input bounds for NaN comparison
minIn = shift;
maxIn = shift + extent;
inside = xIn >= minIn(1) & xIn <= maxIn(1) & ...
         yIn >= minIn(2) & yIn <= maxIn(2) & ...
         zIn >= minIn(3) & zIn <= maxIn(3);

expected = NaN(size(volRot));
expected(inside) = gyroid(xIn(inside), yIn(inside), zIn(inside), s);

mask = ~isnan(expected) & ~isnan(volRot);
maxErr = max(abs(volRot(mask) - expected(mask)));
meanErr = mean(abs(volRot(mask) - expected(mask)));
assert(maxErr < 0.25 && meanErr < 0.05, 'Rotation mismatch in Test 1.');

if doVisualize
    showVolumePair(vol, volRot, 'RotateVolume Test 1');
end

fprintf('All rotateVolume tests passed.\n');

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

function showVolumePair(volIn, volOut, titlePrefix)
if exist('volshow','file') == 2
    showVol(volIn, [titlePrefix ' Input']);
    showVol(volOut, [titlePrefix ' Output']);
else
    warning('volshow not found; skipping visualization.');
end
end

function showVol(volData, figName)
try
    displayData = volData;
    if any(~isfinite(displayData(:)))
        displayData(~isfinite(displayData)) = 0;
    end
    volshow(displayData, 'RenderingStyle', 'MaximumIntensityProjection');
    fig = gcf;
    fig.Name = figName;
    fig.Color = 'w';
catch
    warning('volshow failed in this MATLAB version; skipping visualization.');
end
end
