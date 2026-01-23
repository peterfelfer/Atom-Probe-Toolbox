% Basic tests for volumeReSample

clear; clc;
doVisualize = true;

% Gyroid field for visual inspection (periodic, easy to see binning)
gyroid = @(x,y,z, s) sin(x/s).*cos(y/s) + sin(y/s).*cos(z/s) + sin(z/s).*cos(x/s);

%% Test 1: isotropic resample, no shift
volSize = [20 22 24];
voxelSizeIn = 1; % physical units
binSize = 1;
shift = 0;
s = 3.5;

inCenters = cell(1,3);
for d = 1:3
    inCenters{d} = (0:volSize(d)-1) * voxelSizeIn + voxelSizeIn/2 + shift;
end
[xi, yi, zi] = ndgrid(inCenters{:});
vol = gyroid(xi, yi, zi, s);

volOut = volumeReSample(vol, voxelSizeIn, binSize, shift, 'linear');

extent = volSize .* voxelSizeIn;
expectedSize = ceil(extent ./ binSize);
assert(isequal(size(volOut), expectedSize), 'Size mismatch in Test 1.');

expected = vol;
maxErr = max(abs(volOut(:) - expected(:)));
assert(maxErr < 1e-10, 'Interpolation mismatch in Test 1.');
if doVisualize
    showVolumeTriplet(vol, volOut, expected, 'Test 1');
end

%% Test 2: anisotropic resample with shift
volSize = [21 24 20];
voxelSizeIn = [1 2 3];
binSize = [2 2 2];
shift = [0.5 1.0 1.5];
s = 4.0;

inCenters = cell(1,3);
for d = 1:3
    inCenters{d} = (0:volSize(d)-1) * voxelSizeIn(d) + voxelSizeIn(d)/2 + shift(d);
end
[xi, yi, zi] = ndgrid(inCenters{:});
vol = gyroid(xi, yi, zi, s);

% Inject NaNs to exercise fillmissing
vol(2,3,2) = NaN;

volOut = volumeReSample(vol, voxelSizeIn, binSize, shift, 'linear');

extent = volSize .* voxelSizeIn;
expectedSize = ceil(extent ./ binSize);
assert(isequal(size(volOut), expectedSize), 'Size mismatch in Test 2.');

outCenters = cell(1,3);
for d = 1:3
    outCenters{d} = (0:expectedSize(d)-1) * binSize(d) + binSize(d)/2;
end
[xo, yo, zo] = ndgrid(outCenters{:});
expected = gyroid(xo, yo, zo, s);
maxErr = max(abs(volOut(:) - expected(:)));
meanErr = mean(abs(volOut(:) - expected(:)));
assert(maxErr < 0.25 && meanErr < 0.05, 'Interpolation mismatch in Test 2.');

assert(~any(isnan(volOut(:))), 'NaNs remain after resampling.');
if doVisualize
    showVolumeTriplet(vol, volOut, expected, 'Test 2');
end

%% Test 3: rotation utility
pts = [1 0 0; 0 1 0; 0 0 1];
axis = [0 0 1];
angle = pi/2;
expected = [0 1 0; -1 0 0; 0 0 1];
rot = rotatePointsAxisAngle(pts, axis, angle);
assert(max(abs(rot(:) - expected(:))) < 1e-12, 'Rotation mismatch for z-axis 90deg.');

axis = [1 1 1];
axis = axis / norm(axis);
ptOnAxis = axis;
rot = rotatePointsAxisAngle(ptOnAxis, axis, pi/3);
assert(norm(rot - ptOnAxis) < 1e-12, 'Point on axis should remain unchanged.');

pts = [1 2 3; -2 1 0.5];
rot = rotatePointsAxisAngle(pts, [1 0 0], pi/4);
assert(max(abs(vecnorm(rot,2,2) - vecnorm(pts,2,2))) < 1e-12, ...
    'Rotation should preserve norms.');

fprintf('All volumeReSample tests passed.\n');

function showVolumeTriplet(volIn, volOut, expected, titlePrefix)
if exist('volshow','file') == 2
    showVol(volIn, [titlePrefix ' Input']);
    showVol(volOut, [titlePrefix ' Output']);
    showVol(abs(volOut - expected), [titlePrefix ' |Diff|']);
else
    warning('volshow not found; using slice visualization as fallback.');
    showSlices(volIn, [titlePrefix ' Input']);
    showSlices(volOut, [titlePrefix ' Output']);
    showSlices(abs(volOut - expected), [titlePrefix ' |Diff|']);
end
end

function showVol(volData, figName)
try
    volshow(volData, 'RenderingStyle', 'MaximumIntensityProjection');
    fig = gcf;
    fig.Name = figName;
    fig.Color = 'w';
catch
    warning('volshow failed in this MATLAB version; using slice visualization.');
    showSlices(volData, figName);
end
end

function showSlices(volData, figName)
sz = size(volData);
fig = figure('Name', figName, 'Color', 'w');
slice(volData, round(sz(2)/2), round(sz(1)/2), round(sz(3)/2));
shading interp;
axis tight vis3d off;
colormap(fig, parula);
end
