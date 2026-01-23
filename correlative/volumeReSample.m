function vol = volumeReSample(vol, voxelSizeIn, binSize, shift, method)
% volumeReSample creates a new resampled volume from a 3D volume matrix.
% The voxel sizes are specified in physical units as scalars (isotropic) or
% vectors [xbin ybin zbin]. The output preserves the physical extent of the
% input volume (using ceil for the output grid size). A shift can be applied
% to the input grid (in physical units) to keep a global coordinate system.
% method sets the interpolation method for griddedInterpolant.

if nargin < 4 || isempty(shift)
    shift = 0;
end
if nargin < 5 || isempty(method)
    method = 'linear';
end

numberDims = ndims(vol);
if numberDims < 3
    error('volumeReSample expects a 3D volume.');
end

voxelSizeIn = expandToDims(voxelSizeIn, numberDims, 'voxelSizeIn');
binSize = expandToDims(binSize, numberDims, 'binSize');
shift = expandToDims(shift, numberDims, 'shift');
method = validatestring(method, {'linear','nearest','cubic','spline','makima'}, ...
    mfilename, 'method');

nBinsIn = size(vol);
extent = nBinsIn .* voxelSizeIn;
nBinsOut = ceil(extent ./ binSize);

inCenters = cell(1, numberDims);
outCenters = cell(1, numberDims);
for d = 1:numberDims
    inCenters{d} = (0:nBinsIn(d)-1) * voxelSizeIn(d) + voxelSizeIn(d)/2 + shift(d);
    outCenters{d} = (0:nBinsOut(d)-1) * binSize(d) + binSize(d)/2;
end

% Linear interpolation with linear extrapolation provides gradient-based
% behavior outside the input grid. Any NaNs from missing input data are
% filled afterwards via linear extrapolation along each dimension.
F = griddedInterpolant(inCenters, double(vol), method, 'linear');
[queryGrid{1:numberDims}] = ndgrid(outCenters{:}); %#ok<CCAT>
vol = F(queryGrid{:});

if any(isnan(vol(:)))
    for d = 1:numberDims
        vol = fillmissing(vol, 'linear', d, 'EndValues', 'extrap');
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
