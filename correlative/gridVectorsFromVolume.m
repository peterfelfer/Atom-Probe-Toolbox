function [binCenters, binEdges] = gridVectorsFromVolume(vol, binSize)
% gridVectorsFromVolume creates grid vectors in the form of gv{1}... for
% the volume in vol, with a specified binSize. This can be a scalar or
% vector with the dimension of vol. vol is a scalar field

numberDims = ndims(vol); %number of dimensions
nBins = size(vol); % number of bins in each dimension

if isscalar(binSize)
    binSize = repmat(binSize,numberDims,1);
end

for d = 1:numberDims
    binCenters{d} = (0:nBins(d)-1)*binSize(d) + binSize(d)/2;
    binEdges{d} = (0:nBins(d))*binSize(d);
end