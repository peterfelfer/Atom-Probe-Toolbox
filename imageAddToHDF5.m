function imageAddToHDF5(fileName,image,imageType,contrastMechanism,horizontalPixelScale,verticalPixelScale)
% image is an 2D Array integer
%
%
%
%
%
%
%
%
%


% general base path for images of the sample
BASEPATH = "/sample/images/";

% check how many images are already in the hdf5 file
info = h5info(fileName);
isSample = string(extractfield(info.Groups,'Name')) == '/sample';
sampleGroup = info.Groups(isSample).Groups;
try
    numImages = length(sampleGroup.Datasets);
catch
    numImages = 0;
end
path = BASEPATH + num2str(numImages + 1);

% Writing image at next available index
h5create(fileName,path,size(image));
h5write(fileName,path,image);

% writing metadata
h5writeatt(fileName,path,"imageType",imageType);
h5writeatt(fileName,path,"contrastMechanism",contrastMechanism);
h5writeatt(fileName,path,"horizontalPixelScale",horizontalPixelScale);
h5writeatt(fileName,path,"verticalPixelScale",verticalPixelScale);