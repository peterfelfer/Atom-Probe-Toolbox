function imageAddToHDF5(fileName,image,imageType,contrastMechanism,horizontalPixelScale,verticalPixelScale)
% imageAddToHDF5 adds an image file to the hdf5 file
%
% imageAddToHDF5(fileName,image,imageType,contrastMechanism,horizontalPixelScale,verticalPixelScale);
%
% INPUT
%
% fileName:             name of the already existing .h5 file
% 
% image:                variable form the workspace that contains the image
%                       stored as array integer. The user can use imread() 
%                       to load the image into the workspace
% 
% imageType:            type of the image, for example 'SEM image'
% 
% contrastMechanism:    the contrast mechanism with whom the image was
%                       acquired, e.g. 'secondary electron' for an image
%                       from a SEM
% 
% horizontalPixelScale: is the pixel scale in horizontal direction
% 
% verticalPixelScale:   is the pixel scale in vertical direction



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