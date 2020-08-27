function image = imageFromHDF5(fileName, numImg)
% imageFromHDF5 extract the image with number numImg from the HDF5 file
%
% image = imageFromHDF5(fileName, numImg)
%
% INPUT
% fileName: Name of the HDF5 file
%
% numImg:   Number of the image that should be extracted
%
% OUTPUT
% image:    variable that contains the image


% extract an image from HDF5 file 
numImgStr = num2str(numImg);
addressData = ['/sample/images/' numImgStr];
imageRaw = h5read(fileName, addressData);
image = mat2gray(imageRaw);


end

