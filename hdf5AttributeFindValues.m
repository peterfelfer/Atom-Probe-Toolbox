function values = hdf5AttributeFindValues(fileList,group,attributeName)
% hdf5AttributeFindValues retrives all unique values of the attribute with
% the name attributeName in the given group for all files in the fileList
%
%   values = hdf5AttributeFindValues(fileList,group,attributeName)
%
%   INPUT:
%   fileList        list of hdf5 files that contain the path and attribute
%
%   group           group that contains the attribute. E.g. '/sample/'
%
%   attributeName   name of the attribute to find e.g. 'uniqueIdentifier'