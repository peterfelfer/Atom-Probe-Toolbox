function [meshObject, WCS] = patchExtractFromHDF5(fileName,objectName)
% patchExtractFromHdf5 retrieves a mesh object and in the case of an ROI
% also the corresponding work coordinate system from an hdf5 file. The mesh
% object can be added to the file using hdf5AddMeshRoi.
%
% patchExtractFromHdf5(fileName,objectName)
%
% INPUTS
% fileName:     full file name of existing hdf5 file
%
% objectName:   char array (string) of the object name as saved to
%               the HDF5 file, e.g. 'ROI box'
%
% meshObject:   object stored in the hdf5 file. The object
%               may be an nD mesh (mesh, lines, points) as a patch object
%               The object may contain the following datasets:
%               vertices... Ix3 list of coordinates of the mesh vertices
%               faces... NxM (M = 2..4) list of quadrilaterals, triangles, or
%               lines
%               faceColor... 1x3 (single color) or Nx3 (per face color)
%               faceAlpha... 1x1 (single opacity value) or Nx1 (per face
%               opacity value)
%
% OUTPUT
% meshObject:   struct of mesh object with field of faces and vertices
% WCS:          struct defining a work coordinate system. Fields are
%               'origin' (1x3 vector) and 'coordinateSystemVectors' (3x3 list
%               of 3 basis vectors of coordinate system (exx exy exz; eyx eyy eyz; ezx ezy ezz)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

% get list of properties saved for the object
path = ['/objects/meshObject/' objectName '/'];

info = h5info(fileName,path);

dataList = struct2cell(info.Datasets)';
dataList = dataList(:,1);

% go through list of individual datasets
for d = 1:length(dataList)
    switch dataList{d}
        
        case 'faces'
            meshObject.faces = h5read(fileName, [path dataList{d}]);
        case 'vertices'
            meshObject.vertices = h5read(fileName, [path dataList{d}]);
        case 'faceAlpha'
            meshObject.faceAlpha = h5read(fileName, [path dataList{d}]);
        case 'faceColor'
            meshObject.faceColor = h5read(fileName, [path dataList{d}]);
        case 'origin'
            WCS.origin = h5read(fileName, [path dataList{d}]);
        case 'coordinateSystemVectors'
            WCS.coordinateSystemVectors = h5read(fileName, [path dataList{d}]);
        otherwise
            error(['mesh object property ' dataList{d} ' not recognized']);
    end
    
end
