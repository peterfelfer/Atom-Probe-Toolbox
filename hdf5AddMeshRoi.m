function hdf5AddMeshRoi(fileName,meshOrRoi,objectName)
% hdf5AddMeshRoi adds a mesh such as an isosurface or a ROI to the hdf5
% file. This can be either a graphics object handle or a patch object.
%
% USAGE:    hdf5AddMeshRoi(fileName,meshOrRoi, objectName)
%           hdf5AddMeshRoi(fileName,meshOrRoi)
%
%
% INPUTS:   fileName = full file name of existing hdf5 file
%
%           meshOrRoi = object to be stored in the hdf5 file. The object
%           may be an nD mesh (mesh, lines, points) as a patch object or
%           via matlab graphics handle.
%           patch: matlab struct with meshOrRoi.vertices and
%           meshOrRoi.faces. 
%           graphics handle representing a mesh. If such as in the case of
%           ROIs a work coordinate system exists in
%           meshOrRoi.UserData.ROIxAxis etc. it is stored alongside.
%           
%           objectName = name of a mesh object that parsed as patch or name
%           to override the name in meshOrRoi.DisplayName (property of
%           matlab graphics objects)

% check if mesh object already exists

basePath = '/objects/meshObject/';

try 
    info = h5info(fileName,basePath);
    numMeshObjects = length(info.Groups);
catch
    numMeshObjects = 0;
    fid = H5F.open(fileName);
        groupID = H5G.create(fid,basePath,...
        'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(groupID);
    H5F.close(fid);
end

% index of object to be created
numCurrObj = numMeshObjects + 1;
%objPath = [basePath num2str(numCurrObj)];


% extract data from objects
if isgraphics(meshOrRoi)
    
    % define mesh data
    vertices = meshOrRoi.Vertices;
    faces = meshOrRoi.Faces;
    faceColor = meshOrRoi.FaceColor;
    faceAlpha = meshOrRoi.FaceAlpha;
    
    % get name
    if exist('objectName','var')
        name = objectName;
    else
        name = meshOrRoi.DisplayName;
    end
    
    objPath = [basePath name];
    
        % save to hdf5
    h5create(fileName,[objPath '/vertices'], size(vertices));
    h5write(fileName,[objPath '/vertices'], vertices);
    
    h5create(fileName,[objPath '/faces'], size(faces));
    h5write(fileName,[objPath '/faces'], faces);
    
    h5create(fileName,[objPath '/faceColor'], size(faceColor));
    h5write(fileName,[objPath '/faceColor'], faceColor);
    
    h5create(fileName,[objPath '/faceAlpha'], size(faceAlpha));
    h5write(fileName,[objPath '/faceAlpha'], faceAlpha);
    
    h5writeatt(fileName,objPath,'name',name);
    
    
    % get coordinates for WCS
    if isfield(meshOrRoi.UserData,'ROIxaxis')
        origin = meshOrRoi.UserData.ROIxaxis(1,:);
        coordinateSystemVectors = [meshOrRoi.UserData.ROIxaxis(2,:) - meshOrRoi.UserData.ROIxaxis(1,:);...
            meshOrRoi.UserData.ROIyaxis(2,:) - meshOrRoi.UserData.ROIyaxis(1,:);...
            meshOrRoi.UserData.ROIzaxis(2,:) - meshOrRoi.UserData.ROIzaxis(1,:)];
        
        h5create(fileName,[objPath '/origin'], size(origin));
        h5write(fileName,[objPath '/origin'], origin);
        
        h5create(fileName,[objPath '/coordinateSystemVectors'], size(coordinateSystemVectors));
        h5write(fileName,[objPath '/coordinateSystemVectors'], coordinateSystemVectors);
        
    end
    

    
    
elseif isstruct(meshOrRoi)
    vertices = meshOrRoi.vertices;
    faces = meshOrRoi.faces;
    name = objectName;
    
    objPath = [basePath name];
    
    h5create(fileName,[objPath '/vertices'], size(vertices));
    h5write(fileName,[objPath '/vertices'], vertices);
    
    h5create(fileName,[objPath '/faces'], size(faces));
    h5write(fileName,[objPath '/faces'], faces);
    
    h5writeatt(fileName,objPath,'name',name);
end