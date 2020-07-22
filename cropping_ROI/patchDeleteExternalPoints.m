function fv = patchDeleteExternalPoints(gh,pos,numRims)
% patchDeleteExternalPoints deletes all points of a patch object that are outside the dataset. This
% is done by deleting all points that have no atomic position that is
% closer to them than no other vertex.
% 
% fv = patchDeleteExternalPosnts(fv,pos,numRims)
% 
% INPUT
% gh:       graphic handle that is needed to biult the structure of
%           vertices and faces
%
% pos:      dataset on which the patch object will be fitted
%
% numRims:  numRims is the number of edge loops to be deleted on top of the unused verts
% 
% OUTPUT
% points:   points transformed to the new coordinate system   


%% creating the structure
fv = struct();
fv.vertices = gh.Vertices;
fv.faces = gh.Faces;
%% calculate which vertices are inside the dataset
vertCoords = fv.vertices;
atomCoords = [pos.x, pos.y, pos.z];
closest = dsearchn(vertCoords,atomCoords);
isIn = ismember(1:length(vertCoords),closest);


%% remove unused vertices
fv = patchRemoveVertices(fv,fv.vertices(~isIn,:));

%% remove boundary vertices for n boundary rings
if exist('numRims','var')
    for i=1:numRims
        bnd = computeBoundary(fv);
        fv = patchRemoveVertices(fv,fv.vertices(bnd,:));
    end
end