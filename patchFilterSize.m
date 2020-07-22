function fv = patchFilterSize(fv,Nmin)
% patchFilterSize filters the individually connected parts of a patch and 
% discards patch parts smaller than Nmin. Useful for filtering noise in 
% isosurfaces
% 
% fv = filterPatchSize(fv,Nmin)
% 
% INPUT
% fv =      patch with vertices and faces
% Nmin =    is the boundary, all patches that are smaller than Nmin are
%           discarded
%
% OUTPUT
% fv =      residua patches that are bigger than Nmin

patches = splitFV(fv);
numPatches = length(patches);

%% deleting small patches

isSmall = false(numPatches,1);
 
for p = 1:numPatches
    isSmall(p) = length(patches(p).vertices) < Nmin;
end

patches(isSmall) = [];

%%reassembly of patches
numPatches = sum(~isSmall);

fv.vertices = [];
fv.faces = [];

numVerts = 0;

for p = 1:numPatches
    fv.vertices = [fv.vertices; patches(p).vertices];
    fv.faces = [fv.faces; patches(p).faces + numVerts];
    numVerts = numVerts + length(patches(p).vertices(:,1));
end
    

