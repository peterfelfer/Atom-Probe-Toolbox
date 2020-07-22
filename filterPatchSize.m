function fv = filterPatchSize(fv,Nmin)
% filterPatchSize filters the individually connected parts of a patch and 
% discards patch parts smaller than Nmin. Useful for filtering noise in 
% isosurfaces
% 
% 
% 

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
    

