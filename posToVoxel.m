function vox = posToVoxel(pos,gridVec,species)
% posToVoxel creates a voxelisation of the data in 'pos' based on the bin centers in
% 'gridVec' for the atoms/ions in species. 
% 
% vox = posToVoxel(pos,gridVec,species)
% vox = posToVoxel(pos,gridVec)
%
% INPUT
% pos:     pos file; when input species is given, ranges must be allocated.
%          A decomposed pos file is also possible.
%
% gridVec: are the gridVectors for the grid. If the grid vectors need to be
%          created, use binCenters (output from binVectorsFromDistance
%          function)
%
% species: can be a species list as in {'Fe', 'Mn'}. Ions or atoms 
%          in pos can be parsed. Alternatively, it can be a logical vector 
%          with the same length as the pos file.
%          If no species variable is given, all atoms/ions are taken.
%          To get all ranged atoms, use species = categories(pos.atom)
%          To get all ranged ions, use species = categories(pos.ion)
%
% OUTPUT
% vox:     voxelisation of the point cloud stored in pos
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% Check for species
if ~exist('species','var')
    species = true(height(pos),1);
end


if iscell(species) % create logical vector with true for all elements to be included
    % check if atomic or ionic count is calculated
    if any(ismember(pos.Properties.VariableNames,'atom'))
        species = ismember(pos.atom,species);
        
    elseif any(ismember(pos.Properties.VariableNames,'ion'))
        species = ismember(pos.ion,species);
        
    else
        error('unknown table format');
    end   
end

pos = [pos.x(species), pos.y(species), pos.z(species)];

%% calculating 3d histogram

bin = [gridVec{1}(2)-gridVec{1}(1) gridVec{2}(2)-gridVec{2}(1) gridVec{3}(2)-gridVec{3}(1)];

% calculating bin association (I do not pretend to know what all of this
% does)
for d = 1:3
    % calculate edge vectors from bin centers
    edgeVec{d} = [gridVec{d}-bin(d)/2 gridVec{d}(end)+bin(d)/2];
    
    [useless loc(:,d)] = histc(pos(:,d),edgeVec{d},1);
    sz(d) = length(edgeVec{d})-1;
end
clear useless

% for points on the border
sz = max([sz; max(loc,[],1)]);

% count of atoms in each voxel
vox = accumarray(loc,1,sz);
end
