function patchCreateExcessValue(pos,parentPos,interface,vertices)
% patchCreateExcessValue calculates the average interfacial excess for an 
% interface with an interactive interface for the interfacial excess 
% determination 
% reads a pos file (with x, y, and z) and a vertex file 
% with x, y, and z) and assigns every atom to the closest vertex
%
% patchCreateExcessValue(pos,parentPos,interface)
% patchCreateExcessValue(pos,parentPos,interface,vertices)
%
% INPUT
% pos:          table, pos file of the atom species, of which an interfacial excess 
%               should be calculated
%
% parentPos:    table, pos file of all atom species; optimally all unranged
%               atoms are excluded
%
% interface:    structure with fields of faces (f) and vertices (v)
%
% vertices:     optional, if the interfacial excess of a specific region is
%               wanted
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

addpath('dualMesh');
addpath('patch_normals');
addpath('resources');

normals = patchnormals(interface);

numVerts = length(interface.vertices);
%% tessellation and distance clipping
% finding closest point for each atomic position
closest = dsearchn(interface.vertices,delaunayn(interface.vertices),[parentPos.x, parentPos.y, parentPos.z]);

% calculation of vertex area
[cp,ce,pv,ev] = makedual2(interface.vertices,interface.faces);
[pc,area] = geomdual2(cp,ce,pv,ev);


% if local IE values are calculated, only atomic positions associated with
% the vertices on the vertex list are used
if exist('vertices','var')
    isLocal = ismember(closest,vertices);
    parentPos = parentPos(isLocal,:);
    closest = closest(isLocal,:);
    
    area = area(vertices);
end


distVec = [parentPos.x, parentPos.y, parentPos.z] - interface.vertices(closest,:);
% distance through dot product
dist = sum(normals(closest,:) .* distVec,2);
area = sum(area);






%% calculating cumulative diagram

% indices of atoms that are part of the species in question
idx = ismember([parentPos.x, parentPos.y, parentPos.z],[pos.x, pos.y, pos.z],'rows');
[useless, sortIdx] = sort(dist);
idx = idx(sortIdx);

cumulative = cumsum(idx);

%% index of interface location
interfaceLoc = median(find(abs(dist(sortIdx)) == min(abs(dist))));


excessGUI(cumulative,area,interfaceLoc);



end
