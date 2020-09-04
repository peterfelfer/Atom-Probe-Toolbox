function distance = posDistanceToMesh(pos,interface)
% posDistanceToMesh calculates the distance of individual atoms to a mesh
% representing e.g. an interface or any other object that can be displayed
% as a triangulation 
% fv is a triangulation consisting of faces and vertices (Matlab 'patch')
%
% distance = posDistanceToMesh(pos,interface)
%
% INPUT
% pos:          table, pos file with pos.x, pos.y, and pos.z 
%
% interface:    structure, interface or any other object with fields of faces and vertices
%
% OUTPUT
% distance:     Mx1 array with M as the number of atoms, distance of the atoms in the pos file to the interface mesh

%%
% distances are calculated along vertex normals.
normals = patchnormals(interface);

% tessellation and distance calculation
% finding closest point for each atomic position
closestVertex = dsearchn(interface.vertices,delaunayn(interface.vertices),[pos.x, pos.y, pos.z]);

% distance vector to the nearest interface vertex
distVec = [pos.x, pos.y, pos.z] - interface.vertices(closestVertex,:);

% distance along normal through dot product
distance = sum(normals(closestVertex,:) .* distVec,2);