function inOrIdx = posInVertexVoronoiCell(pos,interface,vertexIdx)
% posInVertexVoronoiCell either gives a true if an atomic position is
% closer to the vertices in vertexIdx than any other vertices or if no
% vertexIdx is parsed, outputs the index of the closest vertex
%
% inOrIdx = posInVertexVoronoiCell(pos,interface,vertexIdx)
% inOrIdx = posInVertexVoronoiCell(pos,interface)
%
% INPUTS
% pos:          pos variable
% interface:    mesh variable (patch) with interface.vertices and
%               interface.faces
% vertexIdx:    indices of vertices in the interface to which the
%               association is tested.
%
% OUTPUTS
% inOrIdx:      if vertexIdx is parsed, it is a logical variable which is
%               true if an atom is closer to a vertex in the vertex list 
%               than any other vertex
%               If no vertexIdx are parsed, it gives the index of the
%               closest vertex to each atomic position.
%
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


if istable(pos)
    xyz = [pos.x, pos.y, pos.z];
end

% finding closest point for each atomic position
try
    Pclosest = dsearchn(interface.vertices,delaunayn(interface.vertices),xyz(:,1:3));
    
catch
    % for very regular meshes, there is no unique solution. Therefore we
    % wash out the interface position
    offset = (rand(size(interface.vertices))*2 - 1) * 0.01;
    interface.vertices = interface.vertices + offset;
    Pclosest = dsearchn(interface.vertices,delaunayn(interface.vertices+offset),xyz(:,1:3));
end

if exist('vertexIdx','var')
    inOrIdx = ismember(Pclosest,vertexIdx);
else
    inOrIdx = Pclosest;
end