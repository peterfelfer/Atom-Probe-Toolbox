function in = posInVertexVoronoiCell(pos,interface,vertexIdx)
% posInVertexVoronoiCell either gives a true if an atomic position is
% closer to the vertices in vertexIdx than any other vertices or if no
% vertexIdx is parsed, outputs the index of the closest vertex

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
    in = ismember(Pclosest,vertexIdx);
else
    in = Pclosest;
end