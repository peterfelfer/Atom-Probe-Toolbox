function map = patchCreate2DdensityMap(pos,fv,d)
% NEEDS DOCUMENTATION
% calculates a pseudo-2D concentration map for a species pos in an atom
% probe dataset parPos
% its like a 2D concentration profile, but for general manifolds and is
% calculated within a clipping distance d. If d is a 2 part vector its from
%d(1) to d(2)

numVerts = length(fv.vertices);

if length(d) == 2
    sort(d);
else
    d = [-d/2, d/2];
end


% variables needed
pos = pos(:,1:3);
normals = patchnormals(fv);



%% calculate per vertex count atoms
% calculate which vertex is closest to an atom
closest = dsearchn(fv.vertices,pos);
distVec = pos(:,1:3) - fv.vertices(closest,:);
dist = sum(normals(closest,:) .* distVec,2);


% clip by distance
isIn = dist > d(1) & dist < d(2);
closest = closest(isIn);



% count per vertex
for v = 1:numVerts
    % raw count of IE atoms per vertex
    solCount(v) = sum(closest == v);
end

densMap = solCount;

%visualise
figure; trisurf(fv.faces,fv.vertices(:,1),fv.vertices(:,2),fv.vertices(:,3),densMap); axis equal; rotate3d on;

map = densMap;
end