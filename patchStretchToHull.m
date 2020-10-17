function fv = patchStretchToHull(fv,boundary)
% NEEDS DOCUMENTATION
%stretches the boundary of a mesh to a specified boundary hull mesh

DEBUG = false;

addpath('patch_normals');
addpath('TriangleRayIntersection');

%% determine which vertices belong to the boundary
% those who belong to an edge that only belongs to one triangle

% produce edge list
ed = [fv.faces(:,[1 2]); fv.faces(:,[2 3]); fv.faces(:,[1 3])];
ed = sort(ed,2);


% filter boundary vertices
[ed,~,c] = unique(ed,'rows');
% The last column lists the counts
cnt =  accumarray(c,1);

% bnd edges are those who only occur once
bndEdg = ed(cnt == 1,:);
bndVerts = unique(bndEdg);

if DEBUG
    figure('Name','boundary');
    plot3(fv.vertices(bndVerts,1),fv.vertices(bndVerts,2),fv.vertices(bndVerts,3),'ok');
    axis equal;
end


%% calcualate displacement direction axis
% the displacement axis is defined by the cross product between the vertex
% normal and the vector between the two neighbouring mesh boundary
% vertices (tangentional displacement)

% calculate vertex normals
normals = patchnormals(fv);

% calculate vector between the neighbouring vertices of each boundary
% vertex

tangent = zeros([length(bndVerts), 3]);
axisVec = zeros(size(tangent));

for v = 1:length(bndVerts)
    vert = bndVerts(v);
    %find edges containg this vertex
    isIn = ~~sum(bndEdg == vert,2);
    idx = unique(bndEdg(isIn,:));
    idx(idx == vert) = []; %indices of the vertices
    
    if length(idx) > 2
        error('invalid boundary');
    end
    
    tangent(v,:) = fv.vertices(idx(1),:) - fv.vertices(idx(2),:);
    
    axisVec(v,:) = cross(normals(v,:),tangent(v,:));
    
end


% scale to unit length
axisVec = axisVec./ repmat(sqrt(sum(axisVec.^2,2)),[1,3]);


if DEBUG
    hold on;
    quiver3(fv.vertices(bndVerts,1),fv.vertices(bndVerts,2),fv.vertices(bndVerts,3),...
        axisVec(:,1),axisVec(:,2),axisVec(:,3));
    axis equal;
end


%% calculate axis (ray) surface intersections
% coose the closest one as the new vertex position

vert1 = boundary.vertices(boundary.faces(:,1),:);
vert2 = boundary.vertices(boundary.faces(:,2),:);
vert3 = boundary.vertices(boundary.faces(:,3),:);

numBndTris = length(boundary.faces(:,1));

origins = fv.vertices(bndVerts,:);

for v = 1:length(bndVerts)
    
    
    [intersect, t] = TriangleRayIntersection(...
        repmat(origins(v,:),[numBndTris,1]) , repmat(axisVec(v,:),[numBndTris,1]),...
        vert1, vert2, vert3);
    
    
    shifts = t(intersect);
    
    shiftVec(v) = shifts(abs(shifts) == min(abs(shifts)));
    
end


% final shift
fv.vertices(bndVerts,:) = fv.vertices(bndVerts,:) + axisVec .* repmat(shiftVec', [1, 3]);

if DEBUG
    figure('Name','stretched mesh');
    patch(fv,'FaceColor',[0 1 1],'FaceAlpha',.6); axis equal; rotate3d on;
    hold on
    patch(boundary,'FaceColor',[1 1 0],'FaceAlpha',.3); axis equal; rotate3d on;
    
end

end
















