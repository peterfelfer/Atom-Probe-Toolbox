function deloc = patchDelocalizeProperty(fv,prop)
% averages the property prop between a vertex and its neighbours of a patch
% fv.


edgeList = tri2edgeList(fv.faces);
numVerts = length(fv.vertices);


if length(prop(1,:))>length(prop(:,1))
    prop = prop';
end

deloc = zeros(size(prop));


for vert = 1:numVerts
    neighbours = sum(edgeList == vert,2);
    neighbours = edgeList(logical(neighbours),:);
    neighbours = unique(neighbours);
    
    
    deloc(vert,:) = sum(prop(neighbours,:))/length(neighbours);
end