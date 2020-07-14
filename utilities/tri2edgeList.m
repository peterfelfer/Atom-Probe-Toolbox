function edgeList = tri2edgeList(tri,m)
%transforms a trinagulation of specified points into an edge list.
% also finds edges based on their manifold count m. m = 1: boundary edges, m >
% 2: non-manifold edges.

% creating all delaunay edges
permut = [1,2;2,3;3,1];

%building adjacency indices (= Delaunay edges) list of all atoms
delAdjacencyIdx = zeros(length(tri)*3,2);
numTet = length(tri);

for perm = 1:3
    %disp(['permutation ' num2str(perm)]);
    p = permut(perm,:);
    
    delAdjacencyIdx(numTet*(perm-1)+1:numTet*perm,:) = tri(:,p);
    
end

delAdjacencyIdx = [delAdjacencyIdx; [delAdjacencyIdx(:,2) delAdjacencyIdx(:,1)]];

edgeList = delAdjacencyIdx;

[edgeList idxa idxb] = unique(edgeList,'rows');



%% calculate edges of manifold# m
if exist('m','var')
    
    
    
    
    
    
end


end