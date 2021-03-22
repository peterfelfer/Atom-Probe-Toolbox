function adjacencyMatrix = delaunayToAdjacencyMat(delTet,qualifiedPointsIdx)
% delaunayToAdjacencyMat transforms a Delaunay trinagulation of specified points into a sparse
% adjacency matrix determining which atoms are clustered after the
% Voronoitesselation and Delaunay triangulation
% 
% adjacencyMatrix = delaunayToAdjacencyMat(delTet,qualifiedPointsIdx)
% 
% INPUT
% delTet:     delaunay triangulation of the clusterPos data set, is a mx4
%             matrix that defines the tetrahedron
% qualifiedPointsIdx:   indices of the points that are clustered after the
%             voronoi tesselation step. All of the points in the list have a Voronoi
%             Volume that is smaller than the volume threshold of the Voronoi
%             tesselation
%
% OUTPUT
% adjacencyMatrix = is a sparse adjacency matrix(Adjazenmatrix)
%
%
% delIsClustered = arrayfun(@(x) sum(x == qualifiedPointsIdx),delTet);
% delIsClustered = ismember(delTet,qualifiedPointsIdx);

% disp('calculating perms');
% creating all delaunay edges
permut = [1,2;2,3;3,4;1,3;1,4;2,4];

%building adjacency indices (= Delaunay edges) list of all atoms
delAdjacencyIdx = zeros(length(delTet)*6,2);
numTet = length(delTet);

for perm = 1:6
    %disp(['permutation ' num2str(perm)]);
    p = permut(perm,:);
    
    delAdjacencyIdx(numTet*(perm-1)+1:numTet*perm,:) = delTet(:,p);
    
end
% delAdjacencyIdx has for all of the tetrahedrons the adjacent
% points
delAdjacencyIdx = [delAdjacencyIdx; [delAdjacencyIdx(:,2) delAdjacencyIdx(:,1)]];

%% removing the edges which contain atoms that are not qualified

% check if the adjacent points are points of the qualifiedPointsIdx, that
% means that their volume is smaller than the voronoi Threshold
%disp('calculating memebership of edges');
%delAdjacencyIdxIsIn = arrayfun(@(x) sum(x == qualifiedPointsIdx),delAdjacencyIdx);
delAdjacencyIdxIsIn = ismember(delAdjacencyIdx,qualifiedPointsIdx);
%size(delAdjacencyIdxIsIn)
%sum(delAdjacencyIdxIsIn)

unqualifiedEdges = ~(delAdjacencyIdxIsIn(:,1) .* delAdjacencyIdxIsIn(:,2));

%sum(unqualifiedEdges)

%length(delAdjacencyIdx)
delAdjacencyIdx(unqualifiedEdges,:) = [];



adjacencyMatrix = sparse(delAdjacencyIdx(:,1),delAdjacencyIdx(:,2),1);


end