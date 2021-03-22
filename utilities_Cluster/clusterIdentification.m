function [clusterIdx, numClusters] = clusterIdentification(clusterPos,volThreshORpos,volVoronoi,posDelaunay)
% clusterIdentification calculates the number of clusters and their
% corresponding cluster Index after Voronoi Tesselation and Dealaunay
% triangulation
% 
% [clusterIdx, numClusters] =
% clusterIdentification(clusterPos,volThreshORpos,volVoronoi,posDelaunay);
% [clusterIdx, numClusters] =
% clusterIdentification(clusterPos,volThreshORpos,volVoronoi);
% [clusterIdx, numClusters] =
% clusterIdentification(clusterPos,volThreshORpos);
%
% INPUT
% clusterPos:       a pos file with the atoms that are in the cluster 
% volThreshORpos:   is the clusterCutoff calculated by the voronoiVolume
%                   Analysis or a pos file containing the entire data set;
% volVoronoi:       the voronoi volume of the atoms from clusterPos
% posDelaunay:      delaunay triangulation of the clusterPos data set, is a mx4
%                   matrix that defines the tetrahedron
%
% OUTPUT
% clusterIdx:       is a vector with the corresponding cluster number for
%                   each clustered atom
% numClusters:      total number of clusters


%% check for inputs

if ~exist('volVoronoi','var')
    volVoronoi = vertexVolume(clusterPos);
end
if~exist('posDelaunay','var')
    posDelaunay = delaunay([clusterPos.x clusterPos.y clusterPos.z]);
end
% patchDelaunay = tetramesh(posDelaunay,[clusterPos.x clusterPos.y
% clusterPos.z]); % visualisation of the triangles

%% checking if the volume threshold is defined. If not, it is calculated
if height(volThreshORpos) == 1
    volThresh = volThreshORpos;
else
    [~, volThresh] = voronoiVolumeAnalysis(clusterPos,volThreshORpos);
end

clear volThreshORpos;




%% defining the clustered atoms and calculating the Delaunay triangulation
clusteredAtomsIndices = find(volVoronoi<=volThresh);


%% determining in the delaunay triangulation which atoms are clustered and to which cluster they belong
delAdjacency = delaunayToAdjacencyMat(posDelaunay,clusteredAtomsIndices);

% Eine Adjazenzmatrix (manchmal auch Nachbarschaftsmatrix) eines Graphen ist
% eine Matrix, die speichert, welche Knoten des Graphen durch eine Kante verbunden sind. 
% Sie besitzt für jeden Knoten eine Zeile und eine Spalte, woraus sich für n Knoten 
% eine n x n-Matrix ergibt. Ein Eintrag in der i-ten 
% Zeile und j-ten Spalte gibt hierbei an, ob eine Kante von dem i-ten zu dem j-ten Knoten führt. 
% Steht an dieser Stelle eine 0, ist keine Kante vorhanden – eine 1 gibt an, 
% dass eine Kante existiert[1], siehe Abbildung rechts. (wikipedia)

%[numClusters clusterIdx] = graphconncomp(delAdjacency,'Directed',false);
[numClusters, clusterIdx] = conncomp(delAdjacency);



end