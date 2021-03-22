function [pass, Nmin, clusterCutoff, clusteredAtoms] = clusterDetermination(clusterPos,pos,Nmin)
%clusterDetermination calculates with the Voronoi tesselation and the
%Delaunay triangulation the cluster of the given dataset
%
% [pass, Nmin, clusterCutoff, clusteredAtoms] =
% clusterDetermination(clusterPos,pos,Nmin);
% [pass, Nmin, clusterCutoff, clusteredAtoms] =
% clusterDetermination(clusterPos,pos);
%
% INPUT
% clusterPos:   pos file that contains the atoms that are in the clusters
% pos:          pos file of the entire dataset
% Nmin:         at this value, the cluster has a higher probability that it
%               it is clustered non-randomly than it's probability to be a random cluster
%
% OUTPUT
% pass:         logical (1 or 0) 
%               0 the data set passes the Kolmogorov-Smirnov test, that
%               means that the voronoi volume distribution deviates from
%               random
%               1 the data set does not pass the Kolmogorov-Smirnov test
% Nmin:         at this value, the cluster has a higher probability that it
%               it is clustered non-randomly than it's probability to be a random cluster
% clusterCutoff:    max volume of the voronoi cell of the clustered atoms
% clusteredAtoms:   pos file of the clustered atoms
%
% For more information read the paper:
% Detecting and extracting clusters in atom probe data: A simple, automated
% method using Voronoi cells, P. Felfer et al, Ultramicroscopy 150 (2015)
% 30-36




%%  check if Nmin is a given value


if exist('Nmin','var')
    NminTmp = Nmin;
end




%% actual Voronoi cluster analysis
figure;
figName = [];
figName = ['Voronoi volume analysis of ' figName];

[numClustered, clusterCutoff, ~, experimentalVolumes, randomVolumes] = ...
    voronoiVolumeAnalysis(clusterPos, pos,true);

figure;
%% analysis of cluster sizes
% experimental
[clusterIdx, numClusters] = clusterIdentification(clusterPos,clusterCutoff,experimentalVolumes.expVol);

% random
[randomClusterIdx, randomNumClusters] = clusterIdentification(randomVolumes, clusterCutoff, randomVolumes.randVol);

clusterSizes = histcounts(clusterIdx, numClusters);
randomClusterSizes = histcounts(randomClusterIdx, randomNumClusters);
Nmin = clusterSizeAnalyse(clusterSizes, randomClusterSizes);

if exist('NminTmp','var')
    disp('Nmin override');
    Nmin = NminTmp;
end

%% Kolmogorov - Smirnov test:
significanceLimit = 1.92 / sqrt(height(clusterPos)) *100;

pass = (numClustered/height(clusterPos)) > significanceLimit;

clusterPct = (numClustered/height(clusterPos)) * 100;
% significanceLimit;

%% identifying clustered atoms
significantClusterIdx = find(clusterSizes >= Nmin);

% actually creating the atomic positions
isClustered = ismember(clusterIdx,significantClusterIdx);
%clusteredAtoms = [clusterPos(isClustered,:), clusterIdx(isClustered)'];


end

