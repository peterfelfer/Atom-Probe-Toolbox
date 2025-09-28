function ptsFiltered = voronoiFilter2D(pointCoords,noiseLimit)
%Voronoi filters a 2D hitmap, so that noise counts are removed. noiseLim is
%the fraction of atoms that is used to determine the limit for the Voronoi
%cell area cutoff. We only use the largest object.
% The scan hull is an alpha Hull that is produced for an atom probe
% detector stack. In order to simplify the topology, it is calculated by
% slicing the dataset into numSeg. 
%
% [hull, Az] = detectorHull(detectorData,alpha,numSeg, nGon, capLoops);
% 
% INPUT
% pointCoords:  2xN variable of point coordinates, e.g. [detx dety]
%
% noiseLimit:   fraction of hits with the larges voronoi areas that are
%               removed for detector hull calculation
%
% OUTPUT
% ptsFiltered:  filtered points with noisLimit fraction of points with the largest Voronoi areas removed
%
% WARNING: using the voronoi filter takes considerable time on larger
% datasets. get lunch.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

DEBUG = false;
if ~exist('noiseLim','var')
    % default noise limit = 5%
    noiseLimit = 0.05;
end
HUGE = 1E6;
% removing duplicate points
pointCoords = unique(pointCoords,'rows');

numPts = length(pointCoords);
pointCoords = double(pointCoords);

%% calculating Voronoi areas
[V, C] = voronoin(pointCoords);

A = zeros(numPts,1);

parfor c = 1:numPts
    if sum(C{c} == 1)
        A(c) = HUGE;
    else
        A(c) = polyarea(V(C{c},1),V(C{c},2));
    end
end

%% determine cutoff

maxIdx = round(numPts * (1-noiseLimit));
Amax = sort(A);
Amax = Amax(maxIdx);

if DEBUG
    figure('Name', 'discarded points');
    triplot(delaunayn(pointCoords),pointCoords(:,1),pointCoords(:,2),'g');
    hold on;
    scatter(pointCoords(A > Amax,1),pointCoords(A > Amax,2),'.k'); axis equal;
end

discard = A > Amax;


%% determine individual clusters from edge connectivity

% calcualte delaunay triangulation
tri = delaunayn(pointCoords);

edg = [tri(:,[1 2]); tri(:, [2 3]); tri(:, [1 3])];


% delete edges that contain points that contain points to be discarded

edgDiscard = ismember(edg, find(discard));
edgDiscard = ~~sum(edgDiscard,2);

sum(edgDiscard);



% draw discarded edges
if false
    %figure('Name', 'endpoints of discarded edges');
    hold on
    idx = unique(edg(edgDiscard,:));
    scatter(pointCoords(idx,1),pointCoords(idx,2),'.r');
    
end

edg = edg(~edgDiscard,:);

edgSp = sparse(edg(:,1),edg(:,2),1);
edgSp(max(size(edgSp)),max(size(edgSp))) = 0;


[~, compIdx] = graphconncomp(edgSp,'Directed',false);


% determine which component is the largest
compSz = tabulate(compIdx);

largest = compSz(compSz(:,2) == max(compSz(:,2)),1);

%% output 
ptsFiltered = pointCoords(compIdx == largest,:);


if DEBUG
    figure('Name','adjacencyb graph');
    scatter(pointCoords(:,1),pointCoords(:,2),'.k'); axis equal
    hold on;
    scatter(ptsFiltered(:,1),ptsFiltered(:,2),'.r'); axis equal
end











