function ptsFiltered = VoronoiFilter2D(pts,noiseLim)
%Voronoi filters a 2D hitmap, so that noise counts are removed. noiseLim is
%the fraction of atoms that is used to determine the limit for the Voronoi
%cell area cutoff. We only use the largest object.

DEBUG = true;


if ~exist('noiseLim','var')
    % default noise limit = 5%
    noiseLim = 0.05;
end


HUGE = 1E6;


% removing duplicate points
pts = unique(pts,'rows');


numPts = length(pts);
pts = double(pts);


%% calculating Voronoi areas
[V, C] = voronoin(pts);

A = zeros(numPts,1);

parfor c = 1:numPts
    if sum(C{c} == 1)
        A(c) = HUGE;
    else
        A(c) = polyarea(V(C{c},1),V(C{c},2));
    end
end


%% determine cutoff

maxIdx = round(numPts * (1-noiseLim));
Amax = sort(A);
Amax = Amax(maxIdx);

if DEBUG
    figure('Name', 'discarded points');
    triplot(delaunayn(pts),pts(:,1),pts(:,2),'g');
    hold on;
    scatter(pts(A > Amax,1),pts(A > Amax,2),'.k'); axis equal;
end

discard = A > Amax;


%% determine individual clusters from edge connectivity

% calcualte delaunay triangulation
tri = delaunayn(pts);

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
    scatter(pts(idx,1),pts(idx,2),'.r');
    
end

edg = edg(~edgDiscard,:);

edgSp = sparse(edg(:,1),edg(:,2),1);
edgSp(max(size(edgSp)),max(size(edgSp))) = 0;


[numComps compIdx] = graphconncomp(edgSp,'Directed',false);


% determine which component is the largest
compSz = tabulate(compIdx);

largest = compSz(compSz(:,2) == max(compSz(:,2)),1);

%% output 
ptsFiltered = pts(compIdx == largest,:);


if DEBUG
    figure('Name','adjacencyb graph');
    scatter(pts(:,1),pts(:,2),'.k'); axis equal
    hold on;
    scatter(ptsFiltered(:,1),ptsFiltered(:,2),'.r'); axis equal
end











