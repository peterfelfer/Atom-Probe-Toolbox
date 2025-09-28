function [hull, Az] = detectorHull(detectorData,alpha,numSeg, nGon, capLoops, useVoronoiFilter, filterNoiseLimit)
% The scan hull is an alpha Hull that is produced for an atom probe
% detector stack. In order to simplify the topology, it is calculated by
% slicing the dataset into numSeg. 
%
% [hull, Az] = detectorHull(detectorData,alpha,numSeg, nGon, capLoops);
% [hull, Az] = detectorHull(detectorData,alpha,numSeg, nGon, capLoops, useVoronoiFilter, filterNoiseLimit);
% 
% INPUT
% detectorData: pos table with detector entries detx and dety or 2xN
%               variable [detx dety] or epos format aray 11xN
%
% alpha:        alpha value for the individual segments
%
% numSeg:       number of segments the dataset is divided into
%
% nGon:         number of vertices of each segment ring
%
% capLoops:     number of edge loops in top and bottom cap
%
% useVoronoiFilter: Using of a Voroin filter to remove detector noise hits
%               if the detector is not fully filled.
%
% filterNoiseLimit: fraction of hits with the larges voronoi areas that are
%               removed for detector hull calculation
%
% OUTPUT
% hull:         array of patches (hull(i).vertices, hull(i).faces)
%               representing top cap, mantle and bottom cap of the detector hull
%
% Az:           numSeg x 2 array with average hit index of slice and effective 
%               detector area of the slice.
%
% WARNING: using the voronoi filter takes considerable time on larger
% datasets. get lunch.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

DEBUG = false;
vis = true;

%noise limit of Voronoi filtering of detector hits
if not(exist('useVoronoiFilter','var'))
    VORONOIFILTER = false;
else
    VORONOIFILTER = useVoronoiFilter;
end

if not(exist("filterNoiseLimit",'var'))
    NOISELIM = 0.10;
else
    NOISELIM = filterNoiseLimit;
end

ZMAX = 100;
RMAX = 100;

if istable(detectorData)
    detectorData = [detectorData.detx, detectorData.dety];
elseif length(detectorData(1,:))>2
    detectorData = detectorData(:,8:9);
end


if ~exist('alpha','var')
    alpha = 3;
end

if ~exist('numSeg','var')
    numSeg = 10;
end

if ~exist('nGon','var')
    nGon = 64;
end

if ~exist('capLoops','var')
    capLoops = 4;
end

if ~exist('incMode','var')
    incMode = 'linear';
end


numAtoms = length(detectorData(:,1));
atPerSlice = numAtoms/numSeg;
addpath('utilities_geom2d','utilities_geom2d/polygons2d','utilities_geom2d/geom2d');
detectorData = double(detectorData);

if DEBUG
    numAtoms = length(detectorData)
    
    alpha
    numSeg
    nGon
    capLoops
end


%% calculating alpha shapes of detector coordinates
wb = waitbar(0,'calculating alpha shapes');
hull = {};

for slice = 1:numSeg
    if DEBUG
        disp(['segment ' num2str(slice) ' of ' num2str(numSeg)]);
    end
    
    

    idxBeg = uint32(floor(atPerSlice*(slice-1)))+1;
    idxEnd = uint32(floor(atPerSlice*slice))-1;
    idxAvg = round((idxBeg+idxEnd)/2);
    
    
    sliceX = detectorData(idxBeg : idxEnd,1:2);
    
    % using Voronoi filter to get rid of noise
    if VORONOIFILTER
        sliceX = VoronoiFilter2D(sliceX,NOISELIM);
    end
    
    [area, alphaShape] = alphavol(sliceX,alpha);
    Az(slice,:) = [idxAvg area];
    %   [V,S] = ALPHAVOL(X,R) outputs a structure S with fields:
    %    S.tri - Triangulation of the alpha shape (Mx3 or Mx4)
    %    S.vol - Area or volume of simplices in triangulation (Mx1)
    %    S.rcc - Circumradius of simplices in triangulation (Mx1)
    %    S.bnd - Boundary facets (Px2 or Px3)
    loop = [sliceX(alphaShape.bnd(:,1),1),sliceX(alphaShape.bnd(:,2),2)];
    loop = [loop; loop(1,:)];
    
    % calculating segmentation (n - gon)
    inc = 2*pi()/nGon;
    
    outline = zeros(nGon,2);
    for n = 1:nGon
        [x,y] = pol2cart(inc*n,RMAX);
        edge = [0 0 x y];
        try
            tmp = intersectEdgePolygon(edge,loop);
            outline(n,:) = tmp(1,:);
            
            
        catch
            
            warning(['error in slice ' num2str(slice) ', segment ' num2str(n)]);
        end
        
    end
    
    
    hull{slice} = [outline, double(ones(length(outline(:,1)),1))*double(idxAvg)];
    
    waitbar(slice/numSeg,wb,['segment ' num2str(slice) ' of ' num2str(numSeg)]);
end
close(wb);


% adding cap edge (upper and lower edge of mantle)
upperEdge = hull{1};
upperEdge(:,3) = ones(length(upperEdge(:,1)),1);

lowerEdge = hull{end};
lowerEdge(:,3) = ones(length(upperEdge(:,1)),1)*numAtoms;

hull = [{upperEdge} hull {lowerEdge}];



%% adding mantle triangulation (separate objects mantle and caps)

% the triangles are foremd by connecting each point with index n to the point of the
% same index n in the next slice s+1 and n+1 in the same slice s and n-1 in the
% next slice s+1 (two triangles: [(n,s)(n,s+1)(n+1,s) |
% (n,s)(n,s+1)(n-1,s+1)] for slices s = 1: second last

mantle.vertices = [];
mantle.faces = [];

for s = 1:length(hull)
    mantle.vertices = [mantle.vertices; hull{s}];
end


for s = 0:length(hull)-2
    for n = 1:nGon
        % create triangles
        mantle.faces = [mantle.faces; n+s*nGon n+1+s*nGon n+(s+1)*nGon; ...
            n+s*nGon  n+(s+1)*nGon n-1+(s+1)*nGon];
        
    end
end

if false
    figure('Name','mantle');
    
    mantleTMP = mantle;
    mantleTMP.vertices(:,3) = mantleTMP.vertices(:,3)/mantleTMP.vertices(end,3) * 100;
    
    patch(mantleTMP,'FaceColor',[0 1 1],'FaceAlpha',.6); axis equal; rotate3d on;
end


%% adding caps
% top cap
rim = hull{1}(:,1:2);
topCap.vertices = [];

%vertex locations
for pt = 1:length(rim(:,1))
    [theta rho] = cart2pol(rim(pt,1),rim(pt,2));
    
    
    for lp = 1:capLoops
        [x,y] = pol2cart(theta,rho/capLoops * lp);
        topCap.vertices(end+1,:) = [x y];
    end
    
end

% triangulation and 3d ifying
topCap.vertices(end+1,:) = [0,0];
numVerts = length(topCap.vertices);

%% top and bottom cap face triangulation
tris = [];
for pt = 1:(nGon-1)
    tris(end+1,:) = [numVerts, 1 + (pt - 1) * capLoops, 1 + pt * capLoops];
    
    for lp = 1:(capLoops - 1)
        tris(end+1,:) = [1+(pt-1)*capLoops+(lp-1), 1+(pt-1)*capLoops+lp, 1+pt*capLoops+(lp-1)];
        tris(end+1,:) = [1+(pt-1)*capLoops+lp, 1+pt*capLoops+lp, 1+pt*capLoops+(lp-1)];
    end
end
% last segment
tris(end+1,:) = [numVerts, 1 + (nGon - 1) * capLoops, 1 + 0 * capLoops];
for lp = 1:(capLoops-1)
    tris(end+1,:) = [1+(nGon-1)*capLoops+(lp-1), 1+(nGon-1)*capLoops+lp, 1+0*capLoops+(lp-1)];
    tris(end+1,:) = [1+(nGon-1)*capLoops+lp, 1+0*capLoops+lp, 1+0*capLoops+(lp-1)];
end


%% top cap
topCap.faces = tris;
topCap.vertices(:,3) = ones(length(topCap.vertices(:,1)),1);








%% bottom cap
rim = hull{end}(:,1:2);
bottomCap.vertices = [];

for pt = 1:length(rim(:,1))
    [theta rho] = cart2pol(rim(pt,1),rim(pt,2));
    for lp = 1:capLoops
        [x,y] = pol2cart(theta,rho/capLoops * lp);
        bottomCap.vertices(end+1,:) = [x y];
    end
    
end

% triangulation and 3d ifying
bottomCap.vertices(end+1,:) = [0,0];
bottomCap.faces = tris;
bottomCap.vertices(:,3) = ones(length(topCap.vertices(:,1)),1)*numAtoms;




%% visulaisation
if vis
    %% drawing the shape
    f = figure('Name','detector hull');
    p = patch(mantle);
    set(p,'FaceColor','g')
    hold on
    tp = patch(topCap);
    set(tp,'FaceColor','r');
    hold on
    bp = patch(bottomCap);
    set(bp,'FaceColor','b');
    
    axis square;
    rotate3d on;
end

clear hull;
hull(1) = mantle;
hull(2) = topCap;
hull(3) = bottomCap;

end






















