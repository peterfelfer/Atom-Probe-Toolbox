function remesh = patchCentroidalVoronoiTessellation(origMesh, CVTsteps, angLim)
% remeshes a patch origMesh by sampling so that the Voronoi cell of each point on
% avg contains avgCount points and performs a relaxation using Lloyds
% algorithm. boundary point where the edge angle is  < angLim are not
% influneced. If CVTsteps < 1, an energy criterion for convergence is used
% and the amount of iterations is automatically determined.

DEBUG = false;
PLOT = false;


%addpath('PointTool','toolbox_graph','fitNormal','patch_normals','resources');

if DEBUG
    disp('debugging mode');
    ph = patch(origMesh,'CData',[0 1 0],'FaceAlpha',0.6,'FaceColor',[0 1 0],'EdgeColor',[0 0 0]);
    ah = gca;
    axis equal
    rotate3d on
    hold on
end




remesh.vertices = origMesh.vertices;
remesh.faces = origMesh.faces;

numVerts = length(remesh.vertices(:,1));

%if DEBUG
%    scatter3(remesh.vertices(:,1),remesh.vertices(:,2),remesh.vertices(:,3),'MarkerFaceColor','r','MarkerEdgeColor','k','SizeData',72);
%end


% recalculating triangulation




%% detect boundary and resample boundary
edgeList = tri2edgeList(remesh.faces);


% detect edge
bndVerts = computeBoundary(remesh);



%% calcualte CVT (Lloyd's algorithm) and perform relaxation (separate for all boundary points)

if DEBUG
    beg = 2844;%randsample(60,1)
    numVerts = beg;
    
    CVTsteps = 1
else
    beg = 1;
end




for i = 1:CVTsteps
    
    normals = patchnormals(remesh);
    
    wb = waitbar(0,['performing CVT step ' num2str(i)]);
    

    for vert = beg:numVerts
        if DEBUG
            vert
            is_boundary = ismember(vert,bndVerts)
        end
        
        
        if ismember(vert,bndVerts)
            
            % get neighbouring vertices
            currentVert = remesh.vertices(vert,:);
            
            neighbours = sum(edgeList == vert,2);
            neighbours = edgeList(logical(neighbours),:);
            neighbours = unique(neighbours);
            neighbours = neighbours(neighbours~=vert);
            
            for n = 1:length(neighbours)
                neighbours(n) = neighbours(n)*ismember(neighbours(n),bndVerts);
                
            end
            neighbours(neighbours == 0) = [];
            
            
            
            neighbourVerts = remesh.vertices(neighbours,:);
            
            if DEBUG
                scatter3(ah,currentVert(:,1),currentVert(:,2),currentVert(:,3),'r','SizeData',128,'MarkerFaceColor','r');
                scatter3(ah,neighbourVerts(:,1),neighbourVerts(:,2),neighbourVerts(:,3),'b','SizeData',128,'MarkerFaceColor','b');
                
            end
            

            vec1 = currentVert - neighbourVerts(1,:);
            vec2 = neighbourVerts(2,:) - currentVert;
            
            d1 = norm(vec1);
            d2 = norm(vec2);
            
            mid = (d1+d2)/2;
            
            % the new point must lie on the longer vector
            
            if d1>d2
                shift3d = vec1/d1 * mid;
                
                cent3 = neighbourVerts(1,:) + shift3d;
            else
                shift3d = vec2/d2 * (mid - d1);
                
                cent3 = currentVert + shift3d;
            end
            
            
            
            ang = radtodeg(acos(dot(vec1,vec2)/(d1+d2)));
            
            
            %plot all 3d properties
            if DEBUG
                scatter3(ah,cent3(:,1),cent3(:,2),cent3(:,3),'k','SizeData',128,'MarkerFaceColor','k');
                
                disp(['edge angle: ' num2str(ang)]);
                
            end
            
            if ang < angLim
                remesh.vertices(vert,:) = cent3;
            end
            
            
            
            
            
            
        else % do not relax boundary vertices
            
            % get polygon around that vertex (all edges that contain the point)
            currentVert = remesh.vertices(vert,:);
            
            neighbours = sum(edgeList == vert,2);
            neighbours = edgeList(logical(neighbours),:);
            neighbours = unique(neighbours);
            neighbours = neighbours(neighbours~=vert);
            
            neighbourVerts = remesh.vertices(neighbours,:);
            
            %poly = remesh.vertices(neighbours,:);
            
            % project polygon into 2D space
            normal = normals(vert,:);
            
            if DEBUG
                
                scatter3(ah,currentVert(:,1),currentVert(:,2),currentVert(:,3),'r','SizeData',128,'MarkerFaceColor','r');
                scatter3(ah,neighbourVerts(:,1),neighbourVerts(:,2),neighbourVerts(:,3),'b','SizeData',128,'MarkerFaceColor','b');
                
            end
            
            R = null(normal);
            
            uv = neighbourVerts * R;
            
            % to be fixed (hull is not generally convex)
            try
                edg = convhull(uv);
            catch
                error(['bad neighbourhood in vertex ' num2str(vert)]);
            end


            
            
            neighbourVerts = neighbourVerts(edg,:); %need to resort too, otherwise the order is screwed up when we change to 3d space again.
            uv = uv(edg,:);
            
            orig = currentVert * R;
            cent = polygonCentroid(uv);
            shift = cent - orig;
            
            
            
            
            % translate the shift to 3D, angle direction needs to be determined
            % since R is not unique (nullspace, singular value decomposition)
            
            %find intersection of shift vector with polygon
            
            HUGE = 1E6;
            
            line = [orig; orig + shift * HUGE];
            
            [xc yc ii] = polyxpoly(uv(:,1),uv(:,2),line(:,1),line(:,2));
            
            ii = ii(1);
            
            % for convex boundary pointd there are two intersections, we
            % entatively take the first one (only for test purposes. need to
            % treat boarder seperately
            if length(xc) > 1
                xc = xc(1);
                yc = yc(1);
                ii = ii(1);
            end
            
            
            % find the distances to the 2d vertices involved
            
            d1 = norm(uv(ii,:) - [xc yc]);
            d2 = norm(uv(ii+1,:) - [xc yc]);
            
            dRatio = d1 / (d1+d2);
            
            
            
            % plot all 2d properties
            if DEBUG
                %figure();
                %voronoi([uv(:,1); orig(1)],[uv(:,2); orig(2)]);
                
                figure();
                plot(uv(:,1),uv(:,2),'-o');
                hold on
                axis equal
                
                scatter(uv(1,1),uv(1,2),'.k');
                scatter(uv(end-1,1),uv(end-1,2),'.g');
                scatter(xc,yc,'.r');
                scatter(cent(1),cent(2),'g');
                scatter(orig(1),orig(2),'r');
                quiver(orig(1),orig(2),shift(1),shift(2));
                
                
                
            end
            
            
            targetPoint = neighbourVerts(ii,:) + (neighbourVerts(ii+1,:) - neighbourVerts(ii,:)) * dRatio;
            
            shift3d = norm(shift) * (targetPoint - currentVert) / norm(targetPoint - currentVert);
            
            cent3 = currentVert + shift3d;
            
            
            %plot all 3d properties
            if DEBUG
                %disp(['angle = ' num2str(rad2deg(ang),3) ' deg']);
                scatter3(ah,neighbourVerts(ii,1),neighbourVerts(ii,2),neighbourVerts(ii,3),'or','SizeData',256);
                scatter3(ah,neighbourVerts(ii+1,1),neighbourVerts(ii+1,2),neighbourVerts(ii+1,3),'ok','SizeData',256');
                
                scatter3(ah,targetPoint(1),targetPoint(2),targetPoint(3),'k','SizeData',128,'MarkerFaceColor','k');
                scatter3(ah,cent3(:,1),cent3(:,2),cent3(:,3),'k','SizeData',128,'MarkerFaceColor','k');
                %quiver3(ah,currentVert(:,1),currentVert(:,2),currentVert(:,3),shift3d(1),shift3d(2),shift3d(3));
            end
            
            
            remesh.vertices(vert,:) = cent3;
            
        end
        
        waitbar(vert/numVerts,wb);
        
    end
    
    close(wb);
end

if PLOT
    figure('Name','before');
    ph = patch(origMesh,'CData',[0 1 0],'FaceAlpha',0.8,'FaceColor',[0 1 1],'EdgeColor',[0 0 0]);
    axis equal
    rotate3d on
    
    
    figure('Name',['after ' num2str(CVTsteps) ' CVT steps']);
    ph = patch(remesh,'CData',[0 1 0],'FaceAlpha',0.8,'FaceColor',[0 1 0],'EdgeColor',[0 0 0]);
    axis equal
    rotate3d on
end





end





function bnd = computeBoundary(fv)
% computes boundary vertices of a mesh
% boundary vertices are all vertices that only belong to edges that only
% belong to one triangle

DEBUG = false;

%% creating edge list
edgList = [fv.faces(:,[1 2]); fv.faces(:, [1 3]); fv.faces(:, [2 3])];
edgList = sort(edgList,2);



%% counting which ones occur once
[edgListF,b,c] = unique(edgList,'rows');

% The last column lists the counts
cnt =  accumarray(c,1);

BNDedg = edgListF(cnt == 1,:);

bnd = unique(BNDedg);


if DEBUG
    figure('Name','boundary of mesh');
    scatter3(fv.vertices(bnd,1),fv.vertices(bnd,2),fv.vertices(bnd,3));
    rotate3d on; axis equal;
end


end




function [centroid area] = polygonCentroid(varargin)
%POLYGONCENTROID Compute the centroid (center of mass) of a polygon
%
%   CENTROID = polygonCentroid(POLY)
%   CENTROID = polygonCentroid(PTX, PTY)
%   Computes center of mass of a polygon defined by POLY. POLY is a N-by-2
%   array of double containing coordinates of vertices.
%
%   [CENTROID AREA] = polygonCentroid(POLY)
%   Also returns the (signed) area of the polygon. 
%
%   See also:
%   polygons2d, polygonArea, drawPolygon
%
%   ---------
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 05/05/2004.
%

% HISTORY
% 2012.02.24 vectorize code

if nargin==1
    var = varargin{1};
    px = var(:,1);
    py = var(:,2);
elseif nargin==2
    px = varargin{1};
    py = varargin{2};
end

% Algorithme P. Bourke, vectorized version
N = length(px);
iNext = [2:N 1];
common = (px .* py(iNext) - px(iNext) .* py);
sx = sum((px + px(iNext)) .* common);
sy = sum((py + py(iNext)) .* common);

area = sum(common) / 2;
centroid = [sx sy] / 6 / area;

end
























