function [vol, errCts] = vertexVolume(pointCloud, serial)
% vertexVolume calculates the volume a point takes up in a point cloud by voronoi
% tesselating the space within the point cloud and calculating the volume
% of the convex hull of every point. Boundary points need special
% treatment: for the volume, a sphere with the radius of the median
% distance point of the voronoi vertices to the actual point is assumed
%
% [vol errCts] = vertexVolume(coords,serial)
% 
% INPUT
% pointCloud: table with the x y and z coordinates of the points in the point cloud 
% serial: true, the calculation is faster and multiple cores are used
% 
% OUTPUT
% vol:      the Voronoi volume for each point in the pointclod
% errCts:   errorCounts 
%
%% was sind genau errCts? brauchen wir das ParforProgress2 überhaupt?



%% Qhull needs double precision
% pointCloudUnique = unique(pointCloud);
% coords = double([pointCloudUnique.x pointCloudUnique.y pointCloudUnique.z]);
pointCloud = unique(pointCloud);
coords = double([pointCloud.x pointCloud.y pointCloud.z]);
%% multicore version
if ~exist('serial','var')
    %addpath('ParforProgress2');
    
    
    numCoords = length(coords(:,1));
    
    
    [Vverts, Vcells] = voronoin(coords(:,1:3));
    
    %ppm = ParforProgressStarter2('calculating voronoi volumes',numCoords,1);
    for vvol = 1:numCoords
        cellVertices = Vverts(Vcells{vvol},:);
        if sum(isinf(cellVertices))
            
            cellVertices = cellVertices(isfinite(cellVertices(:,1)),:);
            cellVertices = cellVertices - repmat(coords(vvol,1:3),length(cellVertices(:,1)),1);
            cellDist = sum(cellVertices.^2,2);
            cellDist = sqrt(median(cellDist));
            
            vol(vvol) = 4/3*pi()*cellDist^3;
        else
            cellVertices = cellVertices(isfinite(cellVertices(:,1)),:);
            try
                [~,vol(vvol)] = convhulln(cellVertices);
            catch
                % if the points are coplanar within the accuracy of
                % conversion from single precision pos file, the volume
                % they occupy is 0 by definition.
                vol(vvol) = 0;
            end
        end
        
        %vol(vvol) = v;
        
        %ppm.increment(vvol);
    end
  
    
%delete(ppm);
else
%% serial version (for reference)        

    
    
    errCts = 0;
    
    numCoords = length(coords(:,1));
    
    [Vverts, Vcells] = voronoin(coords(:,1:3));
    
    
    h = waitbar(0,'calculating voronoi volumes');
    for vvol = 1:numCoords
        
        cellVertices = Vverts(Vcells{vvol},:);
        if sum(isinf(cellVertices))
            
            cellVertices = cellVertices(isfinite(cellVertices(:,1)),:);
            cellVertices = cellVertices - repmat(coords(vvol,1:3),length(cellVertices(:,1)),1);
            cellDist = sum(cellVertices.^2,2);
            cellDist = sqrt(median(cellDist));
            
            v = 4/3*pi()*cellDist^3;
        else
            cellVertices = cellVertices(isfinite(cellVertices(:,1)),:);
            try
                [~,v] = convhulln(cellVertices);
            catch
                v = 0;
                errCts = errCts+1;
            end
        end
        
        vol(vvol) = v;
        
        waitbar(vvol/numCoords);
        
    end
    
    delete(h);
end



end