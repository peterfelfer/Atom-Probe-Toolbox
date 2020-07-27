function points = coordSystemTransform(ROIobject, points)
% coordSystemTransform transforms points from the global coordinate system to a local object defined 
% 
% points = coordSystemTransform(ROIobject, points)
% 
% INPUT
% ROIobject:    object which defines the region of interest  
%
% points:       points to be transformed to the ROI or WCS coordinate system
%               can be a Nx3 array of coordinates as [x y z] or a table
%               with points.x points.y and points.z as columns to be
%               transfomed
% 
% OUTPUT
% points:       points transformed to the ROI or WCS coordinate system.
%               Adheres to the format of the input, i.e. an Nx3 array input
%               produces and Nx3 output, while a table input produces a
%               table output in which only the x, y and z columns have been
%               changed.


isTab = istable(points);

if isTab
    x = points.x;
    y = points.y;
    z = points.z;
else
    x = points(:,1);
    y = points(:,2);
    z = points(:,3);
end

% coordinate origin of the ROI
ROIorigin = ROIobject.UserData.ROIxaxis(1,:); 

% ROI coordinate system basis orientation
ROIbasis = [ROIobject.UserData.ROIxaxis(2,:) - ROIobject.UserData.ROIxaxis(1,:);...
    ROIobject.UserData.ROIyaxis(2,:) - ROIobject.UserData.ROIyaxis(1,:);...
    ROIobject.UserData.ROIzaxis(2,:) - ROIobject.UserData.ROIzaxis(1,:)];
    
pointsTrans = (global2localcoord([x, y, z]','rr',ROIorigin',ROIbasis'))';

if isTab
    points.x = pointsTrans(:,1);
    points.y = pointsTrans(:,2);
    points.z = pointsTrans(:,3);
else
    points = pointsTrans;
end

