function points = coordSystemTransform(ROIobject, points)
% transforms points from the global coordinate system to a local object defined 
% 
% points = coordSystemTransform(ROIobject, points)
% 
% INPUT
% ROIobject:    object which defines the region of interest  
%
% points:       points wihin the ROI shape
% 
% OUTPUT
% points:       points transformed to the new coordinate system    

temp = points;

ROIorigin = ROIobject.UserData.ROIxaxis(1,:);
ROIbasis = [ROIobject.UserData.ROIxaxis(2,:) - ROIobject.UserData.ROIxaxis(1,:);...
    ROIobject.UserData.ROIyaxis(2,:) - ROIobject.UserData.ROIyaxis(1,:);...
    ROIobject.UserData.ROIzaxis(2,:) - ROIobject.UserData.ROIzaxis(1,:)];
    
temp = (global2localcoord([temp.x, temp.y, temp.z]','rr',ROIorigin',ROIbasis'))';

points(:,2:4) = array2table(temp);

