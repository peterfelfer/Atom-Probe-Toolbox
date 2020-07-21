function points = coordSystemTransform(ROIobject, points)
%transforms points from the global coordinate system to a local object defined oe2jhdkvkueqwzfvciuzevwbufihscjkbEJKQWBNC%The basis vectors are stored in the ROI object 'o' as: o.UserData.ROIxaxis
%The basis vectors contain two components: Basis starting and end vector

%If points contains more than 3 colums, only the first three will be treated as coordinates
%(e.g. for pos / epos). Basis vectors need not be of unit length.

temp = points;

ROIorigin = ROIobject.UserData.ROIxaxis(1,:);
ROIbasis = [ROIobject.UserData.ROIxaxis(2,:) - ROIobject.UserData.ROIxaxis(1,:);...
    ROIobject.UserData.ROIyaxis(2,:) - ROIobject.UserData.ROIyaxis(1,:);...
    ROIobject.UserData.ROIzaxis(2,:) - ROIobject.UserData.ROIzaxis(1,:)];
    
temp = (global2localcoord([temp.x, temp.y, temp.z]','rr',ROIorigin',ROIbasis'))';

points(:,2:4) = array2table(temp);

