function lh = drawROIaxes(roi)
% drawROIaxes draws lines denoting the analysis axes for the ROI object 
% red = x, green = y, blue = z
%
% lh = drawROiaxes(roi)
%
% INPUT
% roi: name of roi as given before
%
% OUTPUT
% lh:   handle to the axes of the roi, colors as described before 
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

lh(1) = line(roi.UserData.ROIxaxis(:,1),roi.UserData.ROIxaxis(:,2),roi.UserData.ROIxaxis(:,3),'LineWidth',4,'Color','r');
lh(2) = line(roi.UserData.ROIyaxis(:,1),roi.UserData.ROIyaxis(:,2),roi.UserData.ROIyaxis(:,3),'LineWidth',4,'Color','g');
lh(3) = line(roi.UserData.ROIzaxis(:,1),roi.UserData.ROIzaxis(:,2),roi.UserData.ROIzaxis(:,3),'LineWidth',4,'Color','b');