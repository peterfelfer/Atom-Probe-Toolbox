function lh = drawROIaxes(roi)
% draws lines denoting the analysis axes for the ROI object red = x, 
%green = y, blue = z
%
% bh = drawROiaxes(roi)
%
% INPUTS
% roi: name of roi as given before
%
% OUTPUTS
% bh:   handle to the axes of the roi, colors as described before        
lh(1) = line(roi.UserData.ROIxaxis(:,1),roi.UserData.ROIxaxis(:,2),roi.UserData.ROIxaxis(:,3),'LineWidth',4,'Color','r');
lh(2) = line(roi.UserData.ROIyaxis(:,1),roi.UserData.ROIyaxis(:,2),roi.UserData.ROIyaxis(:,3),'LineWidth',4,'Color','g');
lh(3) = line(roi.UserData.ROIzaxis(:,1),roi.UserData.ROIzaxis(:,2),roi.UserData.ROIzaxis(:,3),'LineWidth',4,'Color','b');