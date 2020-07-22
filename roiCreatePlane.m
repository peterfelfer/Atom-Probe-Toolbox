function ph = roiCreatePlane(dimensions,spacing,location,ax)
% roiCreatePlane creates a plane in the current or parsed axis with 
% specified width and height at location. 
% spacing gives the approximate spacing between vertices for later analysis
% purposes. leave at 0 for simple plane object with no subdivision.
% Output is handle to the object for later manipulation.
%
% ph = roiCreatePlane(dimensions,spacing)
% ph = roiCreatePlane(dimensions,spacing,location)
% ph = roiCreatePlane(dimensions,spacing,location,ax)
%
% INPUT
% dimensions:   1x2 array with x and y dimensions of the plane
%
% spacing:      1x2 array with the spacing value in x and y direction
%
% location:     coordinates of the ROI's center given as [x y] 
%               default is [0 0]
%
% ax:           axes in which the ROI is orientated
%
% OUTPUT
% ph:           handle to the ROIplane

% get axis if necessary
if not(exist('ax','var'))
    ax = gca;
end

% default location is [0, 0, 0]
if ~exist('location','var')
    location = [0, 0, 0];
end

% x and y coordinates of the vertices
xLocs = (0:spacing(:,1):dimensions(:,1)) - dimensions(:,1)/2;
yLocs = (0:spacing(:,2):dimensions(:,2)) - dimensions(:,2)/2;
[tx, ty] = meshgrid(xLocs,yLocs);
numVerts = length(tx(:));
vertices = [tx(:), ty(:)];
vertices(:,3) = 0;
vertices = vertices + repmat(location,numVerts,1);

% delaunay triangulation of vertices
faces = delaunay(vertices(:,1),vertices(:,2));

% plotting of patch object
ph = patch(ax,'Vertices',vertices,'Faces',faces);
ph.FaceColor = [.5 , .5 , .5];
ph.FaceAlpha = 0.5;
ph.DisplayName = 'ROI plane';

% defining reference coordinate system
ph.UserData.ROIzaxis = [location ; location + [0,0,dimensions(1)]];
ph.UserData.ROIyaxis = [location ; location + [0,dimensions(2),0]];
ph.UserData.ROIxaxis = [location ; location + [mean(dimensions),0,0]];