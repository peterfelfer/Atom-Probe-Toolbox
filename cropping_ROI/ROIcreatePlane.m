function ph = ROIcreatePlane(dimensions,spacing,location,ax)
%creates plane in current or parsed axis with specified width and
%height at the location loc. 
%spcing gives the approximate spacing between vertices for later analysis
%purposes. leave at 0 for simple plane object with no subdivision.
%Output is handle to the object for later manipulation.

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

    
% plotting of pathc object
ph = patch(ax,'Vertices',vertices,'Faces',faces);
ph.FaceColor = [.5 , .5 , .5];
ph.FaceAlpha = 0.5;
ph.DisplayName = 'ROI box';

% defining reference coordinate system
ph.UserData.ROIzaxis = [location ; location + [0,0,dimensions(1)]];
ph.UserData.ROIyaxis = [location ; location + [0,dimensions(2),0]];
ph.UserData.ROIxaxis = [location ; location + [mean(dimensions),0,0]];