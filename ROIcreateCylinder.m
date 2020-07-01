function ch = ROIcreateCylinder(radius,height,location,numSegments,ax)
% creates cylinder in current or parsed axis with specified radius and
% height. 
% Output is handle to the object for later manipulation.
%
% ch = ROIcreateCylinder()
% ch = ROIcreateCylinder(radius)
% ch = ROIcreateCylinder(radius,height)
% ch = ROIcreateCylinder(radius,height,location)
% ch = ROIcreateCylinder(radius,height,location,numSegments)
% ch = ROIcreateCylinder(radius,height,location,numSegments,ax)
%
% INPUTS
% radius:       radius of the cylinder, default is 5   
% height:       height of the cylinder, default is 10
% location:     start coordinates of the cylinder, default is [0 0 0]
% numSegments:  number of segments of the cylinder, default is 32
%
% OUTPUTS
% ch:           handle to the ROIcylinder

%% sets default radius to 5
if not(exist('radius','var'))
    radius = 5;
end
%% sets default height to 10
if not(exist('height','var'))
    height = 10;
end
%% sets default location to [0, 0, 0]
if not(exist('location','var'))
    location = [0 0 0];
end
%% sets default numSegments to 32
if not(exist('numSegments','var'))
    numSegments = 32;
end
%% gets axis
if not(exist('ax','var'))
    ax = gca;
end
%% creates the cylinder
[x, y, z] = cylinder(radius,numSegments);

z = z * height - height/2;
x = x';
y = y';
z = z';
x = x(:);
y = y(:);
z = z(:);
x0 = location(:,1); % getting x coordinate of location
y0 = location(:,2); % getting y coordinate of location
z0 = location(:,3); % getting z coordinate of location
x = x + x0; % calculate shifted x coordinates
y = y + y0; % calculate shifted y coordinates
z = z + z0; % calculate shifted z coordinates

fv.vertices = [x,y,z];

faces = convhull(x,y,z);
fv.faces = faces;

%% plotting the patch object
ch = patch(fv);
ch.FaceColor = [.5 , .5 , .5];
ch.FaceAlpha = 0.5;
ch.DisplayName = 'ROI cylinder';
%% defining reference coordinate system
ch.UserData.ROIzaxis = [location ; location + [0,0,height]];
ch.UserData.ROIyaxis = [location ; location + [0,radius,0]];
ch.UserData.ROIxaxis = [location ; location + [radius,0,0]];




