function wcsh = wcsCreate(size,location,ax)
% wcsCreate creates box in current or parsed axis with specified width and
% height at the location loc. 
% Output is handle to the object for later manipulation
%
% wcsh = wcsCreate(size,loc,ax)
% 
% INPUT
% size:     the length, width and height (x,y,z) of the ROI given as 
%           [x y z], default is [10 10 10]
%
% location: the start coordinates of the ROI given as [x y z],
%           default is [0 0 0]
%
% ax:       axes in which the ROI is orientated
% 
% OUTPUT
% points:       points transformed to the new coordinate system
%
%% default size is [10 10]
if not(exist('size','var'))
    size = [10 10 10];
end
%% default location is [0, 0, 0]
if ~exist('location','var')
    location = [0, 0, 0];
end
%% sets axis
if not(exist('ax','var'))
    ax = gca;
end

vertices = [...
    0 0 0;... idx 1 for x axis coloring
    1 0 0;... idx 2 
    1 1 0;... idx 3
    0 1 0;... idx 4 for y axis coloring
    0 0 1;... idx 5 for z axis coloring
    1 0 1;... idx 6
    1 1 1;... idx 7
    0 1 1]; % idx 8
vertices = vertices * size; % scale
vertices = vertices + repmat(location,8,1);

faces = [...
    1 2 3 4;...
    8 4 3 7;...
    1 2 6 5;...
    8 5 6 7;...
    6 7 3 2;...
    4 1 5 8;...
    ];
    

%ch = patch(fv);
wcsh = patch(ax,'Vertices',vertices,'Faces',faces);
wcsh.FaceColor = [.5 , .5 , .5];
wcsh.FaceAlpha = 0;
wcsh.UserData.ROIzaxis = [0,0,0 ; 0,0,size];
wcsh.UserData.ROIyaxis = [0,0,0 ; 0,size,0];
wcsh.UserData.ROIxaxis = [0,0,0 ; size,0,0];
wcsh.DisplayName = 'WCS';

% coloring of vertices
cols = [...
    1, 0, 0;... x axis coloring
    0, 0, 0;... 
    0, 0, 0;...
    0, 1, 0;... y axis coloring
    0, 0, 1;... z axis coloring
    0, 0, 0;...
    0, 0, 0;...
    0, 0, 0;];
wcsh.FaceVertexCData = cols;
wcsh.EdgeColor = 'flat';
wcsh.LineWidth = 4;

wcsh.EdgeAlpha = 'flat';
alpha = [...
    1;
    0;
    0;
    1;
    1;
    0;
    0;
    0;
    ];
wcsh.FaceVertexAlphaData = alpha;
wcsh.LineStyle = '-.';















