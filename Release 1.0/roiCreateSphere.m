function sh = roiCreateSphere(radius,subDivisions,location,ax)
% roiCreateSphere creates sphere in current or parsed axis with specified 
% radius. Output is handle to the object for later manipulation.
% 
% sh = roiCreateSphere()
% sh = roiCreateSphere(radius)
% sh = roiCreateSphere(radius,subDivisions)   
% sh = roiCreateSphere(radius,subDivisions,location)
% sh = roiCreateSphere(radius,subDivisions,location,ax)
%
% INPUT
% radius:       the radius of the ROI, default is 10
%
% subDivisions: fineness of the grid, default is 2 
%
% location:     the start coordinates of the ROI given as [x y z] 
%               default is [0 0 0]
%
% ax:           axes in which the ROI is orientated
%
% OUTPUT
% sh:           handle to the ROIsphere
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

if not(exist('radius','var'))
    radius = 10;
end
if not(exist('subDivisons','var'))
    subDivisions = 2;
end

if not(exist('location','var'))
    location = [0 0 0];
end

if not(exist('ax','var'))
    ax = gca;
end



[fv.vertices fv.faces] = icosphere(subDivisions);
fv.vertices = fv.vertices * radius;

sh = patch(ax,'Vertices',fv.vertices,'Faces',fv.faces);
sh.FaceColor = [.5 , .5 , .5];
sh.FaceAlpha = 0.5;
sh.Vertices(:,1) = sh.Vertices(:,1) + location(:,1);%shifts x coordinates
sh.Vertices(:,2) = sh.Vertices(:,2) + location(:,2);%shifts y coordinates
sh.Vertices(:,3) = sh.Vertices(:,3) + location(:,3);%shifts z coordinates

%% defining reference coordinate system
sh.UserData.ROIzaxis = [location ; location + [0,0,radius]];
sh.UserData.ROIyaxis = [location ; location + [0,radius,0]];
sh.UserData.ROIxaxis = [location ; location + [radius,0,0]];

end
function [vv,ff] = icosphere(varargin)
% ICOSPHERE Generate icosphere.
% Create a unit geodesic sphere created by subdividing a regular
% icosahedron with normalised vertices.
%
%   [V,F] = ICOSPHERE(N) generates two matrices containing vertex and face
%   data so that patch('Faces',F,'Vertices',V) produces a unit icosphere
%   with N subdivisions.
%
%   FV = ICOSPHERE(N) generates an FV structure for using with patch.
%
%   ICOSPHERE(N) and just ICOSPHERE display the icosphere as a patch on the
%   current axes and does not return anything.
%
%   ICOSPHERE uses N = 3.
%
%   ICOSPHERE(AX,...) plots into AX instead of GCA.
%
%   See also SPHERE.
%
%   Based on C# code by Andres Kahler
%   http://blog.andreaskahler.com/2009/06/creating-icosphere-mech-in-code.html
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
% Parse possible axes input
if nargin > 2
    error('Too many input variables, must be 0, 1 or 2.');
end
[cax,args,nargs] = axescheck(varargin{:});
n = 3; % default number of subdivisions
if nargs > 0, n = args{1}; end % override based on input
% generate regular unit icosahedron (20 faced polyhedron)
[v,f] = icosahedron(); % size(v) = [12,3]; size(f) = [20,3];
% recursively subdivide triangle faces
for gen = 1:n
    f_ = zeros(size(f,1)*4,3);
    for i = 1:size(f,1) % for each triangle
        tri = f(i,:);
        % calculate mid points (add new points to v)
        [a,v] = getMidPoint(tri(1),tri(2),v);
        [b,v] = getMidPoint(tri(2),tri(3),v);
        [c,v] = getMidPoint(tri(3),tri(1),v);
        % generate new subdivision triangles
        nfc = [tri(1),a,c;
               tri(2),b,a;
               tri(3),c,b;
                    a,b,c];
        % replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_; % update 
end
% remove duplicate vertices
[v,b,ix] = unique(v,'rows'); clear b % b dummy / compatibility
% reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f),'rows');
switch(nargout)
    case 0 % no output
        cax = newplot(cax); % draw to given axis (or gca)
        chowSphere(cax,f,v);
    case 1 % return fv structure for patch
        vv = struct('Vertices',v,'Faces',f,...
                    'VertexNormals',v,'FaceVertexCData',v(:,3));
    case 2 % return vertices and faces
        vv = v; ff = f;
    otherwise
        error('Too many output variables, must be 0, 1 or 2.');
end
end
function [i,v] = getMidPoint(t1,t2,v)
% GETMIDPOINT calculates point between two vertices
%   Calculate new vertex in subdivision aex
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
% get vertice positions
p1 = v(t1,:); p2 = v(t2,:);
% calculate mid point (on unit sphere)
pm = (p1 + p2) ./ 2;
pm = pm./norm(pm);
% add to vertices list, return index
i = size(v,1) + 1;
v = [v;pm];
end
function [v,f] = icosahedron()
% ICOSAHEDRON creates unit regular icosahedron
%   Returns 12 vertex and 20 face values.
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));
% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20
end
function chowSphere(cax,f,v)
% chOWSPHERE displays icosphere given faces and vertices.
%   Displays a patch surface on the axes, cax, and sets the view.
%   
%   Properties:
%       - vertex normals == vertex vectors
%       - no backface lighting
%       - colour data matches z value of vertex
%       - material properties match default SURF
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK
% set some axes properties if not held
if ~ichold(cax)
    az = -37.5; el = 30;
    view(az,el)
    grid
end
% create patch object on cax
patch('Faces',f,'Vertices',v,...
    'VertexNormals',v,...
    'LineWidth',0.5,'FaceLighting','phong',...
    'BackFaceLighting','unlit',...
    'AmbientStrength',0.3,'DiffuseStrength',0.6,...
    'SpecularExponent',10,'SpecularStrength',0.9,...
    'FaceColor','flat','CData',v(:,3),...
    'Parent',cax,'Tag','Icosphere');
end
