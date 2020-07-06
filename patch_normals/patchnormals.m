function N=patchnormals(FV)
% This function PATCHNORMALS calculates the normals of a triangulated
% mesh. PATCHNORMALS calls the patchnormal_double.c mex function which 
% first calculates the normals of all faces, and after that calculates 
% the vertice normals from the face normals weighted by the angles 
% of the faces.
%
% N=patchnormals(FV);
%
% Inputs,
%   FV : A struct containing FV.faces with a facelist Nx3 and FV.vertices
%        with a Nx3 vertices list. Such a structure is created by Matlab
%        Patch function
% Outputs,
%   N : A Mx3 list with the normals of all vertices
%
% Example,
%   % Compile the c-coded function
%   mex patchnormals_double.c -v
%
%   % Load a triangulated mesh of a sphere
%   load sphere; 
%
%   % Calculate the normals
%   N=patchnormals(FV);
%
%   % Show the normals
%   figure, patch(FV,'FaceColor',[1 0 0]); axis square; hold on;
%   for i=1:size(N,1);
%       p1=FV.vertices(i,:); p2=FV.vertices(i,:)+10*N(i,:);       
%       plot3([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)],'g-');
%   end       
%
% Function is written by D.Kroon University of Twente (June 2009)


sizev=size(FV.vertices);
% Check size of vertice array
if((sizev(2)~=3)||(length(sizev)~=2))
    error('patchnormals:inputs','The vertice list is not a m x 3 array')
end

sizef=size(FV.faces);
% Check size of vertice array
if((sizef(2)~=3)||(length(sizef)~=2))
    error('patchnormals:inputs','The vertice list is not a m x 3 array')
end

% Check if vertice indices exist
if(max(FV.faces(:))>size(FV.vertices,1))
    error('patchnormals:inputs','The face list contains an undefined vertex index')
end

% Check if vertice indices exist
if(min(FV.faces(:))<1)
    error('patchnormals:inputs','The face list contains an vertex index smaller then 1')
end

[Nx,Ny,Nz]=patchnormals_double(double(FV.faces(:,1)),double(FV.faces(:,2)),double(FV.faces(:,3)),double(FV.vertices(:,1)),double(FV.vertices(:,2)),double(FV.vertices(:,3)));

N=zeros(length(Nx),3);
N(:,1)=Nx;
N(:,2)=Ny;
N(:,3)=Nz;

function [Nx,Ny,Nz]=patchnormals_double(Fa,Fb,Fc,Vx,Vy,Vz)
%
%  [Nx,Ny,Nz]=patchnormals_double(Fa,Fb,Fc,Vx,Vy,Vz)
%

FV.vertices=zeros(length(Vx),3);
FV.vertices(:,1)=Vx;
FV.vertices(:,2)=Vy;
FV.vertices(:,3)=Vz;

% Get all edge vectors
e1=FV.vertices(Fa,:)-FV.vertices(Fb,:);
e2=FV.vertices(Fb,:)-FV.vertices(Fc,:);
e3=FV.vertices(Fc,:)-FV.vertices(Fa,:);

% Normalize edge vectors
e1_norm=e1./repmat(sqrt(e1(:,1).^2+e1(:,2).^2+e1(:,3).^2),1,3); 
e2_norm=e2./repmat(sqrt(e2(:,1).^2+e2(:,2).^2+e2(:,3).^2),1,3); 
e3_norm=e3./repmat(sqrt(e3(:,1).^2+e3(:,2).^2+e3(:,3).^2),1,3);

% Calculate Angle of face seen from vertices
Angle =  [acos(dot(e1_norm',-e3_norm'));acos(dot(e2_norm',-e1_norm'));acos(dot(e3_norm',-e2_norm'))]';

% Calculate normal of face
 Normal=cross(e1,e3);

% Calculate Vertice Normals 
VerticeNormals=zeros([size(FV.vertices,1) 3]);
for i=1:size(Fa,1),
    VerticeNormals(Fa(i),:)=VerticeNormals(Fa(i),:)+Normal(i,:)*Angle(i,1);
    VerticeNormals(Fb(i),:)=VerticeNormals(Fb(i),:)+Normal(i,:)*Angle(i,2);
    VerticeNormals(Fc(i),:)=VerticeNormals(Fc(i),:)+Normal(i,:)*Angle(i,3);
end

V_norm=sqrt(VerticeNormals(:,1).^2+VerticeNormals(:,2).^2+VerticeNormals(:,3).^2)+eps;
VerticeNormals=VerticeNormals./repmat(V_norm,1,3);
Nx=VerticeNormals(:,1);
Ny=VerticeNormals(:,2);
Nz=VerticeNormals(:,3);