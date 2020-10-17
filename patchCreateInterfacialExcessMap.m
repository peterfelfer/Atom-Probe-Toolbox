function IEmap = patchCreateInterfacialExcessMap(pos,interface,lim,deloc)
% NEEDS REWORK AND DOCUMENTATION
%   INPLEMENT BACKGROUND PICKING, INDIVIDUAL BACKGROUNDS
% calculates an IE map for the patch 'interface' for the atoms in 'pos'
% within 'lim' nm of the interface
DEBUG = false;

%%reads a pos file [x,y,z] converted to a Matlab variable and a vertex file
%%[x,y,z] and assigns every atom to the closest vertex.

%vap: m y z mc vert# disttovert distalongnormal( or to line element)
%shiftdistance
%vertex: x y z obj# nx ny nz d A(or l)

addpath('patch_normals');
%addpath('./general/resources/pos_tools');
addpath('dualMesh');
%addpath('resources');



normals = patchnormals(interface);


%% tessellation and distance clipping


numAtom = length(pos(:,1));
numVerts = length(interface.vertices);


% finding closest point for each atomic position
closest = dsearchn(interface.vertices,delaunayn(interface.vertices),pos(:,1:3));

distVec = pos(:,1:3) - interface.vertices(closest,:);



%distance through dot product
dist = sum(normals(closest,:) .* distVec,2);
dist = abs(dist);

%pos = pos(dist <= lim);
closest = closest(dist<=lim);


if DEBUG
    hist(dist,20);
end

IEcount = zeros(numVerts,1);

for v = 1:numVerts
    % raw count of IE atoms per vertex
    IEcount(v) = sum(closest == v);
end




%% calculation of vertex area
[cp,ce,pv,ev] = makedual2(interface.vertices,interface.faces);
[pc,area] = geomdual2(cp,ce,pv,ev);

IEmap = IEcount./area;


%% delocalisation step


if exist('deloc','var')
    for d = 1:deloc
        IEmap = delocalizeProperty(interface,IEmap);
    end
end
%IEmap = IEmap';

%% export to *.ply?
%property gets normalized

minIE = min(IEmap);
maxIE = max(IEmap);

limits = ['comment interfacial excess map ranges: ' num2str(minIE,3) ' at/nm2 to ' num2str(maxIE,3) ' at/nm2'];


IEmapPly = IEmap - minIE;
IEmapPly = IEmapPly/(maxIE-minIE);

[file path] = uiputfile('*.ply','save quick IE map as ply','IEmap');

if file
    patch2ply(interface,[IEmapPly, IEmapPly, IEmapPly],[path file], limits);
end

%% visualising the results
interface.facevertexcdata = IEmap;

f = figure('Name','Interfacial excess map');
trisurf(interface.faces,interface.vertices(:,2),interface.vertices(:,1),interface.vertices(:,3),IEmap);
axis equal;
rotate3d on;
shading interp;
colorbar;


end







