function [proxi, binvector] = kwikPointProxigram(pos,parentPos,vertices,bin,vertexIndices)
% calculates a proxigram for the line object 'line' for the atoms in 'pos',
% which are a subset of the atoms in 'parentPos' with a binwidth of bin.
% 'vertexIndices' allows for a local selection

%line is parsed as a struct with line.vertices (Nx3) and line.edges (Mx2)

% distances are calculated perpendicular to the line elements.

vertices = vertices(:,1:3);


%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(vertices,delaunayn(vertices),[parentPos.x, parentPos.y, parentPos.z]);

% vector from atom to closest vertex
vec = [parentPos.x, parentPos.y, parentPos.z] - vertices(closest,1:3); 
dist = sqrt(sum(vec.^2,2)); 

% calculating bin centers
binvector = linspace(0,10000*bin,10001);
binvector = [fliplr(uminus(binvector(2:end))) binvector];
binvector(binvector<min(dist) | binvector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binvector);




%% for element pos files
closestS = dsearchn(vertices,delaunayn(vertices),[pos.x, pos.y, pos.z]);
distVecS = [pos.x, pos.y, pos.z] - vertices(closestS,:);

% vector from atom to closest vertex
vecS = [pos.x, pos.y, pos.z] - vertices(closestS,1:3); 

distS = sqrt(sum(distVecS.^2,2)); 

% number of atoms per bin
proxi = hist(distS,binvector)./posHist;


%% plotting
f = figure;
hold all;

plot(binvector,proxi*100);
set(gcf,'Name','proximity histogram');
set(gcf,'Color',[1 1 1]);
set(get(gca,'XLabel'),'String','distance [nm]');
set(get(gca,'YLabel'),'String','concentration [%]');