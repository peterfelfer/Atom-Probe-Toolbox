function [proxi, binvector] = kwikLineProxigram(pos,parentPos,line,bin,vertexIndices)
% calculates a proxigram for the line object 'line' for the atoms in 'pos',
% which are a subset of the atoms in 'parentPos' with a binwidth of bin.
% 'vertexIndices' allows for a local selection

%line is parsed as a struct with line.vertices (Nx3) and line.edges (Mx2)

% distances are calculated perpendicular to the line elements.
[lineVector, lengthOut] = lineVectors(line.vertices,line.edges);



%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(line.vertices,delaunayn(line.vertices),parentPos(:,1:3));

%vector from atom to closest vertex
vec = parentPos(:,1:3) - line.vertices(closest,1:3); 
vecLen = sqrt(sum(vec.^2,2));
%distance along line vector
distLV = dot(vec, lineVector(closest,:), 2);
%distance normal to line vector
dist = sqrt(vecLen.^2 - distLV.^2); 

% calculating bin centers
binvector = linspace(0,10000*bin,10001);
binvector = [fliplr(uminus(binvector(2:end))) binvector];
binvector(binvector<min(dist) | binvector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binvector);




%% for element pos files
closestS = dsearchn(line.vertices,delaunayn(line.vertices),pos(:,1:3));
distVecS = pos(:,1:3) - line.vertices(closestS,:);

%vector from atom to closest vertex
vecS = pos(:,1:3) - line.vertices(closestS,1:3); 
vecLenS = sqrt(sum(vecS.^2,2));
%distance along line vector
distLVS = dot(vecS, lineVector(closestS,:), 2);
%distance normal to line vector
distS = sqrt(vecLenS.^2 - distLVS.^2); 

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