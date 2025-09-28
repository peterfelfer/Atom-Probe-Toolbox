function [proxi, binVector] = lineCreateProxigram(posSpecies,pos,line,bin)
% lineCreateProxigram calculates a proxigram for the line object 'line' for
% the atoms in 'posSpecies', which are a subset of the atoms in 'pos' with 
% a bin width of 'bin'.
%
% line is parsed as a struct with line.vertices (Nx3) and line.edges (Mx2)
%
% INPUT
% posSpecies:   pos of the desired species, subset of initial pos
%
% pos:          initial pos file with pos.x, pos.y, and pos.z
%
% line:         structure, with fields of vertices (Nx3) and edges (Mx2)
%
% bin:          bin width in nm
%
% OUTPUT
% proxi:        concentration values of the desired species; y-values of
%               the proxigram
%
% binVector:    corresponding distance values; x-values of the proxigram
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


% distances are calculated perpendicular to the line elements.
lineVector = lineVectors(line.vertices,line.edges);



%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(line.vertices,delaunayn(line.vertices),[pos.x, pos.y, pos.z]);

% vector from atom to closest vertex
vec = [pos.x, pos.y, pos.z] - line.vertices(closest,1:3); 
vecLen = sqrt(sum(vec.^2,2));
%distance along line vector
distLV = dot(vec, lineVector(closest,:), 2);
%distance normal to line vector
dist = sqrt(vecLen.^2 - distLV.^2); 

% calculating bin centers
binVector = linspace(0,10000*bin,10001);
binVector = [fliplr(uminus(binVector(2:end))) binVector];
binVector(binVector<min(dist) | binVector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binVector);




%% for element pos files
closestS = dsearchn(line.vertices,delaunayn(line.vertices),[posSpecies.x, posSpecies.y, posSpecies.z]);
distVecS = [posSpecies.x, posSpecies.y, posSpecies.z] - line.vertices(closestS,:);

% vector from atom to closest vertex
vecS = [posSpecies.x, posSpecies.y, posSpecies.z] - line.vertices(closestS,1:3); 
vecLenS = sqrt(sum(vecS.^2,2));
%distance along line vector
distLVS = dot(vecS, lineVector(closestS,:), 2);
%distance normal to line vector
distS = sqrt(vecLenS.^2 - distLVS.^2); 

% number of atoms per bin
proxi = hist(distS,binVector)./posHist;


%% plotting
f = figure;
hold all;

plot(binVector,proxi*100);
set(gcf,'Name','proximity histogram');
set(gcf,'Color',[1 1 1]);
set(get(gca,'XLabel'),'String','distance [nm]');
set(get(gca,'YLabel'),'String','concentration [%]');