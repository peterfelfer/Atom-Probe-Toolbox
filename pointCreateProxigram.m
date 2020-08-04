function [proxi, binVector] = pointCreateProxigram(posSpecies,pos,vertices,bin)
% pointCreateProxigram calculates a proxigram for the input vertices 
% for the atoms in 'posSpecies', which are a subset of the atoms in 'pos' 
% with a bin width of 'bin'.
%
% INPUT
% posSpecies:   pos of the desired species, subset of initial pos
%
% pos:          initial pos file with pos.x, pos.y, and pos.z
%
% vertices:     vertices, m-by-n matrix, with m>3
%
% bin:          bin width in nm
%
% OUTPUT
% proxi:        concentration values of the desired species; y-values of
%               the proxigram
%
% binVector:    corresponding distance values; x-values of the proxigram
%

vertices = vertices(:,1:3);


%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(vertices,delaunayn(vertices),[pos.x, pos.y, pos.z]);

% vector from atom to closest vertex
vec = [pos.x, pos.y, pos.z] - vertices(closest,1:3); 
dist = sqrt(sum(vec.^2,2)); 

% calculating bin centers
binVector = linspace(0,10000*bin,10001);
binVector = [fliplr(uminus(binVector(2:end))) binVector];
binVector(binVector<min(dist) | binVector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binVector);




%% for element pos files
closestS = dsearchn(vertices,delaunayn(vertices),[posSpecies.x, posSpecies.y, posSpecies.z]);
distVecS = [posSpecies.x, posSpecies.y, posSpecies.z] - vertices(closestS,:);

% vector from atom to closest vertex
vecS = [posSpecies.x, posSpecies.y, posSpecies.z] - vertices(closestS,1:3); 

distS = sqrt(sum(distVecS.^2,2)); 

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