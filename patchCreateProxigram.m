function [proxi, binVector] = patchCreateProxigram(posSpecies,pos,interface,bin)
% patchCreateProxigram calculates a proxigram for the patch 'interface' 
% for the atoms in 'posSpecies', which are a subset of the atoms in 'pos' 
% with a bin width of 'bin'.
% 
% [proxi, binVector] = patchCreateProxigram(posSpecies,pos,interface,bin)
%
% INPUT
% posSpecies:   pos of the desired species, subset of initial pos
%
% pos:          initial pos file with pos.x, pos.y, and pos.z
%
% interface:    structure with faces and vertices of the interface
%
% bin:          bin width in nm
%
% OUTPUT
% proxi:        concentration values of the desired species; y-values of
%               the proxigram
%
% binVector:    corresponding distance values; x-values of the proxigram

% distances are calculated along vertex normals.
normals = patchnormals(interface);


%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(interface.vertices,delaunayn(interface.vertices),[pos.x, pos.y, pos.z]);
distVec = [pos.x, pos.y, pos.z] - interface.vertices(closest,:);
% distance through dot product
dist = sum(normals(closest,:) .* distVec,2);


% calculating bin centers
binVector = linspace(0,10000*bin,10001);
binVector = [fliplr(uminus(binVector(2:end))) binVector];
binVector(binVector<min(dist) | binVector>max(dist)) = [];

% number of atoms per bin
posHist = hist(dist,binVector);


%% for element pos files
closestS = dsearchn(interface.vertices,delaunayn(interface.vertices),[posSpecies.x, posSpecies.y, posSpecies.z]);
distVecS = [posSpecies.x, posSpecies.y, posSpecies.z] - interface.vertices(closestS,:);
% distance through dot product
distS = sum(normals(closestS,:) .* distVecS,2);

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

end