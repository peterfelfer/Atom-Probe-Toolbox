function axisSpatialAptify
% axisSpatialAptify takes the current axes and puts them into an APT style 
% display. It sets the axes and background color, the direction of the 
% z-axis and the plot is ready to rotate the reconstructed APT tip.
% 
% axisSpatialAptify
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

axis equal;
rotate3d on;
view(3);

set(gca,'Box','On');
set(gca,'BoxStyle','full');
set(gca,'XColor',[1 0 0]);
set(gca,'YColor',[0 1 0]);
set(gca,'ZColor',[0 0 1]);
set(gca,'Color',[1 1 1]);

set(gca,'ZDir','reverse');

set(gcf,'Color',[1 1 1]);