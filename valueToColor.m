function [vertColors colBar] = valueToColor(values,colorMap, limits)
% NEEDS DOCUMENTATION
% valueToColor maps a range of values to a colormap. This is needed for
% visualisations of property maps on mesh objects
% translates the values in vals (1 x N matrix) into vertex colors using the
% colormap 'cmap'. Default is jet. Also outputs the colorbar colorBar,
% which is a clolrmapped patch that can be used in the actual visualisation
% program.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
if max(values) < 0
    values = uminus(values);
end


normalisedValues = (values - min(limits))/(max(limits) - min(limits));
cLevels = length(colorMap) -1; % 8 bit color is fully sampled
idx = uint8( round( normalisedValues * cLevels ))+1;
vertColors = colorMap(idx,:);

% creating colorbar

D = 1; %width of colorbar
H = 10; %height of colorbar

z = linspace(0,H,cLevels)';

fv.vertices = [zeros(size(z)), zeros(size(z)), z; ones(size(z))*D, zeros(size(z)), z];
fv.faces = delaunay(fv.vertices(:,1),fv.vertices(:,3));

%choice = questdlg('Export color bar to *.ply?','ply export','yes','no','no');

if false
    name = ['colorbar_' num2str(min(values),3) '_to_' num2str(max(values),3) '.ply'];
    [file path] = uiputfile('*.ply','save colorbar to:',name);
    patch2ply(fv,[colorMap; colorMap;],[path file]);
end   

colBar.fv = fv;
colBar.colorValues = [colorMap; colorMap];