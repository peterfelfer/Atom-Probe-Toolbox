function [p, ax, f] = concentrationProfilePlot(conc,excludeList, colorScheme)
% plots a concentration profile
% plotType can be inLine or stacked
% up to 4 volumes with different plot lines
%
% INPUT
% conc:        table that contains the counts or the concentration of 
%              one or more volumes. If multiple volumes are in the 
%              variable, the name will be atom/ion + volume name
%
% excludeList: cell array that contains as character the individual
%              ions that shall not be considered for the plot of the 
%              concentration, unranged atoms appear as 'unranged', if not 
%              parsed, no atoms will be excluded
%
% colorScheme: table with RGB codes assigned to ions
%
% OUTPUT
% p:            bar/pie chart of the volume, with properties (BarLayout, 
%               BarWidth, FaceColor, EdgeColor, BaseValue, XData, YData)
%
% ax:           axes of the plot with properties (XLim, YLim, XScale, 
%               YScale, GridlineStyle, Position, Units)
%
% f:            figure that contains the plot (pie or bar chart) with 
%               properties(Number, Name, Color, Position, Units)
%
%
%
%
% USEFUL Notes:
% If the conc as OUTPUT of posCalculateConcentrationSimple function is used, either
% concentration or counts must be specified in the INPUT argument with the
% following lines as INPUT argument
%       concentration:  conc([conc.format=='concentration'], :)
%       count:         conc([conc.format=='count'], :)
%

if ~exist('excludeList','var')
    excludeList = {};
end

% remove elements on the exclude list
conc = conc(:,~ismember(conc.Properties.VariableNames,excludeList));


%% check for multiple volumes, presence of variance for error bars, options and compatibility
% check if all variables have the same format
if ~xor(any(conc.format == 'concentration'),any(conc.format == 'count'))
    error('only either concentration or count format is allowed as input');
end


% find the individual volumes in the concentration variable
volumes = categories(removecats(conc.volume));
numVolumes = length(volumes);

isColorScheme = exist('colorScheme','var');


%% set up figure and axis
f = figure();
f.Color = 'w';

ax = axes(f);
ax.Box = 'on';
hold on


if numVolumes == 1
    f.Name = [char(conc.volume(1)) ' ' char(conc.type(1)) ' ' char(conc.format(1)) ' plot'];
else
    f.Name = [char(conc.type(1)) ' ' char(conc.format(1)) ' plot'];
end

ax.Title.String = f.Name;
ax.YLabel.String = [char(conc.type(1)) ' ' char(conc.format(1))];
try
    ax.XLabel.String = ['distance [' conc.Properties.VariableUnits{strcmp(conc.Properties.VariableNames,'distance')} ']'];
catch
    ax.XLabel.String = 'distance';
end

% check if we are plotting concentrations, need to adjust data
if any(conc.format == 'concentration')
    pctFac = 100;
    ax.YLabel.String = [ax.YLabel.String ' [%]'];
else
    pctFac = 1;
end



% find which columns are atoms/ions
plotIdx = find(ismember(conc.Properties.VariableDescriptions,{'atom','ion'}));

% all categories that have a value of 0 will be deleted
isNonZero = sum(table2array(conc(:,plotIdx))) ~= 0;
plotIdx = plotIdx(isNonZero);


for pl = 1:length(plotIdx)
    p = plot(conc.distance,...
        table2array(conc(:,plotIdx(pl))) * pctFac,...
        'Color',colorScheme.color(colorScheme.ion == conc.Properties.VariableNames{plotIdx(pl)},:),...
        'DisplayName',conc.Properties.VariableNames{plotIdx(pl)},...
        'LineWidth',1);
    
end

legend();