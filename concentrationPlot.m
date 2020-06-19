function [p, ax, f] = concentrationPlot(conc,excludeList, plotType, colorScheme)
% concentrationPlot plots the concentration given by the variable conc. It
% excludes all the atoms that are on the excludeList. The plot type can be
% a pie or bar chart. The coloring is according to the colorScheme.
%
% [p, ax, f] = concentrationPlot(conc, excludeList, plotType, colorScheme)
% [p, ax, f] = concentrationPlot(conc, excludeList, plotType)
% [p, ax, f] = concentrationPlot(conc, excludeList)
% [p, ax, f] = concentrationPlot(conc)
%
% INPUT
% conc:         is a table that contains the counts or the concentration of 
%               one or more volumes. If multiple volumes are in the 
%               variable, the name will be atom/ion + volume name
%
% excludeList:  is a cell array that contains as character the individual
%               ions that shall not be considered for the plot of the 
%               concentration, unranged atoms appear as 'unranged', if not 
%               parsed, no atoms will be excluded
%
% plotType:     can be a 'pie' or 'bar'
%
% colorScheme:  coloring will be according to colorScheme
%               If multiple volumes are parsed and you want to color by 
%               volume, than don't parse any colorScheme, default will 
%               color the bars by volume.
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
% USEFUL Notes:
% If the conc as OUTPUT of posCalculateConcentrationSimple function is used, either
% concentration or counts must be specified in the INPUT argument with the
% following lines as INPUT argument
%       concentration:  conc([conc.format=='concentration'], :)
%       count:         conc([conc.format=='count'], :)
%
% Display the bar plot with a log scale use
%   ax.YScale = 'log';
%
%
% missing: reduce distance between groups of bars for multiple volumes


%% default plot type is 'bar'
if ~exist('plotType','var')
    plotType = 'bar';
end

if ~exist('excludeList','var')
    excludeList = {};
end
% remove elements on the exclude list
conc = conc(:,~ismember(conc.Properties.VariableNames,excludeList(:)));

%% check for multiple volumes, presence of variance for error bars, options and compatibility
% check if all variables have the same format
if ~xor(any(conc.format == 'concentration'),any(conc.format == 'count'))
    error('only either concentration or count format is allowed as input');
end

volumes = categories(removecats(conc.volume));
numVolumes = length(volumes);

isColorScheme = exist('colorScheme','var');



%% create figure
f = figure;
f.Color = 'w';
ax = axes(f);
if numVolumes == 1
    f.Name = [char(conc.volume(1)) ' ' char(conc.type(1)) ' ' char(conc.format(1)) ' plot'];
else
    f.Name = [char(conc.type(1)) ' ' char(conc.format(1)) ' plot'];
end
% find which columns are atoms/ions
isPlot = ismember(conc.Properties.VariableDescriptions,{'atom','ion'});
plotTable = conc(:,isPlot);


% all categories that have a value of 0 will be deleted
y = table2array(plotTable);
isZero = y(1,:) == 0;
y = y(:,~isZero);


% categories in the categorical need to be ordered, otherwise plot would be
% sorted alphabetical
x = categorical(plotTable.Properties.VariableNames(~isZero(1,:)));
x = reordercats(x, plotTable.Properties.VariableNames(~isZero(1,:)));


%% create actual plots
% bar plot
if strcmp(plotType,'bar')
    ytrans = y'; % bar needs y in a row format
    if conc.format == 'concentration' % display concentration
        p = bar(x,ytrans*100);
        ax.YLabel.String = [char(conc.format(1)) ' [' char(conc.type(1)) ' %]'];
    else
        p = bar(x,ytrans); % display counts
        type = char(conc.type);
        type = type(1:end-2); % loose the 'ic' in 'atomic'
        ax.YLabel.String = [char(conc.format) ' [' type 's]'];
    end

    for pl = 1:length(p)
        p(pl).FaceColor = 'flat';
        for b = 1:length(p(pl).XData) % go through individual bars
            if isColorScheme
                col(b,:) = colorScheme.color(colorScheme.ion == char(p(pl).XData(b)),:);
                p(pl).CData(b,:) = col(b,:);
            end
            p(pl).DisplayName = volumes{pl};

        end
    end

% bar plot
elseif strcmp(plotType,'pie')
    if numVolumes > 1
        delete(f);
        error('pie plot for multiple volumes not possible');
    end
    % pie plot
    p = pie(y);
    for b = 1:length(y) % go through individual bars
        % pie charts are collections of patches. A text and a pie slice for
        % each category, i.e. number of categories * 2 plothandles
        if isColorScheme
            col(b,:) = colorScheme.color(colorScheme.ion == char(x(b)),:);
            p(b*2-1).FaceColor = col(b,:);
        end
        p(b*2-1).DisplayName = char(x(b));
        legend();
    end


end

%% set plot axis title
% needs to be done after plotting, otherwise title will be reset by plot
if length(p) == 1
    ax.Title.String = [char(conc.volume) ' ' char(conc.type) ' ' char(conc.format)];
else
    volNames = char(join(string(volumes)));
    ax.Title.String = [volNames ' ' char(conc.type(1)) ' ' char(conc.format(1))];
end
