function rangeAddAll(spec,colorScheme,rangeMargin,useMin)
% rangeAddAll is a convenience function that guides the user through the
% ranging process for all defined ionic peaks that currently have no range
% associated with them. The focus in the figure is automatically given to
% the current peak and the region around it. 
% When the peak is in focus, press enter for ranging, any other key to skip
% the current peak and the region around it.
%
% rangeAddAll(spec,colorScheme,rangeMargin)
%
% INPUT
% spec:         figure, that holds the mass spectra with the peaks which 
%               will be ranged
%
% colorScheme:  color Scheme with colors for the ranges
%
% rangeMargin:  is the number of Da before and after the peak that are
%               displayed
% useMin:

xPlotLimits = spec.Parent.XLim; % current plot limits
yPlotLimits = spec.Parent.YLim; % current plot limits



%% find individual peaks and determine if a range is defined
% find individual peak locations
ionLocations = [];
peakHeights = [];
plotHandles = spec.Parent.Children;
for pl = 1:length(plotHandles)
    if plotHandles(pl).UserData.plotType == "ion"
        ionLocations = [ionLocations plotHandles(pl).XData];
        peakHeights = [peakHeights plotHandles(pl).YData];
    end
end
ionLocations = ionLocations';
peakHeights = peakHeights';

% compare with ranges
rangeTable = rangesExtractFromMassSpec(spec);

%check which ions currently have no range associated with them
if ~isempty(rangeTable)
    % for any populated range table 
    notInRange = not(any((ionLocations > rangeTable.mcbegin') & (ionLocations < rangeTable.mcend'),2));
else
    %if no ranges exist in the mass spectrum
    notInRange = true(size(ionLocations));
end


%check which peaks are higher than the defined background
isAboveLimit = true(size(ionLocations));
if useMin == true
    limPoints = ginput(); % prompts the user to input points below which peaks are not considered
    limPoints = sortrows(limPoints);
    % extrapolate to limits of mass spectrum as constant
    limPoints = [spec.XData(1) limPoints(1,2) ; limPoints; spec.XData(end) limPoints(end,2)];
    limPlotHandle = line(spec.Parent,limPoints(:,1),limPoints(:,2),'Color',[1 0 0],'LineWidth',2); % line plot of the background selection
    % get limit y values for ion peak locations
    limVals = interp1(limPoints(:,1),limPoints(:,2),ionLocations,'linear');
    % compare to peak heights
    isAboveLimit = peakHeights > limVals;
end

unrangedIonLocations = ionLocations(notInRange & isAboveLimit);
unrangedPeakHeights = peakHeights(notInRange & isAboveLimit);

%% cycle through ions that have no range associated
for ion = 1:length(unrangedIonLocations)
    % put focus on peak
    spec.Parent.XLim = [unrangedIonLocations(ion) - rangeMargin, unrangedIonLocations(ion) + rangeMargin];
    spec.Parent.YLim = [yPlotLimits(1), unrangedPeakHeights(ion)*1.5];
    
    % put a marker there
    
    % range the peak
    rangeAdd(spec,colorScheme);
end

spec.Parent.XLim = xPlotLimits; %revert to old plot limits.
spec.Parent.YLim = yPlotLimits; %revert to old plot limits.
delete(limPlotHandle); %deletes background line
