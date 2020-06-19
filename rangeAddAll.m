function rangeAddAll(spec,colorScheme,rangeMargin)
% rangeAddAll is a convenience function that guides the user through the
% ranging process for all defined ioninc peaks that currently have no range
% associated with them. The focus in the figure is automatically given to
% the current peak and the region around it. 
%
% rangeMargin is the number of Da before and after the peak that are
% displayed
%
% When the peak is in focus, press enter for ranging, any other key to skip
%

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

% compare with ranges
rangeTable = rangesExtractFromMassSpec(spec);

notInRange = not(any((ionLocations > rangeTable.mcbegin') & (ionLocations < rangeTable.mcend'),2));
unrangedIonLocations = ionLocations(notInRange);
unrangedPeakHeights = peakHeights(notInRange);



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
