function rangeAddAll(spec,colorScheme)
% rangeAddAll is a convenience function that guides the user through the
% ranging process for all defined ioninc peaks that currently have no range
% associated with them. The focus in the figure is automatically given to
% the current peak and the region around it. 
%
%

%% find individual peaks and determine if a range is defined
% find individual peak locations
ionLocations = [];
plotHandles = spec.Parent.Children;
for pl = 1:length(plotHandles)
    if plotHandles(pl).UserData.plotType == "ion"
        ionLocations = [ionLocations plotHandles(pl).XData];
    end
end
ionLocations = ionLocations';

% compare with ranges
rangeTable = rangesExtractFromMassSpec(spec);

notInRange = not(any((ionLocations > rangeTable.mcbegin') & (ionLocations < rangeTable.mcend'),2));

%% cycle through ions that have no range associated


disp('fin');