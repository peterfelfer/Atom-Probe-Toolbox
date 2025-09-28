function rangeAddAll(spec,colorScheme,rangeMargin,useMin)
% rangeAddAll is a convenient function that guides the user through the
% ranging process for all defined ionic peaks that currently have no range
% associated with them. The focus in the figure is automatically given to
% the current peak and the region around it. 
% When the peak is in focus, the user needs to range it manually. If the 
% peak should be skipped, enter must be pressed.
%
% rangeAddAll(spec,colorScheme,rangeMargin)
%
% INPUT
% spec:         figure, that holds the mass spectrum with the peaks which 
%               will be ranged
%
% colorScheme:  color scheme with colors for the ranges
%
% rangeMargin:  is the number of Da before and after the peak that are
%               displayed
% useMin:       the user can choose between 'true' and 'false'
%               false:  the user will be taken through every peak
%               true: the user can draw a line over the mass spectrum and
%               only the peaks higher than the line are taken into account
%               for the ranging
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

xPlotLimits = spec.Parent.XLim; % current plot limits
yPlotLimits = spec.Parent.YLim; % current plot limits

%% check for color of the ion in the color Scheme

ionTable = ionsExtractFromMassSpec(spec);
% find categories in ionTable and colorScheme
testIonTable =  categories(ionTable.ionName);
testColorScheme = categories(colorScheme.ion);
% check if ion is in colorScheme
       compIonColor = zeros(size(testIonTable));
        for k = 1:height(testColorScheme)
        compIonColor(strncmp(testIonTable, testColorScheme(k,1), 15)) = k;
        end
        
        % find zeros = no ion defined
        posZeros = find(~compIonColor);
        if posZeros > 0
            error('Ions in the mass spectra have no defined color in the colorScheme. Please add the ion to the colorScheme with colorSchemeIonAdd()');
        end 
        
        




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
    
    if unrangedPeakHeights(ion)*1.5 > yPlotLimits(1) % in case the upper plot y limit is lower than the lower y plot limit
        spec.Parent.YLim = [yPlotLimits(1), unrangedPeakHeights(ion)*1.5];
    else
        spec.Parent.YLim = [unrangedPeakHeights(ion)/10, unrangedPeakHeights(ion)*1.5];
    end
    
    % put a marker there
    
    % range the peak
    rngH = rangeAdd(spec,colorScheme);
    
    %if isempty(rngH)
    %    return
    %end
end

spec.Parent.XLim = xPlotLimits; %revert to old plot limits.
spec.Parent.YLim = yPlotLimits; %revert to old plot limits.
delete(limPlotHandle); %deletes background line
