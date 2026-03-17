function rangeAddAll(spec,colorScheme,rangeMargin,useMin,minPeakSeparation)
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
% minPeakSeparation:
%               optional minimum m/c separation [Da] to treat peaks as
%               separate targets. Peaks closer than this are merged and
%               ranged once. If omitted, one spectrum bin is used.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

if nargin < 2 || isempty(colorScheme)
    if evalin('base', 'exist(''colorScheme'', ''var'')')
        colorScheme = evalin('base', 'colorScheme');
    else
        error('rangeAddAll:missingColorScheme', ...
            'No colorScheme provided and none found in workspace.');
    end
end
if nargin < 3 || isempty(rangeMargin)
    rangeMargin = 1;
end
if nargin < 4 || isempty(useMin)
    useMin = false;
end
if nargin < 5 || isempty(minPeakSeparation)
    minPeakSeparation = estimateSpecBinWidth(spec);
end

useMin = logical(useMin);
rangeMargin = double(rangeMargin);
minPeakSeparation = double(minPeakSeparation);
if ~isfinite(minPeakSeparation) || minPeakSeparation < 0
    minPeakSeparation = estimateSpecBinWidth(spec);
end
rangeTol = max(estimateSpecBinWidth(spec) / 2, eps);

xPlotLimits = spec.Parent.XLim; % current plot limits
yPlotLimits = spec.Parent.YLim; % current plot limits
limPlotHandle = gobjects(0);
peakMarker = gobjects(0);

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
        if ~isempty(posZeros)
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
    notInRange = not(any((ionLocations >= (rangeTable.mcbegin' - rangeTol)) & ...
        (ionLocations <= (rangeTable.mcend' + rangeTol)),2));
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

% Merge near-identical peaks to avoid duplicate prompting.
[unrangedIonLocations, unrangedPeakHeights] = mergeNearbyPeaks( ...
    unrangedIonLocations, unrangedPeakHeights, minPeakSeparation);

% Always progress from low to high mass-to-charge.
[unrangedIonLocations, sortIdx] = sort(unrangedIonLocations, 'ascend');
unrangedPeakHeights = unrangedPeakHeights(sortIdx);

%% cycle through ions that have no range associated
for ion = 1:length(unrangedIonLocations)
    % Skip if the peak is already covered by a range created earlier.
    if isPeakInExistingRange(spec, unrangedIonLocations(ion), rangeTol)
        fprintf('rangeAddAll: skipped already-ranged peak at %.5f Da.\n', unrangedIonLocations(ion));
        continue;
    end

    % put focus on peak
    spec.Parent.XLim = [unrangedIonLocations(ion) - rangeMargin, unrangedIonLocations(ion) + rangeMargin];
    
    if unrangedPeakHeights(ion)*1.5 > yPlotLimits(1) % in case the upper plot y limit is lower than the lower y plot limit
        spec.Parent.YLim = [yPlotLimits(1), unrangedPeakHeights(ion)*1.5];
    else
        spec.Parent.YLim = [unrangedPeakHeights(ion)/10, unrangedPeakHeights(ion)*1.5];
    end
    
    % mark the current peak to avoid ambiguity during manual ranging
    try
        if isgraphics(peakMarker)
            delete(peakMarker);
        end
        markerTop = max(unrangedPeakHeights(ion)*1.1, spec.Parent.YLim(1)*10);
        markerTop = min(markerTop, spec.Parent.YLim(2));
        peakMarker = line(spec.Parent, ...
            [unrangedIonLocations(ion), unrangedIonLocations(ion)], ...
            [spec.Parent.YLim(1), markerTop], ...
            'Color', [0 0 0], 'LineStyle', '--', 'LineWidth', 1.5, ...
            'DisplayName', 'current peak');
    catch
    end
    
    % range the peak
    % Enter in rangeAdd (no picked limits) should skip this peak.
    try
        [rngH, ~] = rangeAdd(spec,colorScheme);
        if isempty(rngH) || ~isgraphics(rngH)
            fprintf('rangeAddAll: skipped peak at %.5f Da.\n', unrangedIonLocations(ion));
            continue;
        end
    catch ME
        if isSkipRangeAddError(ME)
            fprintf('rangeAddAll: skipped peak at %.5f Da.\n', unrangedIonLocations(ion));
            continue;
        end
        rethrow(ME);
    end
end

if isgraphics(peakMarker)
    delete(peakMarker);
end
spec.Parent.XLim = xPlotLimits; %revert to old plot limits.
spec.Parent.YLim = yPlotLimits; %revert to old plot limits.
if isgraphics(limPlotHandle)
    delete(limPlotHandle); %deletes background line
end

end

function tf = isSkipRangeAddError(ME)
% Enter in rangeAdd can produce an "invalid number of range limits" error.
msg = lower(string(ME.message));
tf = contains(msg, "invalid number of range limits");
end

function bw = estimateSpecBinWidth(spec)
bw = NaN;
try
    x = double(spec.XData(:));
    dx = diff(x);
    dx = dx(isfinite(dx) & dx > 0);
    if ~isempty(dx)
        bw = median(dx);
    end
catch
end
if ~isfinite(bw) || bw <= 0
    bw = 1e-6;
end
end

function tf = isPeakInExistingRange(spec, peakX, tol)
tf = false;
try
    rangeTable = rangesExtractFromMassSpec(spec);
    if isempty(rangeTable)
        return;
    end
    tf = any((peakX >= (rangeTable.mcbegin - tol)) & (peakX <= (rangeTable.mcend + tol)));
catch
    tf = false;
end
end

function [xOut, yOut] = mergeNearbyPeaks(xIn, yIn, minSep)
xOut = [];
yOut = [];
if isempty(xIn) || isempty(yIn)
    return;
end

x = double(xIn(:));
y = double(yIn(:));
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if isempty(x)
    return;
end

[x, idx] = sort(x, 'ascend');
y = y(idx);

if minSep <= 0
    [xUnique, ~, g] = unique(x, 'stable');
    yMax = accumarray(g, y, [], @max);
    xOut = xUnique;
    yOut = yMax;
    return;
end

startIdx = 1;
for i = 2:numel(x)
    if (x(i) - x(i-1)) > minSep
        [xSel, ySel] = selectClusterRepresentative(x(startIdx:i-1), y(startIdx:i-1));
        xOut(end+1,1) = xSel; %#ok<AGROW>
        yOut(end+1,1) = ySel; %#ok<AGROW>
        startIdx = i;
    end
end
[xSel, ySel] = selectClusterRepresentative(x(startIdx:end), y(startIdx:end));
xOut(end+1,1) = xSel;
yOut(end+1,1) = ySel;
end

function [xSel, ySel] = selectClusterRepresentative(xCluster, yCluster)
[ySel, idxMax] = max(yCluster);
idxMatches = find(yCluster == ySel);
if numel(idxMatches) > 1
    [~, iMin] = min(xCluster(idxMatches));
    idxMax = idxMatches(iMin);
end
xSel = xCluster(idxMax);
end
