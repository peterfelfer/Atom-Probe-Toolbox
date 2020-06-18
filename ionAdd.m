
function h = ionAdd(spec,ion,chargeState,isotopeTable,colorScheme,sumMargin,minAbundance,maxHeight,maxSeparation)
% ionAdd creates a stem plot in the mass spectrum (spec) for the parsed ion
% and its corresponding charge state. The stem line will be according to
% charge state.
%
% ionAdd(spec,ion,chargeState,isotopeTable,colorScheme,sumMargin,minAbundance,maxHeight,maxSeparation)
%
% INPUT        
% spec:         mass spectrum to whom the stem plot is added to
% ion:          chemical symbol of the ion that will be added,
%               string scalar or string array
% chargeState:  charge state of the ion; if ion is a string array,
%               chargeState can be a scalar, a vector of charge states for
%               all ions (e.g.[1 2 3]) or a vector with the same number 
%               of entries as the ion list
% isotopeTable: the parsed isotope table is the basis of the relative
%               abundances
% colorScheme:  each ion has a different color
% sumMargin:    specifies a margin within which two peaks will be summed up
% minAbundance: is the minimal abundance, value between 0 and 1
% maxHeight:    is the height of the most abundant isotope (counts or
%               relative frequency)
%               default: to the YScale of the plotaxis
%               numeric value: you can type in the height of the highest
%               peak
%               'selection': uses a graphical input to select a peak, to
%               which the nearest isotopic combination will be scaled
%               'most abundant': adjusts the height of the most abundant 
%               peak to the closest peak in the mass spectrum
%               'least squares': adjusts the heights of the peaks to least 
%               squares match the peaks in the mass spectrum. Peaks that 
%               are close to other assigned peaks are not used.
% maxSeparation:is used when peak detection is used. The maximum of the
%               mass spectrum within this range will be used for scaling
% 
% OUTPUT
% h:            optional, is a 1x1 stem plot
% 
% THE FOLLOWING WILL BE STORED IN THE USER DATA SECTION OF THE PLOT:
% plotType = 'ion'
% isotopicCombinations: list of peaks vs. nucleides in the ion and
% charge state
% 
% ToDo:
%   - implement isnumeric('minHeight')
%   - implement 'least squares'
%   - implement checking for duplicate ions
% 

%% check ion name
if isstring(ion)
    ion = char(ion);
end

ax = spec.Parent;

%% get individual isotopic combinations
[ionType, abundance, weight] = ionsCreateIsotopeList(ion, isotopeTable);


%% cluster peaks together if sumMargin > 0
% if peaks are summed up, the dominant peak will be taken
if (sumMargin > 0) & (length(abundance) > 1)
    [ionType, abundance, weight] = ionsMergePeaks(ionType, abundance, weight, sumMargin);
end


%% eliminate peaks that are below minAbundance
weight(abundance <= minAbundance) = [];
ionType(abundance <= minAbundance) = [];
abundance(abundance <= minAbundance) = [];


%% create charge state ("CS")
plotWeight = [];
plotAbundance = [];
ionTypeCS = [];
ionCS = [];
for cs = 1:length(chargeState)
    plotWeight = [plotWeight, weight/chargeState(cs)];
    plotAbundance = [plotAbundance, abundance];
    ionTypeCS = [ionTypeCS ionType];
    ionCS = [ionCS repmat(chargeState(cs), 1, length(abundance))];
end


%% normalise height
if ~exist('maxHeight','var')
    % default is the height of current axis scaling
    maxPeak = max(abundance);
    maxDisp = ax.YLim(2);
    plotAbundance = plotAbundance * maxDisp / maxPeak;
elseif isnumeric(maxHeight)
    % if single value, most abundant isotope will be that height
    % if two element vector, the peak closest to maxHeight(1) will be
    % scaled to maxHeight(2)
    
    if isscalar(maxHeight)
        maxDisp = maxHeight;
        %XXXXX missing implementation
    else
        maxDisp = maxHeight(2);
        dist = maxHeight(1) - plotWeight;
        peakIdx = find(abs(dist) == min(abs(dist)));
        peakIdx = peakIdx(1);
        
        plotAbundance = plotAbundance * maxDisp / plotAbundance(peakIdx);
    end
    
    
    
else
    if strcmp(maxHeight,'selection')
        % selection of individual peak
        rect = getrect(ax);
        mcbegin = rect(1);
        mcend = rect(1) + rect(3);
        xdata = spec.XData;
        ydata = spec.YData;
        
        maxDisp = max(ydata(xdata > mcbegin & xdata < mcend));
        peak = max(plotAbundance((plotWeight > mcbegin) & (plotWeight < mcend)));
        
        if isempty(peak)
            error('peak not within selected range');
        end
        
        plotAbundance = plotAbundance * maxDisp / peak;
        
    elseif strcmp(maxHeight,'most abundant')
        % adjusts the height of the most abundant peak to the closest peak
        % in the mass spectrum
        mostAbundant = find(plotAbundance == max(plotAbundance)); % if multiple chargestates are present, there are as many max height peaks as charge states
        mostAbundant = mostAbundant(1);
        peak = plotAbundance(mostAbundant);
        
        mi = plotWeight(mostAbundant) - maxSeparation;
        mx = plotWeight(mostAbundant) + maxSeparation;
        
        inRange = spec.XData > mi & spec.XData < mx;
        maxDisp = max(spec.YData(inRange));
        
        
        plotAbundance = plotAbundance * maxDisp / peak;
        
    elseif strcmp(maxHeight,'least squares')
        % adjusts the heights of the peaks to least squares that match the
        % peaks in the mass spectrum. Peaks that are close to other 
        % assigned peaks are not used.
        error('least squares fitting not implemented yet');
    end
    
end



%% plot
h = ionStemPlot(ax, plotWeight, plotAbundance, ionTypeCS, chargeState, colorScheme);
