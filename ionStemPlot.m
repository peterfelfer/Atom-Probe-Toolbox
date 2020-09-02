function h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme)
% ionStemPlot plots the relative abundances of an ion in a stem plot
% the information about the ion is given in h.UserData 
% 
% h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme)
%
% INPUT
% ax:           axis of the current mass spectrum
%
% weight:       weight of the ion in amu
%
% adundance:    adundance for the chosen ion
%
% ionList:      list of all possible ions
% 
% chargeStates: charge states for the given ions
%
% colorScheme:  color scheme as provided or self made
%
% OUTPUT
% h:            handle to the stem plot


% generate ion name
ion = ionConvertName(ionList{1}.element);

h = stem(ax,weight, abundance);
try
    h.Color = colorScheme.color(colorScheme.ion == ion,:);
catch
    warning('ion color undefined');
end

if length(chargeStates) == 1
    h.DisplayName = [ion repmat('+',1,chargeStates)];
else
    h.DisplayName = ion;
end
h.LineWidth = 2;
h.UserData.plotType = "ion";
h.UserData.ion = ionList;
h.UserData.chargeState = chargeStates;
h.ButtonDownFcn = @(~,~) disp(h.DisplayName);

% change stem line depending on charge state if only one charge state is given
if length(chargeStates) == 1
    switch chargeStates
        case 1
            h.LineStyle = '--';
        case 2
            h.LineStyle = ':';
        case 3
            h.LineStyle = '-.';
        case 4
            h.LineStyle = '-';
            
    end
end
