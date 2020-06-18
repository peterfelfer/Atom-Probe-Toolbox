function h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme)
% plots the relative abundances of an ion in a stem plot
% the information about the ion is given in h.UserData 


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

%change stem line depending on chargesate if only one charge state is given
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
