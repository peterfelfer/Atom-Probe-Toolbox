function lossTable = intensityLossCalculate(spec)
% intensityLossCalculate puts the intensity values together and calculates
% the intensity loss for each isotope in percentage
 
 
%% reorder of mass spec plot, probably not necessary
massSpecReorderPlot(spec);
 
%% a table of all ion stem plot y values
stemPlots = findobj('Type','Stem');
% get the y coordinates
% for i = 1:numel(stemPlots);
%     yCo = stemPlots(i).YData;
%     yCo = yCo';
%     onlyPart = yCo;
%     allPartsCombined = [onlyPart; yCo];
%     yCo = allPartsCombined;
%     yCo = table(yCo);
%     yCo.Properties.VariableNames{1} = 'y values of ion stem plots';
% end
allYCo = table();  % table with all the Y coordinates
for i = 1:numel(stemPlots)      % get the number of all stem plots
    yCo = stemPlots(i).YData;   % get all the Y coordinates of the corresponding stem plot (can be more than 1)
    yCo = yCo';
    yCo = table(yCo);           % create a table 
    allYCo = [allYCo; yCo] ; % add the new Y Coordinates from one stem plot with the y coordinates of the other stem plots 
end
allYCo.Properties.VariableNames{1} = 'y values of ion stem plots'; % name of the first row


 
%% y value of ion stem plot of the isotope and max y value of corresponding range
ionStem = findobj('Type','Stem');
numberOfIsotopes = numel(ionStem.YData);
for i = 1:numberOfIsotopes;
    heightIonStem = ionStem.YData(i);
        onlyArea = findobj('Type','Area');
    maxRangeValue = max(onlyArea(i).YData);
    lossPercentage = ((heightIonStem-maxRangeValue)/heightIonStem)*100;
    isotope = onlyArea(i).DisplayName;
    isotope = convertCharsToStrings(isotope);
    
    % creation of table with all values and resulting loss in percentage
    lossTable(i,:) = table(isotope,heightIonStem,maxRangeValue,lossPercentage);
    lossTable.Properties.VariableNames{1} = 'isotope';
    lossTable.Properties.VariableNames{2} = 'height of ion stem';
    lossTable.Properties.VariableNames{3} = 'maximum height of corresponding range';
    lossTable.Properties.VariableNames{4} = 'intensity loss [%]';
    
end
 
 
end


