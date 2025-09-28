function chargeStateRatioTable = chargeStateRatioCalculate(element,pos,chargeStates)
% chargeStateRatioTable = chargeStateRatioCalculate(element,pos,chargeStates)
% calculates the chargeStateratio of an element and his chargeStates in the
% pos
%
% chargeStateRatioTable = chargeStateRatioCalculate(element,pos,chargeStates)
%
% INPUT
% 
% element:      character with the name of the element of interesst
%
% pos:          decomposed pos file that contains ion and charge state of
%               the individual atoms
%
% chargeStates: double withe the chargeStates of interest
%
% OUTPUT
% chargeStateRatioTable:    
%               table with total counts of chargeStates and the ratio to
%               the chargeStates
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

% total counts of elements with the right charge state
elementTotalNumber = size(pos(pos.ion == element,:),1);

chargeStatesOverviewDouble = cell(2,size(chargeStates,2));

% get peak ratios
for i = 1:size(chargeStates,2)
    chargeStateNumber = size(pos(pos.ion == element & ismember(pos.chargeState,chargeStates(i)),:),1);
    
    chargeStatesOverviewDouble{1,i} =  chargeStateNumber;
    chargeStatesOverviewDouble{2,i} =  chargeStateNumber./elementTotalNumber;
end

chargeStatesOverviewTable = array2table(chargeStatesOverviewDouble);

% save results in table
% variable names
newNames = append("ChargeState +",string(chargeStates));
chargeStatesOverviewTable = renamevars(chargeStatesOverviewTable,1:width(chargeStates),newNames);

% start tables with columns containing ion name and units
sz = [2 2];
varNames = ["Element","format"];
varTypes = ["string","string"];
chargeStatesRatio = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
chargeStatesRatio.Element = categorical({element; element;});
chargeStatesRatio.format = categorical({'count'; 'ratio'});

% add results (chargeStatesOverviewTable)
chargeStateRatioTable = [chargeStatesRatio chargeStatesOverviewTable];

end