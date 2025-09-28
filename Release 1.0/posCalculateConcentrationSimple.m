function conc = posCalculateConcentrationSimple(pos,detEff,excludeList,volumeName)
% posCalculateConcentrationSimple calculates the concentration of a
% categorical list of atoms or ions
%
% conc = posCalculateConcentrationSimple(pos, detEff, excludeList,volumeName);
% conc = posCalculateConcentrationSimple(pos, detEff, excludeList)
% conc = posCalculateConcentrationSimple(pos, detEff)
%
% INPUT
% pos:          decomposed pos file that contains ion and charge state of
%               the individual atoms
%
% detEff:       detector efficiency of the atom probe, can be parsed as 
%               or as a fraction (for a LEAP 4000X HR it is 0.37)
%               
% excludeList:  cell array that contains as character the individual
%               ions that shall not be considered for the concentration 
%               calculation, unranged atoms appear as 'unranged', if not 
%               parsed, no atoms will be excluded
%
% volumeName:   name of the volume, parsed as character array, if not parsed,
%               the volume will be named after pos
%
% OUTPUT
% conc:         is a table that contains the count, concentration, and 
%               variance for each atom/ion that is not on the excludeList.
%               statistical deviation calculated after Danoix et al., 
%               https://doi.org/10.1016/j.ultramic.2007.02.005
%               variance(conc) = conc*(1-conc)/numAtomsDetected * (1 - detEff)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% detector efficiency
if detEff > 1
    detEff = detEff/100;
end


%% check if atomic or ionic concentration is calculated
if any(ismember(pos.Properties.VariableNames,'atom'))
    type = 'atomic';
    columnType = 'atom';
    atoms = pos.atom;
else
    type = 'ionic';
    columnType = 'ion';
    atoms = pos.ion;
end

if ~exist('volumeName','var')
    volumeName = inputname(1);
end

% need to assign, otherwise not counted in countcats
atoms(isundefined(atoms)) = 'unranged';

%% check for excluded types
cats = categories(atoms);
if exist('excludeList','var')
    isExcluded = ismember(cats,excludeList);
else
    isExcluded = false(size(cats));
end
isExcluded = isExcluded';

%% calculate concentrations for not excluded variables
counts = countcats(atoms);
if iscolumn(counts)
    counts = counts';
end
counts(2,:) = counts./sum(counts(~isExcluded));
counts(2,isExcluded) = 0;
counts(3,:) = counts(2,:).*(1- counts(2,:))./counts(1,:) * (1-detEff);



%% creating output table
conc = array2table(counts,'VariableNames',cats');
conc.Properties.VariableDescriptions = repmat({columnType},size(cats'));

conc = [table(categorical({volumeName;volumeName;volumeName}),'VariableNames',{'volume'})...
    table([0;0;0],'VariableNames',{'distance'} )...
    table(categorical({type;type;type}),'VariableNames',{'type'}),...
    table(categorical({'count'; 'concentration';'variance'}),'VariableNames',{'format'}),...
    conc];

