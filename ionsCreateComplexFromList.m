function ionList = ionsCreateComplexFromList(ions,chargeStates,isotopeTable)
% ionsCreateComplexFromList creates a list of all complex ions and their
% isotopic combinations from a list of known ions (without isotopic
% information). The information about isotope abundance and existance is
% provided via the isotopeTable
% The data for ions and chargeStates can be imported from the created, 
% maintained Excel file in the toolbox via the function ionsImportFromXLS
%
% ionList = ionsCreateComplexFromList(ion,chargeStates,isotopeTable)
%
% INPUT
% ions:         list of ions as a string vector or char cell array, e.g. 
%               ["Ni O"; "Fe2 O3";...] or {'Ni O'; 'Fe2 O3';...}
%
% chargeStates: N x M array of charge states the ion can be found in. M
%               must be the same as the length of the ions list. M is the
%               number of charge states. Entries are true or false for the
%               charge state
%
% isotopeTable: isotopeTable used for the analysis
%
%
% OUTPUT
% ionList:      table that contains the possible individual peaks. 
%               Table with fields:
%               ion, ionIsotopic: Categorical with ion names
%               mc: mass to charge state value of the individual isotopic
%               combination in amu
%
% TODO: permn creates error on large permutations (e.g. Si17)

%% create isotopic combinations from isotope table
ion = {};
ionIsotopic = {};
mc = [];

% include waitbar since this could take a while
wb = waitbar(0,'building isotopic combinations');

for i = 1:size(ions,1)
    ionIsotopicTmp = {};
    ionChargeStateTmp = find(chargeStates(i,:)); % charge states for the particular ion
    
    % create isotope combination list with weights for each elemental permutation
    [isoCombos, ~, weightTmp] = ionsCreateIsotopeList(ions(i), isotopeTable);
    
    for it = 1:length(isoCombos)
        ionIsotopicTmp{it,1} = ionConvertName(isoCombos{it});
    end
    
    ion = [ion; repmat(ions(i),[length(isoCombos) * length(ionChargeStateTmp) 1])];
    
    % create charge states
    ionIsotopicTmpCS = {};
    for cs = 1:length(ionChargeStateTmp)
        ionIsotopicTmpCS = [ionIsotopicTmpCS; cellfun(@(x) [x ' ' repmat('+', [1, ionChargeStateTmp(cs)])], ionIsotopicTmp,'UniformOutput',false)]; 
        mc = [mc; weightTmp'/ionChargeStateTmp(cs)];
    end
    
    ionIsotopic = [ionIsotopic; ionIsotopicTmpCS];
    
    % update waitbar
    waitbar(i/length(ions),wb);
end

delete(wb);


%% creating output table

ion = categorical(ion);
ionIsotopic = categorical(ionIsotopic);

ionList = table(ion,ionIsotopic,mc);
%ionList = unique(ionList); 

ionList = sortrows(ionList,3);

disp('fin');