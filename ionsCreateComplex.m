function ionList = ionsCreateComplex(elements,complexity,isotopeTable,chargeStates)
% ionsCreateComplex creates a list of all complex ions that can be formed
% by the given elements, based on the isotopeTable. 
%
% ionList = ionsCreateComplex(elements,complexity,isotopeTable);
% ionList = ionsCreateComplex(elements,complexity,isotopeTable,chargeStates);
% 
% INPUT
% elements:     list of elements that form the complexes, either as cell,
%               string {'Fe','Cr','O',...} or as vector of atomic numbers
%               HINT: the cell string can be extracted from the decomposed
%               posfile as elements = categories(pos.atom)
%               HINT: to use many elements, standard matlab linear vector
%               generation can be used, e.g. elements = 1:100% 
%
% complexity:   vector of ion complexities to be included, e.g. [1 2 3]
%
% isotopeTable: isotopeTable used for the analysis
%
% chargeStates: vector of charge states to be included, e.g. [1 2 3], 
%               chargeStates are optional, default is +1 - +3
%
% OUTPUT
% ionList:      table with fields:
%               ion, ionIsotopic: categorical with ion names
%               mc: mass-to-charge value of the individual isotopic
%               combinations in amu
%
% WARNING: the number of permutations of ions grows with the Gamma
% function. Beware when creating very large complexities. For the entire
% periodic system, anything above a complexity of 2 will lead to a very
% long list.

%% looking for input chargeStates
if ~exist('chargeStates','var')
    chargeStates = [1, 2, 3];
end


%% create individual permutations of atomic combinations
% if input is cell string, convert to atomic numbers
if iscell(elements)
    for el = 1:numel(elements)
        elementsTmp(el) = symbolConvertAtomicNumber(elements{el});
    end
    elements = elementsTmp;
end

% for each complexity
ionPerm = {};
for comp = 1:length(complexity)
    % create permutations of given elements with repetitions
    elementalPermutations = permn(elements,complexity(comp));
    elementalPermutations = unique(sort(elementalPermutations,2),'rows');
    
    numPerms = size(elementalPermutations,1);
    for elPerm = 1:numPerms
        ionPerm{end+1,1} = ionConvertName((elementalPermutations(elPerm,:)'),NaN,'plain',isotopeTable);
    end
end


%% create isotopic combinations from isotope table
ion = {};
ionIsotopicNoCS = {};
weight = [];

% include waitbar since this could take a while
wb = waitbar(0,'building isotopic combinations');

for i = 1:length(ionPerm)
    ionIsotopicTmp = {};
    % create isotope combination list with weights for each elemental permutation
    [isoCombos, ~, weightTmp] = ionsCreateIsotopeList(ionPerm{i}, isotopeTable);
    for it = 1:length(isoCombos)
        ionIsotopicTmp{it,1} = ionConvertName(isoCombos{it});
    end
    ion = [ion; repmat(ionPerm(i),[length(isoCombos) 1])];
    ionIsotopicNoCS = [ionIsotopicNoCS; ionIsotopicTmp];
    weight = [weight; weightTmp'];
    
    % update waitbar
    waitbar(i/length(ionPerm),wb);
end

delete(wb);


%% add charge states in the mix
ion = repmat(ion,[length(chargeStates) 1]);
mc = [];
ionIsotopic = {};
for cs = 1:length(chargeStates)
    mc = [mc; weight/chargeStates(cs)];
    ionIsotopic = [ionIsotopic; cellfun(@(x) [x ' ' repmat('+', [1, chargeStates(cs)])], ionIsotopicNoCS,'UniformOutput',false)];
end

ion = categorical(ion);
ionIsotopic = categorical(ionIsotopic);

ionList = table(ion,ionIsotopic,mc);
%ionList = unique(ionList); 

ionList = sortrows(ionList,3);

