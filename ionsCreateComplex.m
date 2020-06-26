function ionList = ionsCreateComplex(elements,complexity,isotopeTable,chargeStates)
% ionsCreateComplex creates a list of all complex ions that can be formed
% by the given elements, based on the isotopeTable. chargeStates are
% optional, default is +1 - +3
%
% WARNING: the number of permutations of ions grows with the Gamma
% function. Beware when creating very large complexities
%
% elements is a cellStr in the form {'Fe','Cr','O',...} and can be directly
% extracted e.g. from the decomposed pos file, as a keyword, 'all' can be
% used. This uses all elements up to
%
% complexity is the complexites in vector form, e.g. [1 2 3]


if ~exist('chargeStates','var')
    chargeStates = [1, 2, 3];
end


%% create indivdual permutations of atomic combinations
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
    elementalPermutations = permn(elements,complexity(comp));
    numPerms = size(elementalPermutations,1);
    for elPerm = 1:numPerms
        ionPerm{end+1,1} = ionConvertName(transp(elementalPermutations(elPerm,:)),NaN,'plain',isotopeTable);
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

ionList = sortrows(ionList,3);