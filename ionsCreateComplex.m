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
numElements = numel(elements);
elIdx = 1:numElements; % index of elements in the input vector. not the atomic number!

elementalPermutations = permn(elIdx,complexity);
numPerms = 




%% create isotopic combinations from isotope table



%% add charge atates in the mix