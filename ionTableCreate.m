function list = ionTableCreate(elements,chargeStates,isotopeTable,maxComplexity,complexFormers,abundanceThreshold)
% ionTableCreate creates a table of all possible ions up to a max 
% complexity maxComplexity of elements (can be a list of atomic numbers, 
% symbols or 'all' for a list of chargeStates e.g. [1,2,3]. Complex formers 
% can be specified.
% 
% list = ionTableCreate(elements,chargeStates,isotopeTable,maxComplexity,complexFormers,abundanceThreshold)
%
% INPUTS
% elements:             can be a list of atomic numbers (cell array), 
%                       symbols (string array) or 'all'

% chargeStates:         list of chargeStates e.g. [1,2,3]    

% isotopeTable:         isotopeTable as given in the toolbox

% maxComplexity:        maximal complexity of the created ions

% complexFormers:       if complexformers = 'std' complexformers = H, H2, H3, C, C2, C3,
%                       N, O, O2. If complexFormers = 'and C P ....' it takes the elements
%                       in elements and adds the additional specified elements

% abundanceThreshold:   sets an threshold for used isotopes to build the
%                       list
% 
% OUTPUTS
% list:                 struct with .massToCharge (sorted), ionType{},
%                       relativeAbundance()

%% input interpretation (character inputs ==> list of atomic numbers)
% elements are parsed as a character array or string (whitespace delimited, 'Fe Ni Mo ...')
if isstring(elements)
    elements = char(elements);
end

if strcmp(elements,'all')
    % suggest all elements
    elements = unique(isotopeTable.element);
else
    elements = strread(elements,'%s');
    for i=1:length(elements)
        elements{i} = symbolConvertAtomicNumber(elements{i});
    end
    elements = cell2mat(elements);
end
elements = sort(elements); % final sorted list of atomic numbers of elements given


% include complex formers and interpret keywords 'std', 'same', 'and'
if exist('complexFormers','var')
    % interpretation of keywords
    if strcmp(complexFormers,'std')
        complexFormers = 'H C N O';
    elseif strcmp(complexFormers ,'same')
        complexFormers = elements;
    elseif strcmp(complexFormers(1:3),'and')
        complexFormers = [elements complexFormers(4:end)];
    end
    
    % create list of elements used as complex formers
    complexFormers = strread(complexFormers,'%s');
    for i=1:length(complexFormers)
        complexFormers{i} = number4sym(complexFormers{i});
    end
    complexFormers = cell2mat(complexFormers);
    complexFormers = sort(complexFormers); % final sorted list of atomic numbers of elements given
end

% if ion complexity is not parsed
if ~exist('maxComplexity','var')
    maxComplexity = 1;
end


%% create list of all nucleides used
%nucleides = nucleideList;

if ~strcmp(elements,'all')
    nucleides = [];
    for j = 1:length(elements)
        nucleides = [nucleides; isotopeTable([isotopeTable.atomicNumber==elements(j,1)],:)]; 
    end
end 


%% create list of complex ions

if exist('maxComplexity','var')

if exist('complexFormers','var')

    complexNucleides = isotopeTable;
    complexNucleidesUsed = [];
    for element = 1:length(complexFormers) 
            complexNucleidesUsed = [complexNucleidesUsed; complexNucleides([complexNucleides.atomicNumber == complexFormers(element,1)],:)];
    end
    
        complexNucleides = complexNucleidesUsed;
        
        
    else
        
        complexNucleides = nucleides;
        
        
end    
    %% based on the list 'complexNucleides', the complex ions are
    % calculated this list includes the element list and the list of
    % complex forming ions. This mean ions are formed by nucleides from
    % the list 'nucleidesUsed' + permutations of the nucleides from the
    % list 'complexNucleides' up to maxComplexity -1.
    
    tempIon = {};
    tempMass =[];
    tempAbundance = [];
    
    
    
    for comp = 1:maxComplexity
        tempComplex = createComplexIons(nucleides, complexNucleides, comp);
        tempIon = [tempIon; tempComplex.ionSpecies];
        tempMass = [tempMass; tempComplex.mass];
        tempAbundance = [tempAbundance; tempComplex.relativeAbundance];
        
    end
    
    uniqueIons = zeros(length(tempIon),maxComplexity);
    for ion = 1:length(tempIon)
        uniqueIons(ion,1:length(tempIon{ion}(:,1))) = tempIon{ion}(:,1)';
    end
    
    [uniqueIons idx idxOld] = unique(uniqueIons,'rows');
    
    tempIonID = idxOld;
    
end







%% produce output struct with chargestates list.ionSpecies,
% list.massToCharge, list.relativeAbundance, list.chargeState. list.ionID
% has a unique number for ions of the same type across the chargestates.
list.ionSpecies = {};
list.massToCharge = [];
list.relativeAbundance = [];
list.chargeState = [];
list.ionID = [];

for state = 1:length(chargeStates)
    list.ionSpecies = [list.ionSpecies; tempIon];
    list.massToCharge = [list.massToCharge; tempMass/chargeStates(state)];
    list.relativeAbundance = [list.relativeAbundance; tempAbundance];
    list.chargeState = [list.chargeState; repmat(chargeStates(state),length(tempIon),1)];
    list.ionID = [list.ionID; tempIonID];
    
    
end

% sorted according to mass to charge ratio
[list.massToCharge, idx] = sort(list.massToCharge);
list.ionSpecies = list.ionSpecies(idx);
list.relativeAbundance = list.relativeAbundance(idx);
list.chargeState = list.chargeState(idx);
list.ionID = list.ionID(idx);


%% if a threshold for the ionic abundance is specified
if exist('abundanceThreshold','var')
    if ischar(abundanceThreshold)
        abundanceThreshold = str2num(abundanceThreshold);
    end
    
    idx = list.relativeAbundance >= abundanceThreshold;
    list.ionSpecies = list.ionSpecies(idx);
    list.massToCharge = list.massToCharge(idx);
    list.relativeAbundance = list.relativeAbundance(idx);
    list.chargeState = list.chargeState(idx);
    list.ionID = list.ionID(idx);
end


end









