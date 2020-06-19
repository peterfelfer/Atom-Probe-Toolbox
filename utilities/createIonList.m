function list = createIonList(elements,chargeStates,maxComplexity,complexFormers,abundanceThreshold)
% createIonList creates a list of all possible ions up to a max complexity
% maxComplexity of elements (can be a list of atomic numbers, symbols or
% 'all' for a list of chargeStates e.g. [1,2,3]. Complex formers can be
% specified.
% 
% list = ionTableCreate(elements,chargeStates,isotopeTable,maxComplexity,complexFormers,abundanceThreshold)
%
% INPUTS
% elements:             can be a list of atomic numbers, symbols or 'all'    
% 
% chargeStates:         list of chargeStates e.g. [1,2,3]    
% 
% isotopeTable:         isotopeTable as given in the toolbox
% 
% maxComplexity:        maximal complexity of the created ions
% 
% complexFormers:       if complexformers = 'std' complexformers = H, H2, 
%                       H3, He, B, C, C2, C3, N, O, O2, Ne. 
%                       If complexFormers = 'and C P ....' it takes the 
%                       elements in elements and adds the additional 
%                       specified elements
% 
% abundanceThreshold:   sets an threshold for used isotopes to biuld the
%                       list
% 
% OUTPUTS
% list:                 struct with .massToCharge (sorted), ionSpecies{},
%                       relativeAbundance(), chargeState(), ionID()

%% if complexformers = 'std' complexformers = H, H2, H3, He, B, C, C2, C3,
%% N, O, O2, Ne. If complexFormers = 'and C P ....' it takes the elements
%% in elements and adds the additional specified elements

if exist('complexFormers','var')
    
    if strcmp(complexFormers,'std')
        
        if ischar(elements)
            complexFormers = 'H C N O';
        else
            complexFormers = [1 2 6 7 8];
        end
    end
    
    if strcmp(complexFormers ,'same')
        complexFormers = elements;
    end
    
    
    if ischar(complexFormers)
        tempCF = strread(complexFormers,'%s');
        
        if strcmp(tempCF{1}, 'and')
            
            complexFormers = [elements complexFormers(4:end)];
        end
    end
end
    
    
    
    
    


%% elements can be parsed as a string (whitespace delimited, 'Fe Ni Mo ...') or atomic numbers
%% [26 28 42 ...]
if ischar(elements)
    
    if strcmp(elements,'all')
        
        % suggest all elements
        availableElements = nucleideList;
        elements = unique(availableElements(:,1));
        
    else
        
        elements = strread(elements,'%s');
        for i=1:length(elements)
            elements{i} = number4sym(elements{i});
        end
        
        elements = cell2mat(elements);
    end
    
end

elements = sort(elements);


%% chargeStates are parsed as either a vector [1 2 3] or a string('1 2 3').
if ischar(chargeStates)
    chargeStates = sort(str2num(chargeStates));
end
    


%% same parsing as for elements applies for complexFormers
if exist('complexFormers','var')
    
    if ischar(complexFormers)
        complexFormers = strread(complexFormers,'%s');
        for i=1:length(complexFormers)
            complexFormers{i} = symbolConvertAtomicNumber(complexFormers{i});
        end
        
        complexFormers = cell2mat(complexFormers);
        
    end
    
    complexFormers = sort(complexFormers);
    
end


if ~exist('maxComplexity','var')
    maxComplexity = 1;
end


if ischar(maxComplexity)
    maxComplexity = str2num(maxComplexity);
end








%% create list of all nucleides used
nucleides = nucleideList;

if ~strcmp(elements,'all')
    nucleidesUsed = [];
    for element = 1:length(elements)
        nucleidesUsed = [nucleidesUsed; nucleides(nucleides(:,1) == elements(element),:)];
    end
    
    nucleides = nucleidesUsed;
end









%% create list of complex ions
if exist('maxComplexity','var')
    
    if exist('complexFormers','var')
        
        
        
        complexNucleides = nucleideList;
        
        complexNucleidesUsed = [];
        for element = 1:length(complexFormers)
            complexNucleidesUsed = [complexNucleidesUsed; complexNucleides(complexNucleides(:,1) == complexFormers(element),:)];
        end
        
        complexNucleides = complexNucleidesUsed;
        
        
    else
        
        complexNucleides = nucleides;
        
        
    end
    
    
    
    %% based on the list 'complexNucleides', the complex ions are
    %% calculated this list includes the element list and the list of
    %% complex forming ions. This mean ions are formed by nucleides from
    %% the list 'nucleidesUsed' + permutations of the nucleides from the
    %% list 'complexNucleides' up to maxComplexity -1.
    
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
%% list.massToCharge, list.relativeAbundance, list.chargeState. list.ionID
%% has a unique number for ions of the same type across the chargestates.
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









