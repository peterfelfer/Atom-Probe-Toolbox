function ionList = createComplexIons(baseNucleides, complexNucleides, complexity)
% createcreates a list of complex ions with a certain complexity
%
% ionList = createComplexIons(baseNucleides, complexNucleides, complexity)
% 
% INPUT        
% spec:         mass spectrum to whom the stem plot is added to
% ion:          chemical symbol of the ion that will be added,
%               string scalar or string array
% chargeState:  charge state of the ion; if ion is a string array,
%               chargeState can be a scalar, a vector of charge states for
%               all ions (e.g.[1 2 3]) or a vector with the same number 
%               of entries as the ion list
% isotopeTable: the parsed isotope table is the basis of the relative
%               abundances
% colorScheme:  each ion has a different color
% sumMargin:    specifies a margin within which two peaks will be summed up
% minAbundance: is the minimal abundance, value between 0 and 1
% maxHeight:    is the height of the most abundant isotope (counts or
%               relative frequency)
%               default: to the YScale of the plotaxis
%               numeric value: you can type in the height of the highest
%               peak
%               'selection': uses a graphical input to select a peak, to
%               which the nearest isotopic combination will be scaled
%               'most abundant': adjusts the height of the most abundant 
%               peak to the closest peak in the mass spectrum
%               'least squares': adjusts the heights of the peaks to least 
%               squares match the peaks in the mass spectrum. Peaks that 
%               are close to other assigned peaks are not used.
% maxSeparation:is used when peak detection is used. The maximum of the
%               mass spectrum within this range will be used for scaling
% 
% OUTPUT
% h:            optional, is a 1x1 stem plot
% 
% THE FOLLOWING WILL BE STORED IN THE USER DATA SECTION OF THE PLOT:
% plotType = 'ion'
% isotopicCombinations: list of peaks vs. nucleides in the ion and
% charge state
%% if the complexity is 1, a list of both the base and complex nucleides is
%% returned with the ion species in a cell array

if complexity == 1
    
    ionList.ionSpecies = [baseNucleides(:,1:2); complexNucleides(:,1:2)];
    ionList.mass = [baseNucleides(:,3); complexNucleides(:,3)];
    ionList.relativeAbundance = [baseNucleides(:,4); complexNucleides(:,4)];
        
    %sorting after element and isotope
    
    [ionList.ionSpecies, idx] = sortrows(ionList.ionSpecies, [(-1) (-2)]);
    ionList.mass = ionList.mass(idx);
    ionList.relativeAbundance = ionList.relativeAbundance(idx);
    
    for ion = 1:length(ionList.mass)
        tempIonList{ion} = ionList.ionSpecies(ion,:);
    end
    
    ionList.ionSpecies = tempIonList';
    
    
    return;
    
end



%% check which ions are both base and complex formers and create complex ions of complex formers that are both base ions and
%% complex formers

bothNucleides = intersect(baseNucleides, complexNucleides, 'rows');

numBoth = length(bothNucleides(:,1));

if numBoth
    combinations = nreplacek(numBoth,complexity);
    
    ionSpeciesBoth = cell(length(combinations(:,1)),1);
    massBoth = zeros(length(combinations(:,1)),1);
    relativeAbundanceBoth = zeros(length(combinations(:,1)),1);
    
    
    for ion = 1:length(combinations(:,1))
        
        currentAbundance = 1;
        currentmass = 0;
        
        for atom = 1:(complexity)
            currentIon(atom,:) = bothNucleides(combinations(ion,atom),1:2);
            currentmass = currentmass + bothNucleides(combinations(ion,atom),3);
            currentAbundance = currentAbundance * (bothNucleides(combinations(ion,atom),4)/100);
        end
        
        ionSpeciesBoth{ion} = currentIon;
        massBoth(ion) = currentmass;
        relativeAbundanceBoth(ion) = currentAbundance;
        
    end
    clear currentIon
    
else
    ionSpeciesBoth = {};
    massBoth = [];
    relativeAbundanceBoth = [];
end



%% Create ions from complex formers combined with base ions
%% create permutations of complexNucleides (complexity-1)
complexity = complexity-1;

numComplex = length(complexNucleides(:,1));
numBase = length(baseNucleides(:,1));

% creates sample of unique ions
combinations = nreplacek(numComplex,complexity); 


ionSpecies = cell(length(combinations(:,1)),1);
mass = zeros(length(combinations(:,1)),1);
relativeAbundance = zeros(length(combinations(:,1)),1);


for ion = 1:length(combinations(:,1))
    
    currentAbundance = 1;
    currentmass = 0;
    
    for atom = 1:(complexity)
        currentIon(atom,:) = complexNucleides(combinations(ion,atom),1:2);
        currentmass = currentmass + complexNucleides(combinations(ion,atom),3);
        currentAbundance = currentAbundance * (complexNucleides(combinations(ion,atom),4)/100);
    end
    
    ionSpecies{ion} = currentIon;
    mass(ion) = currentmass;
    relativeAbundance(ion) = currentAbundance;
    
end







%% combine with base ions that are not complex formers

baseOnly = setdiff(baseNucleides,complexNucleides,'rows');
numBaseOnly = length(baseOnly(:,1));
numCombined = length(combinations(:,1)) * numBaseOnly;


ionSpeciesCombined = cell(numCombined,1);
massCombined = zeros(numCombined,1);
relativeAbundanceCombined = zeros(numCombined,1);

k = 1;
for base = 1:numBaseOnly
    
    for ion = 1:length(combinations(:,1))
        
        ionSpeciesCombined{k} = [baseOnly(base,1:2); ionSpecies{ion}];
        massCombined(k) = baseOnly(base,3) + mass(ion);
        relativeAbundanceCombined(k) = baseOnly(base,4)/100 * relativeAbundance(ion);
        
        k = k+1;
    end
        
        
    
end

clear ionSpecies mass relativeAbundance

ionSpecies = [ionSpeciesBoth; ionSpeciesCombined];
mass = [massBoth; massCombined];
relativeAbundance = [relativeAbundanceBoth; relativeAbundanceCombined];









%% renormalize isotopic abundances
numAll = length(mass);
complexity = complexity +1;
%sort ion species ascending according to first element then nucleide

tempIonList = zeros(numAll,(complexity)*2);
for ion = 1:numAll
    ionSpecies{ion} = sortrows(ionSpecies{ion},[-1]);
    tempIonList(ion,:) = [ionSpecies{ion}(:,1)' ionSpecies{ion}(:,2)'];
    
end

[tempIonList, idx] = sortrows(tempIonList,uminus([1:(complexity)*2]));

ionSpecies = ionSpecies(idx);
mass = mass(idx);
relativeAbundance = relativeAbundance(idx);





% find all unique elemental combinations of an ion
[ionType, idx, idxOld] = unique(tempIonList(:,1:(complexity)),'rows','first');
ionType = flipud(ionType);

multiplicity = zeros(length(idx),1);
for type = 1:length(idx)
    multiplicity(length(idx)-type+1) = sum(idxOld == idxOld(idx(type)));    
end



% correction for variations of ions with identical nucleides (e.g.
% 56Fe56Fe12C)
numNucleideVariations = zeros(length(idxOld),1);
for type = 1:length(idxOld)
    currentIon = ionSpecies{type};
    [nucleides nucleideIdx nucleideIdxOld] = unique(currentIon,'rows','first');
    
    nucleideMultiplicity = zeros(length(nucleideIdx(:,1)),1);
    for nuc = 1:length(nucleideIdx)
        nucleideMultiplicity(length(nucleideIdx)-nuc+1) = sum(nucleideIdxOld == nucleideIdxOld(nucleideIdx(nuc)));
        
    end
    
    numNucleideVariations(type) = factorial(complexity) / prod(factorial(nucleideMultiplicity));
    
end

relativeAbundance = relativeAbundance .* numNucleideVariations;


% renormalize 
normalization = zeros(length(idxOld(:,1)),1);
idx = flipud(idx);
for type = 1:length(idx)
    
    sumAbundance = sum(relativeAbundance(idx(type):(idx(type) + multiplicity(type) -1 )));
    
    normalization(idx(type):(idx(type) + multiplicity(type) -1 )) = sumAbundance;
    
end

relativeAbundance = relativeAbundance ./ normalization;


%% parse resulting ion list as a struct
ionList.ionSpecies = ionSpecies;
ionList.mass = mass;
ionList.relativeAbundance = relativeAbundance;
end
















