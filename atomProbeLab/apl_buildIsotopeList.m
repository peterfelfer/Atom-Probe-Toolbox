
function [ionType, abundance, weight] = apl_buildIsotopeList(ionTable, isotopeTable, minAbundance)
% apl_buildIsotopeList builds isotopic combinations for an ion.
% Respects explicitly specified isotopes in ionTable.
%
% ionTable: table with columns element (categorical) and isotope (numeric)
% isotopeTable: natural isotope table
% minAbundance: minimum abundance in percent (optional)

if nargin < 3 || isempty(minAbundance)
    minAbundance = 0.0; % percent
end

cols = apl_isotopeTableColumns(isotopeTable);
% normalize element comparison
isoElemStr = string(isotopeTable.(cols.element));


% normalize abundances to fraction
abund = isotopeTable.(cols.abundance);
if any(abund > 1)
    isotopeTable.(cols.abundance) = abund / 100;
end

numAtoms = height(ionTable);
if numAtoms == 0
    ionType = {};
    abundance = [];
    weight = [];
    return;
end

% Build isotope options per atom
atomIsos = cell(numAtoms,1);
atomAbund = cell(numAtoms,1);
atomMass = cell(numAtoms,1);

for a = 1:numAtoms
    elem = ionTable.element(a);
    iso = ionTable.isotope(a);
    elemMask = strcmpi(isoElemStr, string(elem));
    if ~any(elemMask)
        error('apl_buildIsotopeList:missingElement', 'Element %s not found in isotope table.', char(elem));
    end
    elemRows = isotopeTable(elemMask, :);
    if ~isnan(iso) && iso ~= 0
        isoMask = elemRows.(cols.isotope) == iso;
        if ~any(isoMask)
            % fallback: approximate mass as isotope number
            atomIsos{a} = iso;
            atomAbund{a} = 1;
            atomMass{a} = iso;
        else
            atomIsos{a} = elemRows.(cols.isotope)(isoMask);
            atomAbund{a} = elemRows.(cols.abundance)(isoMask);
            atomMass{a} = elemRows.(cols.mass)(isoMask);
        end
    else
        atomIsos{a} = elemRows.(cols.isotope);
        atomAbund{a} = elemRows.(cols.abundance);
        atomMass{a} = elemRows.(cols.mass);
    end
end

% Combine isotopes across atoms
comboMass = 0;
comboAbund = 1;
comboIsos = cell(1, numAtoms);
comboElems = cell(1, numAtoms);

for a = 1:numAtoms
    newMass = [];
    newAbund = [];
    newIsos = {};
    newElems = {};

    for c = 1:numel(comboAbund)
        for i = 1:numel(atomIsos{a})
            newMass(end+1,1) = comboMass(c) + atomMass{a}(i); %#ok<AGROW>
            newAbund(end+1,1) = comboAbund(c) * atomAbund{a}(i); %#ok<AGROW>
            tmpIsos = comboIsos(c,:);
            tmpElems = comboElems(c,:);
            tmpIsos{a} = atomIsos{a}(i);
            tmpElems{a} = char(ionTable.element(a));
            newIsos(end+1,:) = tmpIsos; %#ok<AGROW>
            newElems(end+1,:) = tmpElems; %#ok<AGROW>
        end
    end

    comboMass = newMass;
    comboAbund = newAbund;
    comboIsos = newIsos;
    comboElems = newElems;
end

% Apply abundance threshold
if minAbundance > 0
    mask = comboAbund >= (minAbundance / 100);
    comboMass = comboMass(mask);
    comboAbund = comboAbund(mask);
    comboIsos = comboIsos(mask,:);
    comboElems = comboElems(mask,:);
end

% Normalize abundances
if ~isempty(comboAbund)
    comboAbund = comboAbund / sum(comboAbund);
end

% Build ionType output
ionType = cell(1, numel(comboAbund));
for c = 1:numel(comboAbund)
    elem = categorical(comboElems(c,:)');
    iso = zeros(numAtoms,1);
    for a = 1:numAtoms
        iso(a) = comboIsos{c,a};
    end
    ionType{c} = table(elem, iso, 'VariableNames', {'element','isotope'});
end

abundance = comboAbund(:)';
weight = comboMass(:)';
end
