
function OP = apl_overlapProblemBuilder(rangeTable, varargin)
% apl_overlapProblemBuilder builds overlap problem struct from rangeTable
%
% OP = apl_overlapProblemBuilder(rangeTable, 'isotopeTable', tbl, 'addIons', {...}, 'minAbundance', 0.21)

p = inputParser;
addParameter(p, 'isotopeTable', []);
addParameter(p, 'ions', []);
addParameter(p, 'ionTable', []);
addParameter(p, 'massSpec', []);
addParameter(p, 'addIons', {});
addParameter(p, 'minAbundance', 0.21);
addParameter(p, 'ciSamples', 300);
parse(p, varargin{:});

isotopeTable = apl_isotopeTableLoad(p.Results.isotopeTable);
minAb = p.Results.minAbundance;

if isempty(rangeTable) || ~istable(rangeTable)
    error('apl_overlapProblemBuilder:invalidRangeTable', 'rangeTable must be a table.');
end

% Ion list
if ~isempty(p.Results.massSpec)
    ionTable = ionsExtractFromMassSpec(p.Results.massSpec);
    ions = apl_ionListFromIonTable(ionTable, p.Results.addIons);
elseif ~isempty(p.Results.ionTable)
    ions = apl_ionListFromIonTable(p.Results.ionTable, p.Results.addIons);
elseif isempty(p.Results.ions)
    ions = apl_ionListFromRangeTable(rangeTable, p.Results.addIons);
else
    if istable(p.Results.ions)
        ions = p.Results.ions;
    else
        ions = apl_ionListFromRangeTable(table(), p.Results.ions);
    end
end

% Build fraction matrix
[F, rangeCenters] = apl_buildFractionMatrix(rangeTable, ions, isotopeTable, minAb);

% Overlap groups
[groupRanges, groupIons, rangeToGroup] = apl_findOverlapGroups(F);

% Build overlap tables
overlapTable = cell(1, numel(groupRanges));
overlapIonNames = cell(1, numel(groupRanges));
rankDefFlag = false(1, numel(groupRanges));

for g = 1:numel(groupRanges)
    rIdx = groupRanges{g};
    iIdx = groupIons{g};
    overlapTable{g} = [rangeCenters(rIdx), F(rIdx, iIdx)];
    overlapIonNames{g} = cellstr(ions.ionLabel(iIdx));
    C = F(rIdx, iIdx);
    rankDefFlag(g) = (rank(C) < numel(iIdx));
end

% range elements from ion tables
rangeElements = {};
for i = 1:height(ions)
    ionTbl = ions.ionTable{i};
    if istable(ionTbl) && any(strcmp('element', ionTbl.Properties.VariableNames))
        rangeElements = [rangeElements; cellstr(categories(ionTbl.element))]; %#ok<AGROW>
    end
end
rangeElements = unique(rangeElements, 'stable');

OP = struct();
OP.ions = ions;
OP.ionNames = cellstr(ions.ionLabel);
OP.rangeTable = rangeTable;
OP.rangeElements = rangeElements;
OP.rangeNames = rangeTable.rangeName;
OP.rangeNamedIons = ismember(lower(string(rangeTable.rangeName)), ["background","not assigned","unknown"]);
OP.abundanceTable = [rangeCenters, F];
OP.fractionMatrix = F;
OP.overlapIons = groupIons;
OP.overlapTable = overlapTable;
OP.overlapRangeIdx = groupRanges;
OP.overlapIonNames = overlapIonNames;
OP.rangeToGroup = rangeToGroup;
OP.rankDefFlag = rankDefFlag;
OP.ciSamples = p.Results.ciSamples;
end
