
function [posOut, rangeTableOut, info] = apl_localDecompAtom(pos, rangeTable, varargin)
% apl_localDecompAtom performs per-atom MLE deconvolution on overlaps
%
% OPTIONS (name-value):
%   'isotopeTable'   isotope table (table or .mat path)
%   'minAbundance'   minimum isotope abundance in percent (default 0.21)
%   'method'         'sphere' or 'knn' (default 'sphere')
%   'linkSize'       sphere radius [nm] or kNN count (default 2)
%   'maxmz'          max m/z for synthetic peaks (default max(rangeTable.mcend))
%   'seed'           rng seed (default [])
%   'maxAtoms'       limit number of overlapped atoms to process (default 0 = all)
%

p = inputParser;
addParameter(p, 'isotopeTable', []);
addParameter(p, 'minAbundance', 0.21);
addParameter(p, 'ionTable', []);
addParameter(p, 'massSpec', []);
addParameter(p, 'method', 'sphere');
addParameter(p, 'linkSize', 2);
addParameter(p, 'maxmz', []);
addParameter(p, 'seed', []);
addParameter(p, 'maxAtoms', 0);
parse(p, varargin{:});

if ~istable(pos)
    error('apl_localDecompAtom:invalidPos', 'pos must be a table');
end

if isempty(p.Results.seed)
    rng('shuffle');
else
    rng(p.Results.seed);
end

if isempty(p.Results.maxmz)
    maxmz = max(rangeTable.mcend);
else
    maxmz = p.Results.maxmz;
end

OP = apl_overlapProblemBuilder(rangeTable, 'isotopeTable', p.Results.isotopeTable, 'minAbundance', p.Results.minAbundance, 'ionTable', p.Results.ionTable, 'massSpec', p.Results.massSpec);
F = OP.fractionMatrix;

% Identify overlap ranges (more than one ion contributes)
overlapRanges = find(sum(F > 0, 2) > 1);

posOut = pos;
rangeIdx = apl_assignRangeIdx(pos.mc, rangeTable);

% mask for atoms in overlapping ranges
maskOverlap = ismember(rangeIdx, overlapRanges);
idxOverlap = find(maskOverlap);

if isempty(idxOverlap)
    rangeTableOut = rangeTable;
    info = struct('overlapProblem', OP, 'assigned', 0);
    return;
end

if p.Results.maxAtoms > 0 && numel(idxOverlap) > p.Results.maxAtoms
    idxOverlap = idxOverlap(1:p.Results.maxAtoms);
end

% Build KDTree on overlapped atoms only
X = [pos.x(maskOverlap), pos.y(maskOverlap), pos.z(maskOverlap)];
idxGlobal = find(maskOverlap);

Mdl = KDTreeSearcher(X);

method = lower(string(p.Results.method));
linkSize = p.Results.linkSize;

numAssign = numel(idxOverlap);
ionLabels = OP.ionNames;
numIons = numel(ionLabels);

% map range -> group and local index
rangeToGroup = OP.rangeToGroup;
rangeToLocal = cell(1, numel(OP.overlapRangeIdx));
for g = 1:numel(OP.overlapRangeIdx)
    rIdx = OP.overlapRangeIdx{g};
    map = zeros(height(rangeTable),1);
    map(rIdx) = 1:numel(rIdx);
    rangeToLocal{g} = map;
end

% prepare outputs
mcDeconv = pos.mc;
ionDeconv = strings(height(pos),1);
if ismember('ion', pos.Properties.VariableNames)
    try
        ionDeconv = string(pos.ion);
    catch
    end
end

% new m/z for each ion
newMz = maxmz + (1:numIons)';

r = rand(numAssign,1);

for k = 1:numAssign
    idx = idxOverlap(k);
    % position in KDTree index
    idxLocal = find(idxGlobal == idx, 1, 'first');
    if isempty(idxLocal)
        continue;
    end

    if method == "sphere"
        [nbrIdx,~] = rangesearch(Mdl, X(idxLocal,:), linkSize);
        nbrIdx = nbrIdx{1};
    else
        if linkSize < 1
            linkSize = 1;
        end
        [nbrIdx,~] = knnsearch(Mdl, X(idxLocal,:), 'K', linkSize+1);
    end

    % remove self
    nbrIdx = nbrIdx(nbrIdx ~= idxLocal);
    if isempty(nbrIdx)
        continue;
    end

    % determine overlap group
    rIdx = rangeIdx(idx);
    g = rangeToGroup(rIdx);
    if g == 0
        continue;
    end

    % map neighbor ranges to local indices
    neighborGlobalIdx = idxGlobal(nbrIdx);
    neighborRanges = rangeIdx(neighborGlobalIdx);
    localMap = rangeToLocal{g};
    localRanges = localMap(neighborRanges);
    localRanges = localRanges(localRanges > 0);
    if isempty(localRanges)
        continue;
    end

    numRanges = numel(OP.overlapRangeIdx{g});
    rangeCounts = accumarray(localRanges, 1, [numRanges, 1]);

    overlapTable = OP.overlapTable{g};
    [ionicRangeComp, failedIons] = apl_peakDecompMLE_ranges(rangeCounts, overlapTable, []);
    if any(failedIons)
        continue;
    end

    % local index of current range
    localIdx = localMap(rIdx);
    if localIdx == 0
        continue;
    end

    probs = ionicRangeComp(localIdx, :);
    if sum(probs) <= 0
        continue;
    end
    probs = probs / sum(probs);
    cumProbs = cumsum(probs);
    ionChoice = find(r(k) <= cumProbs, 1, 'first');
    if isempty(ionChoice)
        continue;
    end

    ionIdx = OP.overlapIons{g}(ionChoice);
    mcDeconv(idx) = newMz(ionIdx);
    ionDeconv(idx) = ionLabels{ionIdx};
end

posOut.mcDeconv = mcDeconv;
posOut.ionDeconv = ionDeconv;

% Build output rangeTable with extra deconv ranges
rangeTableOut = rangeTable;
numAdd = numIons;

% prepare add table with same variables
varNames = rangeTable.Properties.VariableNames;
add = table();
for v = 1:numel(varNames)
    add.(varNames{v}) = repmat(rangeTable.(varNames{v})(1,:), 0, 1); %#ok<AGROW>
end

addRows = table();
addRows.rangeName = categorical(string(ionLabels(:)));
addRows.chargeState = ionsChargeStates(OP.ions);
addRows.mcbegin = newMz - 0.1;
addRows.mcend = newMz + 0.1;
addRows.volume = zeros(numAdd,1);
addRows.ion = cell(numAdd,1);
for i = 1:numAdd
    addRows.ion{i} = OP.ions.ionTable{i};
end
addRows.color = repmat([0 0 0], numAdd, 1);

% align variables
for v = 1:numel(varNames)
    vn = varNames{v};
    if ~ismember(vn, addRows.Properties.VariableNames)
        % create missing column with matching type
        example = rangeTable.(vn);
        if isnumeric(example)
            addRows.(vn) = NaN(numAdd, size(example,2));
        elseif islogical(example)
            addRows.(vn) = false(numAdd, size(example,2));
        elseif iscategorical(example)
            addRows.(vn) = categorical(repmat(missing, numAdd, size(example,2)));
        elseif isstring(example)
            addRows.(vn) = strings(numAdd, size(example,2));
        elseif iscell(example)
            addRows.(vn) = cell(numAdd, size(example,2));
        else
            addRows.(vn) = repmat(example(1,:), numAdd, 1);
        end
    end
end
addRows = addRows(:, varNames);

rangeTableOut = [rangeTableOut; addRows];

info = struct();
info.overlapProblem = OP;
info.newMz = newMz;
info.assigned = numAssign;
end

function cs = ionsChargeStates(ions)
    cs = ions.chargeState;
end
