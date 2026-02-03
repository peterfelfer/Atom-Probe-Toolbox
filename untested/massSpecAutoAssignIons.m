function [ionTable, info] = massSpecAutoAssignIons(peaks, elements, options)
% MASSSPECAUTOASSIGNIONS Auto-assign ions to detected mass spectrum peaks.
%
% ionTable = massSpecAutoAssignIons(peaks, elements)
% [ionTable, info] = massSpecAutoAssignIons(peaks, elements, ...)
%
% INPUT
% peaks:     table from massSpecFindPeaks (must contain variable 'mc' and
%            preferably 'prominence' or 'height')
% elements:  list of allowed elements (string/cell or atomic numbers)
%
% OPTIONS
%   'maxAtoms'        - max atoms per ion (default: 2)
%   'maxCharge'       - max charge state (default: 3)
%   'chargeStates'    - explicit charge state list (default: 1:maxCharge)
%   'toleranceMode'   - 'relative' or 'fixed' (default: 'relative')
%   'toleranceValue'  - value for tolerance (default: 1e-3)
%   'conservative'    - conservative assignment (default: true)
%   'minIsotopeMatches' - minimum isotopes to match (default: auto)
%   'minPatternScore' - minimum pattern score (default: auto)
%   'requireMajorPeak' - require most abundant isotope (default: true)
%   'isotopeTable'    - isotope table (default: isotopeTable_naturalAbundances)
%
% OUTPUT
% ionTable: table with columns ionName, chargeState, ion, color, isTracer
% info:     struct with candidate matches and settings
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    peaks table
    elements
    options.maxAtoms (1,1) double {mustBePositive, mustBeInteger} = 2
    options.maxCharge (1,1) double {mustBePositive, mustBeInteger} = 3
    options.chargeStates double = []
    options.toleranceMode (1,1) string = "relative"
    options.toleranceValue (1,1) double {mustBePositive} = 1e-3
    options.conservative (1,1) logical = true
    options.minIsotopeMatches (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    options.minPatternScore (1,1) double = NaN
    options.requireMajorPeak (1,1) logical = true
    options.isotopeTable = []
end

if ~istable(peaks) || ~ismember('mc', peaks.Properties.VariableNames)
    error('massSpecAutoAssignIons:invalidPeaks', ...
        'peaks must be a table from massSpecFindPeaks (with mc column).');
end

if isempty(options.chargeStates)
    chargeStates = 1:options.maxCharge;
else
    chargeStates = options.chargeStates(:)';
end

if isempty(options.isotopeTable)
    data = load('isotopeTable_naturalAbundances.mat', 'isotopeTable');
    isotopeTable = data.isotopeTable;
else
    isotopeTable = options.isotopeTable;
end

if options.minIsotopeMatches == 0
    if options.conservative
        minIsotopeMatches = 2;
    else
        minIsotopeMatches = 1;
    end
else
    minIsotopeMatches = options.minIsotopeMatches;
end

if isnan(options.minPatternScore)
    if options.conservative
        minPatternScore = 0.6;
    else
        minPatternScore = 0.3;
    end
else
    minPatternScore = options.minPatternScore;
end

baseIons = buildBaseIons(elements, options.maxAtoms, isotopeTable);

peakMc = peaks.mc(:);
if ismember('prominence', peaks.Properties.VariableNames)
    peakIntensity = peaks.prominence(:);
elseif ismember('height', peaks.Properties.VariableNames)
    peakIntensity = peaks.height(:);
else
    peakIntensity = ones(size(peakMc));
end

candidate = struct('ion', {}, 'chargeState', {}, 'score', {}, ...
    'matchFraction', {}, 'matchedIdx', {}, 'predictedMc', {}, ...
    'predictedAbund', {}, 'matchedMc', {}, 'matchedIntensity', {});

for i = 1:numel(baseIons)
    ion = baseIons{i};
    [~, abund, weight] = ionsCreateIsotopeList(ion, isotopeTable);
    if isempty(weight)
        continue;
    end
    for cs = chargeStates
        predMc = weight / cs;
        predAbund = abund(:);
        [matchIdx, score, matchFraction, matchedMc, matchedInt, matchedPredMask] = matchIsotopicPattern( ...
            predMc, predAbund, peakMc, peakIntensity, options, minIsotopeMatches, minPatternScore);

        if isempty(matchIdx)
            continue;
        end

        if options.requireMajorPeak
            [~, idxMajor] = max(predAbund);
            if idxMajor > numel(matchedPredMask) || ~matchedPredMask(idxMajor)
                continue;
            end
        end

        candidate(end+1).ion = ion; %#ok<AGROW>
        candidate(end).chargeState = cs;
        candidate(end).score = score;
        candidate(end).matchFraction = matchFraction;
        candidate(end).matchedIdx = matchIdx;
        candidate(end).predictedMc = predMc;
        candidate(end).predictedAbund = predAbund;
        candidate(end).matchedMc = matchedMc;
        candidate(end).matchedIntensity = matchedInt;
    end
end

if isempty(candidate)
    ionTable = table();
    info = struct('candidates', [], 'assigned', [], 'settings', options);
    return;
end

% Greedy assignment
scores = [candidate.score]';
[~, order] = sort(scores, 'descend');
assigned = candidate([]);
peakAssigned = false(numel(peakMc), 1);

for k = 1:numel(order)
    c = candidate(order(k));
    matchIdx = c.matchedIdx;
    newMask = ~peakAssigned(matchIdx);
    if options.conservative && sum(newMask) < min(1, minIsotopeMatches)
        continue;
    end
    if options.conservative && sum(newMask) == 0
        continue;
    end
    assigned(end+1) = c; %#ok<AGROW>
    peakAssigned(matchIdx) = true;
end

% Build ionTable for ionAdd
numAssign = numel(assigned);
ionName = strings(numAssign, 1);
chargeState = zeros(numAssign, 1);
ion = cell(numAssign, 1);
color = nan(numAssign, 3);
isTracer = false(numAssign, 1);

for i = 1:numAssign
    ion{i} = assigned(i).ion;
    ionName(i) = string(assigned(i).ion);
    chargeState(i) = assigned(i).chargeState;
end

ionName = categorical(ionName);
ionTable = table(ionName, chargeState, ion, color, isTracer);

info = struct();
info.candidates = candidate;
info.assigned = assigned;
info.settings = options;
end

function baseIons = buildBaseIons(elements, maxAtoms, isotopeTable)
    if isstring(elements)
        elements = cellstr(elements);
    end

    if iscell(elements)
        elemNums = zeros(numel(elements), 1);
        for i = 1:numel(elements)
            elem = elements{i};
            if isstring(elem)
                elem = char(elem);
            end
            if strcmpi(elem, 'Cb')
                elem = 'Nb';
            end
            elemNums(i) = symbolConvertAtomicNumber(elem);
        end
    else
        elemNums = elements(:);
    end

    baseList = {};
    for comp = 1:maxAtoms
        combos = permn(elemNums, comp);
        combos = unique(sort(combos, 2), 'rows');
        for k = 1:size(combos, 1)
            baseList{end+1, 1} = ionConvertName(combos(k, :)', NaN, 'plain', isotopeTable); %#ok<AGROW>
        end
    end

    baseIons = unique(baseList);
end

function [matchIdx, score, matchFraction, matchedMc, matchedInt, matchedPredMask] = matchIsotopicPattern( ...
        predMc, predAbund, peakMc, peakIntensity, options, minIsotopeMatches, minPatternScore)

    predMc = predMc(:);
    predAbund = predAbund(:);
    numPred = numel(predMc);
    matchedIdx = zeros(numPred, 1);
    matchedMc = nan(numPred, 1);
    matchedInt = zeros(numPred, 1);

    for i = 1:numPred
        if options.toleranceMode == "relative"
            tol = predMc(i) * options.toleranceValue;
        else
            tol = options.toleranceValue;
        end
        idx = find(abs(peakMc - predMc(i)) <= tol);
        if isempty(idx)
            continue;
        end
        [~, relIdx] = min(abs(peakMc(idx) - predMc(i)));
        pick = idx(relIdx);
        matchedIdx(i) = pick;
        matchedMc(i) = peakMc(pick);
        matchedInt(i) = peakIntensity(pick);
    end

    matchedMask = matchedIdx > 0;
    matchCount = sum(matchedMask);
    if matchCount == 0
        matchIdx = [];
        score = 0;
        matchFraction = 0;
        matchedMc = [];
        matchedInt = [];
        matchedPredMask = matchedMask;
        return;
    end

    matchFraction = matchCount / numPred;
    if matchCount < min(minIsotopeMatches, numPred)
        matchIdx = [];
        score = 0;
        matchedMc = [];
        matchedInt = [];
        matchedPredMask = matchedMask;
        return;
    end

    expNorm = predAbund / sum(predAbund);
    obs = matchedInt;
    obsNorm = obs / max(sum(obs), eps);
    tvd = 0.5 * sum(abs(obsNorm - expNorm));
    score = matchFraction * (1 - tvd);

    if score < minPatternScore
        matchIdx = [];
        matchedMc = [];
        matchedInt = [];
        matchedPredMask = matchedMask;
        return;
    end

    matchIdx = unique(matchedIdx(matchedMask));
    matchedMc = matchedMc(matchedMask);
    matchedInt = matchedInt(matchedMask);
    matchedPredMask = matchedMask;
end
