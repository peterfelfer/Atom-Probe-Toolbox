function [matchedIons, ionList] = ionMatchPeaks(peakTable, elements, isotopeTable, options)
% IONMATCHPEAKS Match candidate ions to detected mass spectrum peaks.
%
% [matchedIons, ionList] = ionMatchPeaks(peakTable, elements, isotopeTable)
%
% Generates all plausible ions from the given elements using
% ionsCreateComplex, then matches the most abundant isotope of each ion
% to detected peaks.
%
% INPUT
%   peakTable    - table from massSpecFindPeaks (needs: mc, height, prominence)
%   elements     - cell array of element symbols, e.g. {'Pd','H','C','N','O'}
%   isotopeTable - isotope reference table
%
% OPTIONS
%   'complexity'          - vector of complexities (default: [1 2 3])
%   'chargeStates'        - vector of charge states (default: [1 2])
%   'tolerance'           - mc matching tolerance in Da (default: 0.3)
%   'minPeakProminence'   - min peak prominence to consider (default: 100)
%   'minAbundance'        - min relative abundance to keep isotope (default: 0.01)
%
% OUTPUT
%   matchedIons - table of ions matched to peaks, sorted by peakHeight
%   ionList     - full ion list from ionsCreateComplex (before matching)
%
% USES: ionsCreateComplex, ionsCreateIsotopeList, ionConvertName
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    peakTable table
    elements cell
    isotopeTable table
    options.complexity (1,:) double {mustBePositive, mustBeInteger} = [1 2 3]
    options.chargeStates (1,:) double {mustBePositive, mustBeInteger} = [1 2]
    options.tolerance (1,1) double {mustBePositive} = 0.3
    options.minPeakProminence (1,1) double = 100
    options.minAbundance (1,1) double = 0.01
end

% Filter peaks
peaks = peakTable(peakTable.prominence >= options.minPeakProminence, :);
if isempty(peaks)
    warning('ionMatchPeaks:noPeaks', 'No peaks above prominence threshold.');
    matchedIons = table();
    ionList = table();
    return;
end

% Generate all candidate ions using the toolbox function
fprintf('Generating candidate ions (complexity %s, charge states %s)...\n', ...
    mat2str(options.complexity), mat2str(options.chargeStates));
ionList = ionsCreateComplex(elements, options.complexity, isotopeTable, options.chargeStates);
fprintf('Generated %d isotopic ion candidates.\n', height(ionList));

% For matching, we need the most abundant isotope per ion+chargeState.
% ionList has columns: ion (categorical), ionIsotopic (categorical), mc (double)
% Group by ion name and find the entry closest to integer mc (most abundant)
uniqueIons = unique(ionList.ionIsotopic);

% Match each ion candidate against detected peaks
candIonName = {};
candMc = [];
candPeakMc = [];
candPeakHeight = [];
candPeakProminence = [];

for i = 1:height(ionList)
    mc_cand = ionList.mc(i);

    % Find closest detected peak
    [minDist, pIdx] = min(abs(peaks.mc - mc_cand));
    if minDist <= options.tolerance
        candIonName{end+1, 1} = char(ionList.ionIsotopic(i));
        candMc(end+1, 1) = mc_cand;
        candPeakMc(end+1, 1) = peaks.mc(pIdx);
        candPeakHeight(end+1, 1) = peaks.height(pIdx);
        candPeakProminence(end+1, 1) = peaks.prominence(pIdx);
    end
end

if isempty(candIonName)
    fprintf('No ions matched any detected peaks.\n');
    matchedIons = table();
    return;
end

matchedIons = table(candIonName, candMc, candPeakMc, candPeakHeight, candPeakProminence, ...
    'VariableNames', {'ionName', 'theoreticalMc', 'peakMc', 'peakHeight', 'peakProminence'});

% Sort by peak height
matchedIons = sortrows(matchedIons, 'peakHeight', 'descend');

% For each detected peak, mark the best match (closest theoretical mc)
matchedIons.isBestMatch = false(height(matchedIons), 1);
uniquePeaks = unique(matchedIons.peakMc);
for p = 1:numel(uniquePeaks)
    groupMask = matchedIons.peakMc == uniquePeaks(p);
    groupIdx = find(groupMask);
    [~, bestInGroup] = min(abs(matchedIons.theoreticalMc(groupIdx) - uniquePeaks(p)));
    matchedIons.isBestMatch(groupIdx(bestInGroup)) = true;
end

nMatched = sum(matchedIons.isBestMatch);
fprintf('Matched %d ions to %d peaks (total candidates: %d).\n', ...
    height(matchedIons), nMatched, height(ionList));

end
