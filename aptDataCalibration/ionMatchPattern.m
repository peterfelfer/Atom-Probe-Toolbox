function [matchedIons, scores] = ionMatchPattern(peakTable, elements, isotopeTable, options)
% IONMATCHPATTERN Identify ions by matching isotopic patterns to detected peaks.
%
% [matchedIons, scores] = ionMatchPattern(peakTable, elements, isotopeTable)
%
% For each candidate ion (generated via ionsCreateComplex), computes the
% full isotopic pattern and correlates it with the observed peak heights.
% Ions whose isotopic fingerprint matches the spectrum pattern score high.
%
% This is more robust than simple closest-mc matching because it uses the
% characteristic shape (relative isotope abundances) to distinguish e.g.
% 110Pd+ from 104Pd2 12C++ — only the correct ion will have a pattern
% that matches across all its isotopes.
%
% INPUT
%   peakTable    - table from massSpecFindPeaks (needs: mc, height, prominence)
%   elements     - cell array of element symbols, e.g. {'Pd','H','C','N','O'}
%   isotopeTable - isotope reference table
%
% OPTIONS
%   'complexity'        - vector of complexities (default: [1 2 3])
%   'chargeStates'      - vector of charge states (default: [1 2])
%   'tolerance'         - mc matching tolerance in Da (default: 0.15)
%   'minPeakProminence' - min peak prominence (default: 100)
%   'minScore'          - min pattern match score to keep (default: 0.5)
%   'minAbundance'      - min isotope abundance fraction to include (default: 0.005)
%   'binWidth'          - spectrum bin width for height lookup (default: 0.01)
%   'showPlot'          - show diagnostic plot (default: false)
%
% OUTPUT
%   matchedIons - table: ionName, chargeState, score, complexity,
%                        peakPositions, theoreticalPattern, observedPattern
%   scores      - struct with all candidate scores for inspection
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
    options.tolerance (1,1) double {mustBePositive} = 0.15
    options.minPeakProminence (1,1) double = 100
    options.minScore (1,1) double = 0.5
    options.minAbundance (1,1) double = 0.005
    options.binWidth (1,1) double = 0.01
    options.showPlot (1,1) logical = false
end

peaks = peakTable(peakTable.prominence >= options.minPeakProminence, :);
if isempty(peaks)
    matchedIons = table();
    scores = struct();
    return;
end

% Generate candidate ions
fprintf('Generating candidate ions...\n');
ionList = ionsCreateComplex(elements, options.complexity, isotopeTable, options.chargeStates);
fprintf('Generated %d isotopic combinations.\n', height(ionList));

% Get unique ion species (not isotopic variants)
uniqueIonNames = unique(ionList.ion);
nUnique = numel(uniqueIonNames);
fprintf('Unique ion species: %d\n', nUnique);

% For each unique ion × charge state, compute pattern match score
resultIonName = {};
resultChargeState = [];
resultScore = [];
resultComplexity = [];
resultNIsotopes = [];
resultPeakPositions = {};
resultTheoPattern = {};
resultObsPattern = {};

for cs = options.chargeStates
    for u = 1:nUnique
        ionName = uniqueIonNames(u);

        % Get all isotopic variants for this ion
        mask = ionList.ion == ionName;
        mcValues = ionList.mc(mask);

        % Convert to mc for this charge state
        mcCS = mcValues / cs;

        % Get abundances via ionsCreateIsotopeList
        try
            [~, abundances, weights] = ionsCreateIsotopeList(char(ionName), isotopeTable);
        catch
            continue;
        end

        if isempty(abundances)
            continue;
        end

        % Filter by min abundance
        keepMask = abundances >= options.minAbundance;
        if sum(keepMask) < 1, continue; end

        theoWeights = weights(keepMask) / cs;
        theoAbund = abundances(keepMask);
        theoAbund = theoAbund / max(theoAbund); % normalise to max = 1

        % Match each theoretical isotope to a detected peak
        nIsotopes = numel(theoWeights);
        obsHeights = zeros(1, nIsotopes);
        matchedPeakMc = zeros(1, nIsotopes);
        nMatched = 0;

        for iso = 1:nIsotopes
            [minDist, pIdx] = min(abs(peaks.mc - theoWeights(iso)));
            if minDist <= options.tolerance
                obsHeights(iso) = peaks.height(pIdx);
                matchedPeakMc(iso) = peaks.mc(pIdx);
                nMatched = nMatched + 1;
            end
        end

        % Need at least 1 matched isotope
        if nMatched == 0, continue; end

        % Compute pattern match score
        % Normalise observed heights to max = 1
        if max(obsHeights) > 0
            obsNorm = obsHeights / max(obsHeights);
        else
            continue;
        end

        % Score: cosine similarity between theoretical and observed patterns
        dotProd = sum(theoAbund .* obsNorm);
        normTheo = sqrt(sum(theoAbund.^2));
        normObs = sqrt(sum(obsNorm.^2));

        if normTheo == 0 || normObs == 0
            cosineSim = 0;
        else
            cosineSim = dotProd / (normTheo * normObs);
        end

        % Penalty for unmatched isotopes (expected but not found)
        fractionMatched = nMatched / nIsotopes;

        % Confidence bonus: ions with more matched isotopes provide
        % stronger evidence. A 6-isotope match is far more informative
        % than a 1-isotope match.
        confidenceWeight = 1 - exp(-nMatched / 2);  % 0.39 for 1, 0.63 for 2, 0.95 for 6

        % Compute complexity (number of atoms in ion)
        try
            ionTable = ionConvertName(char(ionName));
            comp = height(ionTable);
        catch
            comp = 1;
        end

        % Complexity penalty: prefer simpler ions (fewer atoms)
        complexityPenalty = 1 / comp;  % 1.0 for atomic, 0.5 for diatomic, 0.33 for triatomic

        score = cosineSim * fractionMatched * confidenceWeight * complexityPenalty;

        % Build charge state string
        chargeStr = repmat('+', 1, cs);
        fullName = [char(ionName) ' ' chargeStr];

        resultIonName{end+1, 1} = fullName;
        resultChargeState(end+1, 1) = cs;
        resultScore(end+1, 1) = score;
        resultComplexity(end+1, 1) = comp;
        resultNIsotopes(end+1, 1) = nIsotopes;
        resultPeakPositions{end+1, 1} = matchedPeakMc;
        resultTheoPattern{end+1, 1} = theoAbund;
        resultObsPattern{end+1, 1} = obsNorm;
    end
end

if isempty(resultIonName)
    fprintf('No ions matched.\n');
    matchedIons = table();
    scores = struct();
    return;
end

matchedIons = table(resultIonName, resultChargeState, resultScore, ...
    resultComplexity, resultNIsotopes, ...
    resultPeakPositions, resultTheoPattern, resultObsPattern, ...
    'VariableNames', {'ionName', 'chargeState', 'score', 'complexity', ...
                      'nIsotopes', 'peakPositions', 'theoPattern', 'obsPattern'});

% Sort by score (descending), then by complexity (ascending) as tiebreaker
matchedIons = sortrows(matchedIons, {'score', 'complexity'}, {'descend', 'ascend'});

% Filter by minimum score
matchedIons = matchedIons(matchedIons.score >= options.minScore, :);

scores = struct('allScores', matchedIons);

fprintf('Ions with score >= %.2f: %d\n', options.minScore, height(matchedIons));

% Optional diagnostic plot
if options.showPlot && height(matchedIons) > 0
    nShow = min(20, height(matchedIons));
    figure('Name', 'Ion Pattern Match Scores', 'NumberTitle', 'off');
    barh(1:nShow, matchedIons.score(nShow:-1:1));
    set(gca, 'YTick', 1:nShow, 'YTickLabel', matchedIons.ionName(nShow:-1:1));
    xlabel('Pattern match score (cosine similarity)');
    title(sprintf('Top %d ion identifications by isotopic pattern', nShow));
    xlim([0 1.05]);
    grid on;
end

end
