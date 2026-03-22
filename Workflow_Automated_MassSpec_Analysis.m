%% Automated Mass Spectrum Analysis
% This workflow performs automated peak detection, ion identification via
% isotopic pattern matching, and ranging on a calibrated atom probe mass
% spectrum.
%
% Steps:
%   1. Load data and create mass spectrum
%   2. Peak detection with adjustable parameters
%   3. Specify elements and generate candidate ions
%   4. Match candidates to detected peaks via isotopic pattern matching
%   5. Review and edit ion assignments
%   6. Create ranges using auto EER algorithm
%
% Prerequisites: calibrated pos table with mc column (from .epos, .pos,
% or pyccapt data after tofToMassToCharge + voltage/bowl correction).

%% 1. Setup and Load Data
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme.mat');

% Load dataset — edit path or use file dialog
pos = posLoad();

%% 2. Peak Detection
% Create initial mass spectrum and detect peaks.
% Adjust the parameters below to control sensitivity.

% --- Peak detection parameters (adjust these) ---
binWidth            = 0.01;   % Da — histogram bin width
minProminenceFactor = 6;      % prominence must exceed factor x local noise
minDistanceDa       = 0.3;    % Da — minimum distance between peaks
baselineQuantile    = 0.1;    % quantile for baseline estimation
smoothSpan          = 0.1;    % Da — smoothing window for residual

% Run peak detection
[peakTable, peakInfo] = massSpecFindPeaks(pos.mc, ...
    'binWidth', binWidth, ...
    'minProminenceFactor', minProminenceFactor, ...
    'minDistanceDa', minDistanceDa, ...
    'baselineQuantile', baselineQuantile, ...
    'smoothSpan', smoothSpan);

fprintf('Detected %d peaks.\n', height(peakTable));

%% 2b. Diagnostic Plot — Inspect Peak Detection
% Visualise the spectrum, baseline, noise floor, and detected peaks.
% If too many or too few peaks, adjust the parameters above and re-run.

% --- Minimum prominence filter for display ---
minProminenceDisplay = 500;

sigPeaks = peakTable(peakTable.prominence > minProminenceDisplay, :);
fprintf('Significant peaks (prominence > %d): %d\n', minProminenceDisplay, height(sigPeaks));

figure('Name', 'Peak Detection', 'NumberTitle', 'off', 'Position', [100 100 1200 500]);
semilogy(peakInfo.centers, peakInfo.counts, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
hold on;
semilogy(peakInfo.centers, peakInfo.baseline, 'b-', 'LineWidth', 1.5);
semilogy(peakInfo.centers, peakInfo.baseline + peakInfo.noise, 'g--', 'LineWidth', 1);
semilogy(sigPeaks.mc, sigPeaks.height, 'rv', 'MarkerSize', 5, 'MarkerFaceColor', 'r');

for i = 1:height(sigPeaks)
    if sigPeaks.height(i) > 3000
        text(sigPeaks.mc(i), sigPeaks.height(i) * 1.8, sprintf('%.1f', sigPeaks.mc(i)), ...
            'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', 'r');
    end
end

xlim([0 min(200, max(pos.mc))]);
ylim([1 max(peakInfo.counts) * 3]);
xlabel('Mass-to-charge (Da)');
ylabel('Counts per bin');
legend('Spectrum', 'Baseline', 'Noise threshold', 'Detected peaks', 'Location', 'northeast');
title(sprintf('Peak detection: %d significant peaks', height(sigPeaks)));
hold off;

%% 2c. Peak Table — Review Detected Peaks
disp(sigPeaks(:, {'mc', 'height', 'prominence'}));

%% 3. Specify Elements and Ion Generation Parameters
% List the chemical elements present in your specimen.
% The algorithm generates all plausible ionic and molecular species
% from these elements (via ionsCreateComplex), computes their isotopic
% fingerprints, and matches the patterns against the detected peaks.

% --- Edit these for your specimen ---
elements = {'Pd', 'H', 'C', 'N', 'O'};

% --- Ion generation parameters ---
complexity     = [1 2 3];  % atom counts per ion (1=atomic, 2=diatomic, 3=triatomic)
chargeStates   = [1 2];    % charge states to consider
matchTolerance = 0.15;     % Da — how close each isotope must be to a peak

%% 4. Isotopic Pattern Matching
% For each candidate ion, the full isotopic pattern (all isotope positions
% and relative abundances) is compared against the observed peak heights.
% Ions with many isotopes matching the correct ratios score highest.
%
% Scoring:
%   score = cosine_similarity x fraction_matched x confidence x (1/complexity)
%
% Cosine similarity measures how well relative peak heights match the
% expected natural abundance ratios. Confidence rewards ions with more
% matched isotopes (a 6-isotope match like Pd is far more informative
% than a 1-isotope match like H). Complexity penalises molecular ions.

[matchedIons, ~] = ionMatchPattern(peakTable, elements, isotopeTable, ...
    'complexity', complexity, ...
    'chargeStates', chargeStates, ...
    'tolerance', matchTolerance, ...
    'minPeakProminence', minProminenceDisplay, ...
    'minScore', 0.05, ...
    'showPlot', true);

fprintf('\nTop ion identifications:\n');
fprintf('%-25s  CS  Score  Cmpl  nIso\n', 'Ion');
nShow = min(30, height(matchedIons));
for i = 1:nShow
    fprintf('%-25s  %d   %.3f   %d     %d\n', ...
        matchedIons.ionName{i}, matchedIons.chargeState(i), ...
        matchedIons.score(i), matchedIons.complexity(i), ...
        matchedIons.nIsotopes(i));
end

%% 4b. Select Ions to Keep
% Review the ranked list above and choose a score threshold.
% Ions below the threshold are discarded.
% You can also manually add or remove ions after this step.

% --- Adjust this threshold ---
scoreThreshold = 0.15;

selectedIons = matchedIons(matchedIons.score >= scoreThreshold, :);
fprintf('\nSelected %d ions with score >= %.2f:\n', height(selectedIons), scoreThreshold);
disp(selectedIons(:, {'ionName', 'score', 'chargeState', 'complexity', 'nIsotopes'}));

%% 4c. Manual Edits (optional)
% Remove incorrect assignments:
%   selectedIons(strcmp(selectedIons.ionName, 'O N ++'), :) = [];
%
% The selected ions will be used for ranging in the next step.

%% 5. Create Mass Spectrum with Ion Stems
% Plot the calibrated mass spectrum and overlay ion stems for the
% selected ions.

spec = massSpecPlot(pos.mc, binWidth, 'normalised');
xlim([0 min(150, max(pos.mc))]);

% Add stems for each unique element at each charge state
stemmed = struct();
for i = 1:height(selectedIons)
    ionName = selectedIons.ionName{i};
    cs = selectedIons.chargeState(i);

    % Extract the base element name (first word before space or +)
    tokens = regexp(ionName, '(\S+)\s*\++$', 'tokens');
    if isempty(tokens), continue; end
    formula = tokens{1}{1};

    % Only add stems for single-element ions (ionAdd handles the isotopes)
    formulaParts = strsplit(strtrim(formula));
    if numel(formulaParts) > 1, continue; end

    % Check it's a pure element (no digits = no isotope prefix)
    elemName = regexprep(formula, '^\d+', '');
    key = matlab.lang.makeValidName(sprintf('%s_%d', elemName, cs));

    if ~isfield(stemmed, key)
        try
            ionAdd(spec, elemName, cs, isotopeTable, colorScheme, 0, 0.005, 'most abundant', 0.5);
            stemmed.(key) = true;
        catch
        end
    end
end

title('Mass spectrum with matched ion stems');

%% 6. Auto-Range with EER Algorithm
% Collect all peak mc positions from the selected ions and compute
% optimal range boundaries using the Equal Error Rate criterion.

% Gather peak positions from each ion's isotope pattern
allPeakMc = [];
for i = 1:height(selectedIons)
    positions = selectedIons.peakPositions{i};
    positions = positions(positions > 0);  % remove unmatched zeros
    allPeakMc = [allPeakMc, positions]; %#ok<AGROW>
end
allPeakMc = unique(allPeakMc);

fprintf('Ranging %d unique peak positions.\n', numel(allPeakMc));

[eerRanges, eerInfo] = rangeAutoEER(pos.mc, allPeakMc, ...
    'binWidth', binWidth, 'showPlot', true);

fprintf('EER ranges computed: %d\n', height(eerRanges));

%% 7. Apply Ranges to Mass Spectrum
% For each EER range, find the best-matching ion and apply it with rangeAdd.
% Molecular ions are added to the colorScheme automatically.

% Pre-add molecular ions to colorScheme
for i = 1:height(selectedIons)
    ionName = selectedIons.ionName{i};
    tokens = regexp(ionName, '(.+?)\s*\++$', 'tokens');
    if isempty(tokens), continue; end
    formula = strtrim(tokens{1}{1});
    formulaParts = strsplit(formula);
    if numel(formulaParts) > 1
        catName = categorical(string(formula));
        if ~any(colorScheme.ion == catName)
            newColor = rand(1, 3) * 0.5 + 0.25;
            colorScheme = [colorScheme; table(catName, newColor, 'VariableNames', {'ion', 'color'})]; %#ok<AGROW>
        end
    end
end

% Build a lookup: for each peak mc, find the best ion name
peakToIon = containers.Map('KeyType', 'double', 'ValueType', 'char');
for i = 1:height(selectedIons)
    positions = selectedIons.peakPositions{i};
    theoPattern = selectedIons.theoPattern{i};
    ionName = selectedIons.ionName{i};

    for p = 1:numel(positions)
        if positions(p) == 0, continue; end
        pk = round(positions(p), 4);

        % Only assign if this ion has higher score than any existing assignment
        if ~peakToIon.isKey(pk)
            peakToIon(pk) = ionName;
        end
    end
end

fprintf('\nApplying ranges:\n');
nOk = 0;
for i = 1:height(eerRanges)
    lo = eerRanges.mcbegin(i);
    hi = eerRanges.mcend(i);
    pk = round(eerRanges.peakMc(i), 4);

    if peakToIon.isKey(pk)
        ionName = peakToIon(pk);
    else
        % Try to find closest key
        allKeys = cell2mat(peakToIon.keys());
        [minDist, closestIdx] = min(abs(allKeys - eerRanges.peakMc(i)));
        if minDist < 0.3
            ionName = peakToIon(allKeys(closestIdx));
        else
            fprintf('  [%7.2f, %7.2f] Da — no ion assignment, skipping\n', lo, hi);
            continue;
        end
    end

    % Extract just the ion name for rangeAdd (strip trailing spaces)
    ionName = strtrim(ionName);

    try
        rangeAdd(spec, colorScheme, ionName, [lo hi]);
        nOk = nOk + 1;
        fprintf('  %-20s [%7.2f, %7.2f] Da\n', ionName, lo, hi);
    catch e
        fprintf('  %-20s [%7.2f, %7.2f] — %s\n', ionName, lo, hi, e.message);
    end
end

fprintf('\n%d of %d ranges applied.\n', nOk, height(eerRanges));

%% 8. Extract Range Table
rangeTable = rangesExtractFromMassSpec(spec);
fprintf('\nFinal range table (%d ranges):\n', height(rangeTable));
disp(rangeTable(:, {'rangeName', 'mcbegin', 'mcend', 'chargeState'}));

title('Automated mass spectrum analysis — fully ranged');

%% 9. Summary
totalRanged = 0;
for i = 1:height(rangeTable)
    totalRanged = totalRanged + sum(pos.mc >= rangeTable.mcbegin(i) & pos.mc <= rangeTable.mcend(i));
end
fprintf('\n========================================\n');
fprintf('  Elements: %s\n', strjoin(elements, ', '));
fprintf('  Peaks detected: %d (significant: %d)\n', height(peakTable), height(sigPeaks));
fprintf('  Ions matched: %d (above threshold: %d)\n', height(matchedIons), height(selectedIons));
fprintf('  Ranges applied: %d\n', nOk);
fprintf('  Ions ranged: %d / %d (%.1f%%)\n', totalRanged, height(pos), 100*totalRanged/height(pos));
fprintf('========================================\n');
