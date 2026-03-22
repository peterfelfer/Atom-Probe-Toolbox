%% Automated Mass Spectrum Analysis
% This workflow performs automated peak detection, ion identification, and
% ranging on a calibrated atom probe mass spectrum.
%
% Steps:
%   1. Load data and create mass spectrum
%   2. Peak detection with adjustable parameters
%   3. Specify elements and generate candidate ions
%   4. Match candidates to detected peaks
%   5. Create ranges using auto EER algorithm
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
minProminenceFactor = 6;      % prominence must exceed factor × local noise
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
% Inspect the detected peaks. Adjust minProminenceDisplay above if needed.
disp(sigPeaks(:, {'mc', 'height', 'prominence'}));

%% 3. Specify Elements and Ion Generation Parameters
% List the chemical elements present in your specimen.
% The algorithm will generate all plausible ionic and molecular species
% from these elements and match them against the detected peaks.

% --- Edit these for your specimen ---
elements = {'Pd', 'H', 'C', 'N', 'O'};

% --- Ion generation parameters ---
complexity     = [1 2 3];  % atom counts per ion (1=atomic, 2=diatomic, 3=triatomic)
chargeStates   = [1 2];    % charge states to consider
matchTolerance = 0.3;      % Da — how close theoretical mc must be to a detected peak

%% 4. Match Candidate Ions to Detected Peaks
% Uses ionsCreateComplex to generate all possible ions from the element
% list, then matches their isotopic mc values against detected peaks.

[matchedIons, ionList] = ionMatchPeaks(peakTable, elements, isotopeTable, ...
    'complexity', complexity, ...
    'chargeStates', chargeStates, ...
    'tolerance', matchTolerance, ...
    'minPeakProminence', minProminenceDisplay);

fprintf('\nGenerated %d isotopic candidates from %d elements.\n', height(ionList), numel(elements));

% Show best matches (one per peak)
bestMatches = matchedIons(matchedIons.isBestMatch, :);
fprintf('\nBest ion assignment per peak (%d):\n', height(bestMatches));
disp(bestMatches(:, {'ionName', 'theoreticalMc', 'peakMc', 'chargeState', 'peakHeight', 'abundance'}));

%% 4b. Review and Edit Ion Assignments
% The automatic matching may not be perfect. Review the table above and
% remove incorrect assignments or add missing ones.
%
% To remove an ion: delete the row from bestMatches
%   bestMatches(strcmp(bestMatches.ionName, 'C O+'), :) = [];
%
% To add an ion manually:
%   newRow = table({"N2+"}, 28.006, 1, 28.01, 30000, 99.6, {"N2"}, true, ...
%       'VariableNames', bestMatches.Properties.VariableNames);
%   bestMatches = [bestMatches; newRow];

%% 5. Create Mass Spectrum with Ion Stems
% Plot the mass spectrum and overlay ion stems for all matched ions.

spec = massSpecPlot(pos.mc, binWidth, 'normalised');
xlim([0 min(150, max(pos.mc))]);

% Add stems for each unique element at each matched charge state
stemmedChargeStates = struct();
for i = 1:height(bestMatches)
    formula = bestMatches.formula{i};
    q = bestMatches.chargeState(i);

    % Only add atomic ion stems (molecular ions can't go through ionAdd)
    % Check if formula is a single element
    tokens = strsplit(strtrim(formula));
    if numel(tokens) == 1 && ~any(isstrprop(formula, 'digit'))
        key = sprintf('%s_%d', formula, q);
        if ~isfield(stemmedChargeStates, matlab.lang.makeValidName(key))
            try
                ionAdd(spec, formula, q, isotopeTable, colorScheme, 0, 0.005, 'most abundant', 0.5);
                stemmedChargeStates.(matlab.lang.makeValidName(key)) = true;
            catch
            end
        end
    end
end

title('Mass spectrum with matched ion stems');

%% 6. Auto-Range with EER Algorithm
% Use the matched peak positions as input to rangeAutoEER.

peakMcForRanging = bestMatches.peakMc;
[eerRanges, eerInfo] = rangeAutoEER(pos.mc, peakMcForRanging, ...
    'binWidth', binWidth, 'showPlot', true);

fprintf('\nEER ranges computed for %d peaks.\n', height(eerRanges));

%% 7. Apply Ranges to Mass Spectrum
% Apply each EER range using rangeAdd with the matched ion name.

% Ensure molecular ions are in colorScheme
molecularIons = {};
for i = 1:height(bestMatches)
    formula = bestMatches.formula{i};
    tokens = strsplit(strtrim(formula));
    if numel(tokens) > 1
        catName = categorical(string(formula));
        if ~any(colorScheme.ion == catName)
            molecularIons{end+1} = formula;
            newColor = rand(1, 3) * 0.6 + 0.2;
            colorScheme = [colorScheme; table(catName, newColor, 'VariableNames', {'ion', 'color'})];
        end
    end
end
if ~isempty(molecularIons)
    fprintf('Added %d molecular ions to colorScheme: %s\n', ...
        numel(molecularIons), strjoin(molecularIons, ', '));
end

fprintf('\nApplying ranges:\n');
nOk = 0;
for i = 1:height(bestMatches)
    ionName = bestMatches.ionName{i};
    targetMc = bestMatches.peakMc(i);

    [dist, matchIdx] = min(abs(eerRanges.peakMc - targetMc));
    if dist > 1.0
        fprintf('  %-20s mc=%6.2f — no EER range\n', ionName, targetMc);
        continue;
    end

    lo = eerRanges.mcbegin(matchIdx);
    hi = eerRanges.mcend(matchIdx);

    try
        rangeAdd(spec, colorScheme, ionName, [lo hi]);
        nOk = nOk + 1;
        fprintf('  %-20s [%7.2f, %7.2f] Da\n', ionName, lo, hi);
    catch e
        fprintf('  %-20s [%7.2f, %7.2f] — %s\n', ionName, lo, hi, e.message);
    end
end

fprintf('\n%d of %d ranges applied.\n', nOk, height(bestMatches));

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
fprintf('  Peaks detected: %d\n', height(peakTable));
fprintf('  Ions matched: %d\n', height(bestMatches));
fprintf('  Ranges applied: %d\n', nOk);
fprintf('  Ions ranged: %d / %d (%.1f%%)\n', totalRanged, height(pos), 100*totalRanged/height(pos));
fprintf('========================================\n');
