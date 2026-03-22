%[text] # Automated Mass Spectrum Analysis
%[text] This workflow performs automated peak detection, ion identification via isotopic pattern matching, and ranging on a calibrated atom probe mass spectrum.
%[text] **Steps:**
%[text] 1. Load data and create mass spectrum
%[text] 2. Peak detection with adjustable parameters
%[text] 3. Specify elements and generate candidate ions
%[text] 4. Match candidates to detected peaks via isotopic pattern matching
%[text] 5. Review and edit ion assignments
%[text] 6. Add ion stems to mass spectrum
%[text] 7. Create ranges using auto EER algorithm
%[text] 8. Apply ranges (auto-detect ions from stems) \

%%
%[text] ## 1. Setup and Load Data

setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme.mat');

%%
%[text] Load dataset — edit the path below or leave empty for a file dialog.

pos = posLoad();

%%
%[text] ## 2. Peak Detection
%[text] Adjust the parameters below to control sensitivity. Re-run this section until the detection captures all relevant peaks without too much noise.

binWidth            = 0.01;   % Da — histogram bin width
minProminenceFactor = 6;      % prominence must exceed factor x local noise
minDistanceDa       = 0.3;    % Da — minimum distance between peaks
baselineQuantile    = 0.1;    % quantile for baseline estimation
smoothSpan          = 0.1;    % Da — smoothing window for residual
minProminenceDisplay = 500;   % minimum prominence for a peak to be considered

[peakTable, peakInfo] = massSpecFindPeaks(pos.mc, ...
    'binWidth', binWidth, ...
    'minProminenceFactor', minProminenceFactor, ...
    'minDistanceDa', minDistanceDa, ...
    'baselineQuantile', baselineQuantile, ...
    'smoothSpan', smoothSpan);

sigPeaks = peakTable(peakTable.prominence > minProminenceDisplay, :);

fprintf('Detected %d peaks (%d above prominence threshold).\n', ...
    height(peakTable), height(sigPeaks))

%%
%[text] ## 2b. Diagnostic Plot
%[text] Inspect the spectrum with baseline, noise floor, and detected peaks. If too many or too few peaks are shown, adjust the parameters above and re-run section 2.

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

%%
%[text] ## 3. Specify Elements
%[text] List the chemical elements present in your specimen. The algorithm generates all plausible ionic and molecular species (via `ionsCreateComplex`), computes their isotopic fingerprints, and matches against detected peaks.
%[text] Physically impossible charge states (e.g. $H^{2+}$, $H_2^{2+}$) are automatically filtered.

elements = {'Ti', 'Al', 'N', 'O', 'H', 'C', 'Si'};

complexity     = [1 2 3];  % atom counts per ion (1=atomic, 2=diatomic, 3=triatomic)
chargeStates   = [1 2 3];  % charge states to consider
matchTolerance = 0.15;     % Da — isotope-to-peak matching tolerance

%%
%[text] ## 4. Isotopic Pattern Matching
%[text] For each candidate ion, the full isotopic pattern (all isotope positions and relative abundances) is compared against the observed peak heights using cosine similarity.
%[text] **Scoring:** $\text{score} = \text{cosine similarity} \times \text{fraction matched} \times \text{confidence} \times (1 / \text{complexity})$
%[text] Ions with many isotopes matching the correct ratios (e.g. Ti with 5 isotopes) score highest. Monoisotopic ions score lower. Molecular ions are penalised by complexity.

[matchedIons, ~] = ionMatchPattern(peakTable, elements, isotopeTable, ...
    'complexity', complexity, ...
    'chargeStates', chargeStates, ...
    'tolerance', matchTolerance, ...
    'minPeakProminence', minProminenceDisplay, ...
    'minScore', 0.05, ...
    'showPlot', true);

%%
%[text] ### Top Ion Identifications

nShow = min(40, height(matchedIons));
matchedIons(1:nShow, {'ionName', 'chargeState', 'score', 'complexity', 'nIsotopes'})

%%
%[text] ## 5. Select Ions to Keep
%[text] Choose a score threshold. Ions below the threshold are discarded. Review the list and manually remove any incorrect assignments in the next section.

scoreThreshold = 0.15;

selectedIons = matchedIons(matchedIons.score >= scoreThreshold, :);

fprintf('Selected %d ions with score >= %.2f\n', height(selectedIons), scoreThreshold)

%%
%[text] ### Manual Edits (optional)
%[text] Remove incorrect assignments by uncommenting and editing lines like:

% selectedIons(strcmp(selectedIons.ionName, 'Si C H ++'), :) = [];

%%
%[text] ## 6. Add Ion Stems to Mass Spectrum
%[text] Creates a fresh mass spectrum and adds stems for all selected ions (both atomic and molecular). `rangeAdd` will auto-detect ion identity from these stems when applying ranges.

close all;
spec = massSpecPlot(pos.mc, binWidth, 'normalised');
xlim([0 min(150, max(pos.mc))]);

stemmed = {};
for i = 1:height(selectedIons)
    ionName = selectedIons.ionName{i};
    cs = selectedIons.chargeState(i);
    tokens = regexp(ionName, '(.+?)\s*\++$', 'tokens');
    if isempty(tokens), continue; end
    formula = strtrim(tokens{1}{1});
    key = sprintf('%s_%d', formula, cs);
    if ismember(key, stemmed), continue; end
    try
        ionAdd(spec, formula, cs, isotopeTable, colorScheme, 0, 0.005, 'most abundant', 0.5);
        stemmed{end+1} = key; %#ok<SAGROW>
    catch
    end
end

fprintf('%d ion stems added to spectrum.\n', numel(stemmed))

%%
%[text] ## 7. Compute EER Ranges
%[text] Collects all peak $m/c$ positions from the selected ions' isotope patterns and computes optimal range boundaries using the Equal Error Rate criterion. The EER boundary is where missed signal equals included background.

allPeakMc = [];
for i = 1:height(selectedIons)
    positions = selectedIons.peakPositions{i};
    positions = positions(positions > 0);
    allPeakMc = [allPeakMc, positions]; %#ok<AGROW>
end
allPeakMc = unique(allPeakMc);

fprintf('Computing EER ranges for %d peak positions...\n', numel(allPeakMc))

[eerRanges, eerInfo] = rangeAutoEER(pos.mc, allPeakMc, 'binWidth', binWidth, 'showPlot', true);

fprintf('EER ranges computed: %d\n', height(eerRanges))

%%
%[text] ## 8. Apply Ranges
%[text] Each EER range is applied using `rangeAdd`, which auto-detects the ion from the stems on the spectrum.
%[text] - If **one** ion stem is in the range, it is assigned automatically.
%[text] - If **multiple** stems overlap, a selection dialog appears showing *"Select ion for range at XX.X Da"*. Pick the correct ion.
%[text] - If **no** stem is in range, the range is skipped. You can add it manually afterward.
%[text] - Narrow ranges (< 2 bins) are skipped automatically. \

minRangeWidth = 2 * binWidth;

nOk = 0;
nNoIon = 0;
nSkip = 0;
for i = 1:height(eerRanges)
    lo = eerRanges.mcbegin(i);
    hi = eerRanges.mcend(i);
    if (hi - lo) < minRangeWidth
        nSkip = nSkip + 1;
        continue;
    end
    try
        rangeAdd(spec, colorScheme, [], [lo hi]);
        nOk = nOk + 1;
    catch e
        if contains(e.message, 'no ion defined')
            nNoIon = nNoIon + 1;
        else
            nSkip = nSkip + 1;
        end
    end
end

fprintf('%d ranges applied, %d no ion in range, %d skipped.\n', nOk, nNoIon, nSkip)

%%
%[text] ## 9. Summary

ax = spec.Parent;
plots = ax.Children;
nRanges = 0;
totalRanged = 0;
for pl = 1:numel(plots)
    try
        if plots(pl).UserData.plotType == "range"
            nRanges = nRanges + 1;
            totalRanged = totalRanged + sum( ...
                pos.mc >= min(plots(pl).XData) & pos.mc <= max(plots(pl).XData));
        end
    catch
    end
end

fprintf('Elements: %s\n', strjoin(elements, ', '))
fprintf('Peaks detected: %d (significant: %d)\n', height(peakTable), height(sigPeaks))
fprintf('Ions matched: %d (above threshold: %d)\n', height(matchedIons), height(selectedIons))
fprintf('Ion stems added: %d\n', numel(stemmed))
fprintf('EER ranges computed: %d\n', height(eerRanges))
fprintf('Ranges applied: %d\n', nRanges)
fprintf('Ions ranged: %d / %d (%.1f%%)\n', totalRanged, height(pos), 100*totalRanged/height(pos))

%%
%[text] ## 10. Next Steps
%[text] After ranging, continue with the standard toolbox workflow:

%   rangeTable = rangesExtractFromMassSpec(spec);
%   pos = posAllocateRange(pos, rangeTable, 'decompose');
%   conc = posCalculateConcentrationSimple(pos, detEff, {'unranged'});
%   scatterPlotPosWidget(pos, colorScheme);

%%
%[text] ---
%[text] *Part of the [Atom Probe Toolbox](README.md). (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg*

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
