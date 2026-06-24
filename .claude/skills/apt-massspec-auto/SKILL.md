---
name: apt-massspec-auto
description: Automated mass spectrum analysis — peak detection, ion identification by isotopic pattern matching, and auto-ranging (EER). Use when the user wants ions and ranges found automatically instead of picking them by hand (contrast with the manual /apt-spectrum).
argument-hint: "[elements]"
---

Run fully automated peak detection, ion identification via isotopic pattern matching, and EER auto-ranging on a calibrated mass spectrum. Use the MATLAB MCP server. This is the automated counterpart to the manual `/apt-spectrum` workflow.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data with a calibrated `mc` column (from `/apt-load`)
- `isotopeTable` and `colorScheme` — from `/apt-setup` (or `load('isotopeTable_naturalAbundances.mat')` and `load('colorScheme.mat')`)

If `pos` is missing, tell the user to run `/apt-load` first. If the lookup tables are missing, run `/apt-setup`.

## Steps

### 1. Peak detection (interactive GUI)

```matlab
[peakTable, peakInfo, detParams] = tunePeakDetection(pos, 'mcRange', [0 min(200, max(pos.mc))]);
```

This opens an interactive figure. Tell the user to adjust the sliders and press **Accept** when satisfied:
- **Prominence factor** — how far a peak must rise above the noise
- **Min distance** — minimum spacing between peaks (Da)
- **Min prominence** — display threshold for significant peaks
- **Baseline quantile** — quantile for background estimation
- **Smooth span** — smoothing window for the residual signal

Do not try to automate the slider interaction. Defaults: `binWidth=0.01`, `minProminenceFactor=6`, `minDistanceDa=0.3`, `baselineQuantile=0.1`, `smoothSpan=0.1`, `minProminence=500`.

Then extract parameters and the significant peaks:
```matlab
binWidth             = detParams.binWidth;
minProminenceDisplay = detParams.minProminence;
sigPeaks = peakTable(peakTable.prominence > minProminenceDisplay, :);
```

### 2. Specify elements

```matlab
elements       = {'Ti', 'Al', 'N', 'O', 'H', 'C', 'Si'};  % from $ARGUMENTS or ask the user
complexity     = [1 2 3];   % atom counts per ion (1=atomic, 2=diatomic, 3=triatomic)
chargeStates   = [1 2 3];   % charge states to consider
matchTolerance = 0.15;      % Da — isotope-to-peak matching tolerance
```

If the user provided elements in `$ARGUMENTS`, use them. Otherwise ask which elements are present in the specimen — this list drives candidate ion generation.

### 3. Isotopic pattern matching

```matlab
[matchedIons, ~] = ionMatchPattern(peakTable, elements, isotopeTable, ...
    'complexity', complexity, ...
    'chargeStates', chargeStates, ...
    'tolerance', matchTolerance, ...
    'minPeakProminence', minProminenceDisplay, ...
    'minScore', 0.05, ...
    'showPlot', true);
nShow = min(40, height(matchedIons));
matchedIons(1:nShow, {'ionName', 'chargeState', 'score', 'complexity', 'nIsotopes'})
```

Each candidate ion's full isotopic pattern is compared against observed peak heights by cosine similarity; molecular ions are penalised by complexity.

### 4. Select ions to keep

```matlab
scoreThreshold = 0.15;  % adjust this value
selectedIons = matchedIons(matchedIons.score >= scoreThreshold, :);
fprintf('Selected %d ions with score >= %.2f\n', height(selectedIons), scoreThreshold)
selectedIons(:, {'ionName', 'chargeState', 'score', 'complexity', 'nIsotopes'})
```

Show the list and let the user raise/lower the threshold or manually drop wrong assignments, e.g.:
```matlab
% selectedIons(strcmp(selectedIons.ionName, 'Si C H ++'), :) = [];
```

### 5. Add ion stems to the spectrum

```matlab
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
```

`rangeAdd` later auto-detects the ion from these stems.

### 6. Compute EER ranges

```matlab
allPeakMc = [];
for i = 1:height(selectedIons)
    positions = selectedIons.peakPositions{i};
    positions = positions(positions > 0);
    allPeakMc = [allPeakMc, positions]; %#ok<AGROW>
end
% Merge with ALL significant detected peaks for full coverage
allPeakMc = unique([allPeakMc, sigPeaks.mc']);
fprintf('Computing EER ranges for %d peak positions...\n', numel(allPeakMc))

[eerRanges, eerInfo] = rangeAutoEER(pos.mc, allPeakMc, 'binWidth', binWidth, 'showPlot', true);
fprintf('EER ranges computed: %d\n', height(eerRanges))
```

The Equal Error Rate boundary is where missed signal equals included background.

### 7. Apply ranges

First pre-register any molecular-ion formulas in `colorScheme` so `rangeAdd` can colour them:
```matlab
for i = 1:height(selectedIons)
    tokens = regexp(selectedIons.ionName{i}, '(.+?)\s*\++$', 'tokens');
    if isempty(tokens), continue; end
    formula = strtrim(tokens{1}{1});
    catName = categorical(string(formula));
    if ~any(colorScheme.ion == catName)
        firstElem = regexp(formula, '^[A-Z][a-z]?', 'match', 'once');
        if ~isempty(firstElem) && any(colorScheme.ion == categorical(string(firstElem)))
            baseColor = colorScheme.color(colorScheme.ion == categorical(string(firstElem)), :);
            newColor = min(1, baseColor(1,:) + 0.15 * randn(1,3));
        else
            newColor = rand(1,3) * 0.5 + 0.25;
        end
        colorScheme = [colorScheme; table(catName, newColor, 'VariableNames', {'ion','color'})]; %#ok<AGROW>
    end
end

minRangeWidth = 2 * binWidth;
nOk = 0; nNoIon = 0; nSkip = 0;
for i = 1:height(eerRanges)
    lo = eerRanges.mcbegin(i);
    hi = eerRanges.mcend(i);
    if (hi - lo) < minRangeWidth, nSkip = nSkip + 1; continue; end
    try
        rangeAdd(spec, colorScheme, [], [lo hi]);
        nOk = nOk + 1;
    catch e
        if contains(e.message, 'no ion defined'), nNoIon = nNoIon + 1; else, nSkip = nSkip + 1; end
    end
end
fprintf('%d ranges applied, %d no ion in range, %d skipped.\n', nOk, nNoIon, nSkip)
```

`rangeAdd(spec, colorScheme, [], [lo hi])` auto-assigns the ion when exactly one stem falls in the range. If **multiple** stems overlap, a *"Select ion for range at XX.X Da"* dialog appears — tell the user to pick the correct ion. Ranges with no stem are skipped.

## Report

Summarise: elements supplied, peaks detected (and significant), ions matched / above threshold, ion stems added, EER ranges computed, ranges applied, and percent of atoms ranged.

The spectrum `spec` now carries ions and ranges, ready for downstream skills:
- `/apt-composition` — extract ranges and compute concentration
- Manual continuation:
  ```matlab
  rangeTable = rangesExtractFromMassSpec(spec);
  pos = posAllocateRange(pos, rangeTable, 'decomposed');
  conc = posCalculateConcentrationSimple(pos, detEff, {'unranged'}, '', 'mode', 'atomic');
  scatterPlotPosWidget(pos, colorScheme);
  ```
