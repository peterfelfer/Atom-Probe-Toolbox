---
name: apt-peak-background
description: Background-corrected peak counting — fit a linear background under a single mass-spectrum peak from interactively brushed bins and get its corrected counts/percentage. Use when the user wants accurate peak counts with background subtraction, signal-to-background, or to build a table of corrected peaks.
---

Compute background-corrected counts for individual mass-spectrum peaks by fitting a linear background from manually selected background bins. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data in RAW (undecomposed) format (from `/apt-load`)
- `spec` — an open, ranged mass spectrum (from `/apt-spectrum` or `/apt-massspec-auto`)
- `rangeTable` — extracted from the spectrum. If missing, create it:
  ```matlab
  rangeTable = rangesExtractFromMassSpec(spec);
  ```

If `pos` is decomposed, tell the user this workflow needs the raw pos file. If no ranged spectrum exists, run `/apt-spectrum` or `/apt-massspec-auto` first.

## Steps

### 1. Background-correct one peak (interactive)

The mass spectrum figure must be the current figure and zoomed to the peak of interest. Because the function requires interactive point selection, the user must run this line in the MATLAB **command window** (it will not work fired blindly from a script):

```matlab
[peakData] = peakBGCorrectedCount(pos, rangeTable, ionName);
```

`ionName` must be the specific-isotope form, e.g. `'14N 12C2 +'`. Tell the user to select **4 points** in the spectrum — two before and two after the peak — marking the bins used for the linear background fit. The function plots the peak with its fitted background and returns the `peakData` struct (fields include `counts`, `pct`, `loc`).

Optional arguments (confirmed signature `peakBGCorrectedCount(pos, rangeTable, ionName, boundary, options)`):
- `boundary` — Da before/after the peak for the background range, scalar or `[lo hi]` (default `0.5`)
- a `figure` flag (`true`/`false`) to suppress the plot
- options `'r'` (add a fictional range and report signal-to-background) or `'a'` (auto-create that range by balancing missed atoms vs. background atoms)

### 2. Start a results table (first peak)

```matlab
rangeTable = rangesExtractFromMassSpec(spec);
peakName        = rangeTable.rangeName(rangeTable.mcbegin < peakData.loc & rangeTable.mcend > peakData.loc);
peakChargeState = rangeTable.chargeState(rangeTable.mcbegin < peakData.loc & rangeTable.mcend > peakData.loc);
peakConc = table(peakName, peakChargeState, peakData.loc, peakData.counts, peakData.pct, ...
    'VariableNames', {'peakName', 'peakChargeState', 'peakLoc', 'peakCounts', 'peakPct'});
clearvars peakName peakChargeState
```

### 3. Append each additional peak

Re-run step 1 on the next peak, then append:

```matlab
rangeTable = rangesExtractFromMassSpec(spec);
peakName        = rangeTable.rangeName(rangeTable.mcbegin < peakData.loc & rangeTable.mcend > peakData.loc);
peakChargeState = rangeTable.chargeState(rangeTable.mcbegin < peakData.loc & rangeTable.mcend > peakData.loc);
peakNew = table(peakName, peakChargeState, peakData.loc, peakData.counts, peakData.pct, ...
    'VariableNames', {'peakName', 'peakChargeState', 'peakLoc', 'peakCounts', 'peakPct'});
if ~exist('peakConc','var'); peakConc = table(); end
peakConc = [peakConc; peakNew];
clearvars peakName peakChargeState
```

Repeat steps 1 and 3 for every peak the user wants. To compute a ratio between two elements, divide their `peakCounts` rows in `peakConc`.

## Report

Show `peakConc` with each peak's name, charge state, location (Da), background-corrected counts, and percentage. Note that counts here are background-subtracted, so element ratios are more reliable than raw range counts.

Follow-ups:
- Export: `writetable(peakConc, 'peakBackgroundCorrected.csv')`
- For bulk background-removed composition over the whole spectrum, see `posCalculateConcentrationBackgroundRemoved` (note: in `untested/`).
- Related helpers: `tofSpecBackgroundDetermination`, `backgroundCorrectedPeakCount`.
