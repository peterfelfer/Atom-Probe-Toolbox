---
name: apt-load-rrng
description: Load an IVAS-style .rrng (or .rng) range file onto an existing mass spectrum, recreating all its ion ranges (and their colors) on the plot. Use when the user has a .rrng/.rng file and wants to apply those ranges instead of defining them by hand.
argument-hint: "<file.rrng>"
---

Load a `.rrng` / `.rng` range file and add every range it defines to the current mass spectrum `spec`. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `spec` — an open mass spectrum figure (from `/apt-spectrum`, e.g. `spec = massSpecPlot(pos, 0.01)`). The ranges are drawn onto this figure.
- `colorScheme` and `isotopeTable` — loaded by `/apt-setup` (or `load('colorScheme.mat')` and `load('isotopeTable_naturalAbundances.mat')`).

If `spec` does not exist, tell the user to run `/apt-spectrum` first.

## Steps

### 1. Locate the range file

`$ARGUMENTS` is the path to the `.rrng`/`.rng` file. If a path is given, use it directly; otherwise let the user browse:

```matlab
[file, path] = uigetfile('*.rrng');
```

The file-picker dialog is a manual GUI action — tell the user to select the file in the dialog. (If `$ARGUMENTS` is supplied, instead set `path = ''; file = '<file.rrng>';`.)

Make sure the lookup tables are available:
```matlab
load('colorScheme.mat')
load('isotopeTable_naturalAbundances.mat')
```

### 2. Parse elements and build the candidate ion list

```matlab
elements = elementsExtractFromText(fileread([path file]));
ionList  = ionsCreateComplex(elements, [1 2], isotopeTable, [1 2 3]);
```

`ionsCreateComplex(elements, complexity, isotopeTable, chargeStates)` here builds mono- and di-atomic complexes (`complexity = [1 2]`) for charge states 1-3 (`[1 2 3]`).

### 3. Split the file into range lines

```matlab
rrng = string(splitlines(fileread([path file])));
rrng(rrng == "") = [];
useRrngColor = true;     % use the colors stored in the .rrng file (set false to auto-assign)
```

### 4. Add each range to the spectrum

```matlab
for i = 1:length(rrng)
    [mcBegin, mcEnd, ionName, ionVolume, ionColor] = rangeInfoFromRRNG(rrng(i));
    if not(isempty(mcBegin))
        [h, txt, colorScheme] = rangeAddFromRangeInfo(spec, colorScheme, isotopeTable, ...
            ionList, mcBegin, mcEnd, ionName, ionVolume, ionColor, useRrngColor);
    end
end
```

`rangeInfoFromRRNG` returns empty `mcBegin` for non-range lines (headers, comments), which are skipped. `rangeAddFromRangeInfo` draws each range and may extend `colorScheme` with any new ion colors, so reassign it as shown.

## Report

Tell the user how many ranges were added to `spec` and confirm `colorScheme` may have gained new ions from the file. Note that `useRrngColor = true` preserves the original IVAS colors.

Suggest follow-ups:
- Extract the ranges into a table for allocation: `rangeTable = rangesExtractFromMassSpec(spec);`
- Compute composition: `/apt-composition`
- Tweak or add more ranges manually with `rangeAdd(spec, colorScheme)`.

## Toolbox functions used

`elementsExtractFromText`, `ionsCreateComplex`, `rangeInfoFromRRNG`, `rangeAddFromRangeInfo`; `uigetfile`, `fileread`, `splitlines` (MATLAB built-ins). Follow-ups: `rangesExtractFromMassSpec`, `rangeAdd`.
