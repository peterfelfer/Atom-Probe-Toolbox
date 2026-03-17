---
name: apt-analyze
description: Run a full atom probe analysis pipeline on a data file. Loads data, creates mass spectrum, and guides through ion identification, ranging, and composition calculation. Use for end-to-end APT analysis.
argument-hint: "<file.pos|file.epos>"
---

Run a complete atom probe analysis pipeline. Use the MATLAB MCP server.

The file to analyze is: `$ARGUMENTS`

If no file is provided, ask the user for the file path.

## Full Pipeline

### Step 1: Setup (if not already done)

```matlab
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme_default.mat');
cfg = APTConfig.getInstance();
```

### Step 2: Ask about the instrument

Before proceeding, ask the user which instrument acquired this data:
- LEAP 4000X HR (detEff = 0.37)
- LEAP 5000 XR (detEff = 0.52)
- LEAP 5000 XS (detEff = 0.80)
- EIKOS (detEff = 0.50)

Apply the preset:
```matlab
cfg.applyInstrumentPreset('LEAP5000XR');  % adjust based on answer
```

### Step 3: Load data

```matlab
pos = posLoad('$ARGUMENTS');
fprintf('Loaded %d atoms\n', height(pos));
```

### Step 4: Mass spectrum

```matlab
spec = massSpecPlot(pos, 0.01, 'normalised');
```

### Step 5: Ask about expected elements

Ask the user what elements they expect in their sample (e.g., "Fe, C, Mn" for steel). Then add ion markers:

```matlab
ionAdd(spec, 'Fe', 2, isotopeTable, colorScheme);
ionAdd(spec, 'C', 1, isotopeTable, colorScheme);
% ... for each element at likely charge states
```

### Step 6: Define ranges

Tell the user to check the spectrum in MATLAB and confirm the ion identification looks correct. Then:

```matlab
rangeAddAll(spec, colorScheme, 0.04, true);
```

### Step 7: Extract ranges and allocate

```matlab
rangeTable = rangesExtractFromMassSpec(spec);
pos = posAllocateRange(pos, rangeTable, 'decomposed');
```

### Step 8: Composition

```matlab
detEff = cfg.reconstruction.detectionEfficiency;
conc = posCalculateConcentrationSimple(pos, detEff, {'unranged'}, '', 'mode', 'atomic');
unc = concentrationUncertainty(conc.counts);
disp(conc);
```

### Step 9: Visualization

```matlab
scatterPlotPosWidget(pos, colorScheme);
```

### Step 10: Report

Present a summary:
- Dataset: file name, atom count
- Instrument: name and detection efficiency used
- Composition table with uncertainties
- Suggest next steps: ROI analysis, proxigram, cluster analysis
