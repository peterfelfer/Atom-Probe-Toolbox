---
name: apt-spectrum
description: Create a mass spectrum from loaded atom probe data and optionally identify ions. Use when the user wants to see or analyze a mass spectrum.
argument-hint: "[element1 element2 ...]"
---

Create a mass spectrum from the `pos` variable in the MATLAB workspace. Use the MATLAB MCP server.

## Steps

### 1. Plot the mass spectrum

```matlab
spec = massSpecPlot(pos, 0.01, 'normalised');
```

Use bin width 0.01 Da and normalised mode by default. If the user requests different settings, adjust accordingly.

### 2. Identify ions (if elements provided)

If the user specified elements as `$ARGUMENTS`, add ion markers for each element. For each element, try common charge states (typically 1+ and 2+ for most elements, 3+ for transition metals):

```matlab
ionAdd(spec, 'Fe', 2, isotopeTable, colorScheme);
ionAdd(spec, 'Fe', 1, isotopeTable, colorScheme);
```

If no elements are specified, just show the spectrum and ask the user which elements they expect in their sample. Suggest looking at the major peaks for identification.

### 3. Report

Tell the user:
- The spectrum is displayed in the MATLAB figure
- Which ions were identified (if any)
- Suggest next steps: identify more ions with `/apt-spectrum element1 element2`, or define ranges with `rangeAddAll`

## Notes

- `isotopeTable` and `colorScheme` must be in the workspace (from `/apt-setup`)
- If they're missing, load them first:
  ```matlab
  load('isotopeTable_naturalAbundances.mat');
  load('colorScheme_default.mat');
  ```
