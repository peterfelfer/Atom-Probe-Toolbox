---
name: apt-profile-1d
description: Compute and plot a 1D concentration profile (composition vs. distance) from ranged/allocated atom probe data by voxelising along one axis. Use when the user wants a concentration profile, line profile, or composition-along-z plot.
argument-hint: "[bin-nm]"
---

Build a 1D concentration profile by slicing the reconstructed tip into bins along one distance axis and computing the composition of each slice. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — ranged/allocated atom probe data with an `ion`/`atom` column (from `/apt-composition`)
- `colorScheme` — ion color table (from `/apt-setup`)

If `pos` has no `ion` column, tell the user to run `/apt-composition` first (ranging + `posAllocateRange`). If `colorScheme` is missing, run `/apt-setup`.

## Steps

### 1. Choose the distance axis and create bin vectors

```matlab
dist = pos.z;                 % distance variable to bin along (z is the default profile direction)
bin  = 1;                     % bin width in nm
mode = 'distance';
[binCenters, binEdges] = binVectorsFromDistance(dist, bin, mode);
```

`bin` is the slice thickness in nm. Default 1 nm; use 0.25 nm near interfaces, 0.5 nm for a typical profile. Parse `$ARGUMENTS` for a bin width if the user supplied one. To profile along a different direction use `pos.x` or `pos.y`, or supply an ROI/interface distance from `/apt-roi`. There is also a `'count'` mode (equal atoms per bin, e.g. `bin = 390000`) used mainly for proxigrams — only use it if the user explicitly asks.

### 2. Voxelise the data

```matlab
vox = posNdBin(pos, dist, binEdges);   % cell array: one pos subset per slice
```

The last slice is usually thinner / has fewer atoms than the rest — this is expected.

### 3. Define the concentration kernel (set detection efficiency)

```matlab
detEff          = 37;             % detector efficiency, percent or decimal both accepted
excludeListConc = {'unranged'};   % species ignored in the composition, optional
volumeName      = '';             % volume label, optional
concentrationKernel = @(pos) posCalculateConcentrationSimple(pos, detEff, excludeListConc, volumeName);
```

The detection efficiency is critical for correct composition and must match the instrument: 37 (LEAP 4000X HR), 52 (LEAP 5000 XR), 80 (LEAP 5000 XS), 50 (EIKOS). If the user has not confirmed their instrument, ask before proceeding. Read it from config if available: `cfg = APTConfig.getInstance(); detEff = cfg.reconstruction.detectionEfficiency*100;`.

### 4. Apply the kernel to every slice

```matlab
conc1D = binApplyConcentrationKernel(vox, binCenters, concentrationKernel, {'nm'});
```

This returns one table holding the concentration, counts and variance per slice, tagged with the bin-center distance in nm. Composition is reported with counting-based variance, so each point carries a statistical uncertainty (see `concentrationUncertainty` below).

### 5. Plot the profile

```matlab
excludeListPlot = {'unranged'};   % species hidden from the plot, optional
[p, ax, f] = concentrationProfilePlot(conc1D(conc1D.format=='concentration', :), excludeListPlot, colorScheme);
```

This plots concentration vs. distance for all ranged species (up to 4 volumes can be overlaid). `ax` is left in the workspace for further tweaking.

### 6. (Optional) Extract values for one species

```matlab
conc = cellfun(concentrationKernel, vox, 'UniformOutput', false);
concWantedSpecies = cellfun(@(x) x.Fe(x.format == 'concentration'), conc);   % replace Fe with element/complex ion
```

`concWantedSpecies` is a column vector of the fraction of that species per slice, aligned with `binCenters{1}`.

### 7. (Optional) Counting uncertainty

```matlab
unc = concentrationUncertainty(conc1D.counts);
```

Reports the statistical (counting) uncertainty on the profile values.

## Report

Tell the user the bin width used, the distance axis, and the detection efficiency applied, and that a 1D profile figure was created. Mention that each point has counting uncertainty (smaller bins = noisier profiles). Suggest:
- `/apt-roi` to define a profile direction or sub-region (e.g. across an interface) before re-running
- `/apt-profile-2d` for a 2D concentration map
- Export: `writetable(conc1D, 'profile_1d.csv')` or save the figure with `exportgraphics(ax, 'profile_1d.png')`
