---
name: apt-profile-2d
description: Compute and display a 2D concentration map (heat map of one species over a plane) from ranged/allocated atom probe data by voxelising along two axes. Use when the user wants a 2D concentration map, composition heat map, or chemical map.
argument-hint: "[species]"
---

Build a 2D concentration map by binning the tip into pixels over two distance axes and computing the composition of each pixel, then displaying a heat map for one species. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — ranged/allocated atom probe data with an `ion`/`atom` column (from `/apt-composition`)

If `pos` has no `ion` column, tell the user to run `/apt-composition` first (ranging + `posAllocateRange`).

## Steps

### 1. Choose the two distance axes and create bin vectors

```matlab
dist = [pos.y, pos.z];        % the two coordinates to bin (one column each)
bin  = [1 1];                 % bin width [dim1 dim2] in nm; non-isometric bins allowed
mode = 'distance';            % 'distance' is the only supported mode in 2D
[binCenters, binEdges] = binVectorsFromDistance(dist, bin, mode);
```

`bin` sets the pixel resolution and strongly affects compute time — 1 to 2 nm is reasonable. Pick the axis pair the user wants: `[pos.x pos.y]`, `[pos.x pos.z]`, or `[pos.y pos.z]`.

### 2. Voxelise the data

```matlab
startIdx = 1;                 % first ion hit index to include
endIdx   = height(pos);       % last ion hit index (must be > startIdx)
pix = posNdBin(pos(startIdx:endIdx,:), dist(startIdx:endIdx,:), binEdges);  % cell array of pos subsets (pixels)
```

By default include all atoms (`startIdx = 1`, `endIdx = height(pos)`). Use a narrower index range only if the user wants to restrict the hit range.

### 3. Define the concentration kernel (set detection efficiency)

```matlab
detEff          = 37;             % detector efficiency, percent or decimal both accepted
excludeListConc = {'unranged'};   % species ignored in the composition, optional
volumeName      = '';             % volume label, optional
concentrationKernel = @(pos) posCalculateConcentrationSimple(pos, detEff, excludeListConc, volumeName);
```

The detection efficiency is critical for correct composition and must match the instrument: 37 (LEAP 4000X HR), 52 (LEAP 5000 XR), 80 (LEAP 5000 XS), 50 (EIKOS). If the user has not confirmed their instrument, ask before proceeding. Read it from config if available: `cfg = APTConfig.getInstance(); detEff = cfg.reconstruction.detectionEfficiency*100;`.

### 4. Compute concentration for every pixel

```matlab
conc = cellfun(concentrationKernel, pix, 'UniformOutput', false);
```

This applies the kernel to each pixel, producing a concentration (plus count and variance) table per pixel. Composition therefore carries counting statistics per pixel (see `concentrationUncertainty`).

### 5. Build and display the heat map for one species

```matlab
wantedSpecies = cellfun(@(pos) pos.Al(pos.format == 'concentration'), conc);  % replace Al with element/complex ion
figure;
concentrationMap = imagesc(binCenters{1}, binCenters{2}, wantedSpecies);
axis equal;
```

Set the species to the element or complex ion the user asked for (parse `$ARGUMENTS`); default to a major solute if unspecified, and tell the user which one you picked. Add `colorbar; xlabel('nm'); ylabel('nm');` if helpful.

### 6. (Optional) Counting uncertainty

```matlab
counts = cellfun(@(x) sum(x.counts), conc);
unc = concentrationUncertainty(counts);
```

Sparse pixels (few atoms) have large statistical uncertainty — increase `bin` to average more atoms per pixel if the map is noisy.

## Report

Tell the user the axis pair and bin sizes used, the detection efficiency applied, and which species is mapped, and that a 2D heat map figure was created. Note that pixel composition carries counting uncertainty (larger bins = smoother, less noisy maps). Suggest:
- `/apt-roi` to restrict the map to a region of interest before re-running
- `/apt-profile-1d` for a 1D concentration profile along a single axis
- Re-run step 5 with a different `wantedSpecies` to map another element
- Export the figure with `exportgraphics(gca, 'concentration_map.png')`
