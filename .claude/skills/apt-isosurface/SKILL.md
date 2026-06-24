---
name: apt-isosurface
description: Voxelise atom probe data, compute a concentration field, and render a 3D isosurface at a chosen at.% threshold (e.g. for interfaces, precipitates, or grain boundaries). Use when the user wants an isoconcentration surface or isosurface from their APT data.
argument-hint: "[isovalue]"
---

Voxelise the reconstruction, compute the concentration of one or more species per voxel, and draw an isosurface at a chosen at.% level. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data with ions allocated (has `ion` column, from `/apt-composition`)
- `colorScheme` — color scheme (from `/apt-setup`)

If ions have not been allocated yet, tell the user to run `/apt-composition` first.

## Steps

### 1. Create grid vectors

Build a 3D binning grid from the atom coordinates. `bin` is the voxel size in nm and may be isotropic `[2 2 2]` or anisotropic `[1 2 2]`. Mode is always `'distance'` for 3D voxelisation:

```matlab
dist = [pos.x pos.y pos.z];
bin  = [2 2 2];            % voxel size in nm (1 = fine detail, 2 = general, 3 = fast overview)
mode = 'distance';
[binCenters, binEdges] = binVectorsFromDistance(dist, bin, mode);
```

### 2. Voxelise and compute concentration

Voxelise all atoms, then voxelise only the species of interest, and take the ratio to get the per-voxel concentration:

```matlab
gridVec = binCenters;
vox     = posToVoxel(pos, gridVec);            % counts of all species per voxel
species = {'Al', 'Ti'};                        % species defining the isosurface; ask the user
voxIon  = posToVoxel(pos, gridVec, species);   % counts of selected species per voxel
conc    = voxIon ./ vox;                        % concentration field (fraction, 0..1)
```

Ask the user which species the isosurface should be built on. `conc` is a fraction; an isovalue of `8` below means 8 at.%.

### 3. Create and render the isosurface

```matlab
isovalue = 8;   % at.% threshold (typical 1-40 depending on the system)
fv = isosurface(gridVec{2}, gridVec{1}, gridVec{3}, conc, isovalue);
RGB = [1 1 0];  % face color
p = patch(fv, 'FaceColor', RGB); axis equal; rotate3d on;
axisSpatialAptify;   % switch to the standard APT viewing orientation
```

Use `$ARGUMENTS` for the isovalue if the user gave one; otherwise start at `8` and suggest trying several values to best capture the feature.

**Critical axis-ordering note:** MATLAB's `isosurface` treats the first array dimension as Y and the second as X, so the grid vectors are passed as `gridVec{2}, gridVec{1}, gridVec{3}` (X and Y swapped). Do not "fix" this to `{1},{2},{3}` — if the surface looks mirrored or rotated, this swap is what keeps it aligned with the data.

If `conc` is single precision and `isosurface` complains, cast to double: `conc = double(conc);`.

### 4. (Optional) Interactive isovalue tuning

Live scripts cannot drive the app, so tell the user to paste this into the MATLAB **command window** themselves (manual GUI step — do not automate). It opens a slider to sweep the isovalue live, recolor patches, and split/export patches to the workspace:

```matlab
species    = {'Cr'};                   % numerator species
speciesAll = {'Al','Ni','Co','Cr'};    % denominator species (must include numerator); {} = use raw counts
bin        = 2;                        % voxel size in nm
isosurfaceWidget(pos, species, speciesAll, bin, colorScheme, ax);
```

### 5. (Optional) Overlay ions with the isosurface

Plot ions first, then add the isosurface patch into the same axes:
```matlab
species = {'C'};
sample  = 1;     % fraction (<1) or absolute number (>1) of atoms to plot
size    = 1;     % marker size
pPlot = scatterPlotPosData(pos, species, sample, colorScheme, size, axes(figure()));
ax = pPlot.Parent;
patch(ax, fv, 'DisplayName', 'Isosurface');
```
Or add ions into an existing isosurface plot: `ax = gca; scatterPlotPosData(pos, {'Mo','B'}, sample, colorScheme, size, ax);`

### 6. (Optional) Export / import the surface

Export to OBJ (e.g. for Blender):
```matlab
patchToObj(fv, 'grain boundary', 'fileName.obj');
```
Re-import an edited OBJ (call with no argument to pick a file via dialog):
```matlab
[obj, gr] = objToPatch('isosurface.obj');
```

## Report

Tell the user:
- The voxel size, species, and isovalue used.
- That `fv` (faces/vertices) and `conc` (concentration field) are in the workspace.
- Suggest trying other isovalues, and the `isosurfaceWidget` for interactive tuning.

Suggested follow-ups: `/apt-proxigram` (proxigram across this isosurface) or exporting via `patchToObj` for external rendering.
