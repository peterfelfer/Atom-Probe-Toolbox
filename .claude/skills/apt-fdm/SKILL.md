---
name: apt-fdm
description: Create a field desorption map (detector hit-density / concentration image) and optionally a movie of detector activity over the course of the experiment. Use when the user wants an FDM, desorption map, detector image, or to watch features evolve during the measurement.
argument-hint: "[species]"
---

Create a field desorption map (FDM) from the detector hit coordinates, and optionally a movie of how the detector pattern evolves during the experiment. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — atom probe data **with detector coordinates** `detx`, `dety` (in nm). These exist only in extended formats (`.epos` / `.apt` / pyccapt HDF5), not plain `.pos`. If only a `.pos` is loaded, tell the user to reload via `/apt-load` with an `.epos`/`.apt`/HDF5 file.
- For a **species-specific** FDM only: `pos` must also be allocated (have an `ion`/`atom` column) — run `/apt-composition` first so a `spec` with ranges exists, then allocate. A plain hit-density FDM needs no ranging.

## Steps

### 1. Build the bin vectors (detector grid)

```matlab
dist = [pos.detx, pos.dety];
bin  = [0.07 0.07];          % bin width in nm for x and y -> FDM resolution
mode = 'distance';           % only valid 2D/3D mode
[binCenters, binEdges] = binVectorsFromDistance(dist, bin, mode);
clearvars bin mode
```

Smaller `bin` = finer FDM but more pixels. 0.07 nm is the workflow default.

### 2a. Plain hit-density FDM (no ranging needed)

Fastest option — a 2D histogram of all detector hits. `from`/`to` select an ion-index window, which corresponds to a slice of the experiment in time (lower indices = earlier in the run):

```matlab
from = 1;
to   = height(pos);          % use all hits; narrow this to look at a time window
FDM  = hist3([pos.detx(from:to), pos.dety(from:to)], binCenters);
figure; imagesc(FDM); axis equal;
clearvars from to
```

### 2b. Species concentration FDM (requires allocated, ranged data)

Voxelise the detector plane and compute the local concentration in each pixel, then display one species. `$ARGUMENTS` is the species to map (e.g. `Al`):

```matlab
startIdx = 1;
endIdx   = height(pos);
pix = posNdBin(pos(startIdx:endIdx,:), dist(startIdx:endIdx,:), binEdges);

detEff          = 82;                       % detector efficiency (% or decimal both accepted)
excludeListConc = {'unranged', 'Ga'};       % species to ignore (optional)
volumeName      = '';                       % volume name, optional
concentrationKernel = @(p) posCalculateConcentrationSimple(p, detEff, excludeListConc, volumeName);
conc = cellfun(concentrationKernel, pix, 'UniformOutput', false);

wantedSpecies = cellfun(@(p) p.Al(p.format == 'concentration'), conc);   % replace Al with $ARGUMENTS
figure;
imagesc(binCenters{1}, binCenters{2}, wantedSpecies);
axis equal;
```

Confirm the detection efficiency with the user (depends on instrument) and which species column to extract. Bin size strongly affects compute time here.

### 3. (Optional) FDM movie over the experiment

Stack frames sampled across the run into an `.avi`. Uses all species (hit density). Plays at 30 fps:

```matlab
frames   = 300;        % number of frames (movie length = frames/30 s)
sample   = 10000;      % ion hits per frame
RES      = 100;        % FDM resolution RES x RES (>1)
fileName = "FDM";      % output .avi saved in current folder
movieFDM = movieCreateFieldDesorptionMap(pos.detx, pos.dety, frames, sample, RES, fileName);
```

`frames`, `sample`, `RES`, `fileName` are all optional in the function. The movie reveals grain boundaries, poles, and phases "moving" as the tip is consumed.

## Interactive cropping (manual GUI action)

Spatial or temporal cropping is done by hand:
- **Temporal**: change the `from`/`to` ion-index window (step 2a) or the frame sampling (step 3) to focus on a part of the run.
- **Interactive spatial crop**: the toolbox function `plotFieldDesorptionMap(pos)` opens a figure with an interactive circular spatial crop. Tell the user to drag the circle in the MATLAB figure to keep only the hits inside it — do not try to automate the mouse interaction.

## Report

Tell the user which FDM was produced (plain hit-density vs. species concentration), the bin size/resolution used, and where any movie file was saved. Note that the FDM axis units are arbitrary for the `hist3` version and scale with bin width. Suggest:
- Changing `bin`/`RES` for more or less detail.
- Mapping a different species (step 2b) or a different time window (step 2a `from`/`to`).
- `plotFieldDesorptionMap(pos)` for interactive inspection and cropping.

## Toolbox functions used

`binVectorsFromDistance`, `posNdBin`, `posCalculateConcentrationSimple`, `movieCreateFieldDesorptionMap`, `plotFieldDesorptionMap` (optional, interactive); `hist3`, `cellfun`, `imagesc` (MATLAB built-ins).
