---
name: apt-cluster-voronoi
description: Voronoi-cell cluster determination (Felfer et al. 2015) — detect and extract clusters via Voronoi-cell volumes, Delaunay triangulation, and a Kolmogorov-Smirnov test against a randomised reference. Use when the user wants the Voronoi/clusterDetermination method specifically (this is NOT DBSCAN — for density-based DBSCAN use /apt-cluster instead).
argument-hint: "[ion]"
---

Determine clusters using the Voronoi-cell volume method (`clusterDetermination`) from Felfer, Ceguerra, Ringer, Cairney, Ultramicroscopy 150 (2015) 30-36. Use the MATLAB MCP server. This is distinct from the DBSCAN approach in `/apt-cluster`.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — a DECOMPOSED, ranged pos file (has an `ion` column; from `/apt-composition`)
- `colorScheme` — from `/apt-setup` (only needed for the optional visualisation)

If ions are not allocated, run `/apt-composition` first.

**Data types:** these geometry functions need double-precision coordinates. If you hit a type error, wrap arguments in `double()`, e.g. `clusterDetermination(double(clusterPos), double(pos))`.

## Steps

### 1. Select the clustering species and run cluster determination

```matlab
clusterPos = pos(pos.ion == 'D H', :);   % choose the solute ion — from $ARGUMENTS or ask the user
[clusterParameter, clusteredAtoms, clusterPosVol, randomVolumes] = clusterDetermination(clusterPos, pos);
```

`clusterDetermination(clusterPos, pos, Nmin)` returns four tables:
- `clusterParameter` — K-S test result (1 = passed, 0 = failed), cluster percentage, `Nmin`, and the cluster cutoff volume
- `clusteredAtoms` — pos rows of the atoms assigned to clusters
- `clusterPosVol` — `clusterPos` plus a Voronoi-cell volume per cell
- `randomVolumes` — randomised reference volumes used for the test

`Nmin` (minimum atoms to count as a cluster) is determined automatically; pass a third argument to override.

### 2. (Optional) Combine runs across several species

```matlab
clusterElement = "D H";
clusterParameter = addvars(clusterParameter, clusterElement);
% clusterParaSample = clusterParameter;            % run after the FIRST determination
clusterParaSample = [clusterParaSample; clusterParameter];  % run after each later determination
```

### 3. (Optional) Step through the internals for more detail

```matlab
% Voronoi volume analysis
[numClustered, clusterCutoff, histCounts, experimentalVolumes, randomVolumes] = ...
    voronoiVolumeAnalysis(clusterPos, pos);

% Identify clusters from the volume threshold
volThreshORpos = clusterCutoff;
[clusterIdx, numClusters] = clusterIdentification(clusterPos, volThreshORpos, experimentalVolumes.expVol);
[randomClusterIdx, randomNumClusters] = clusterIdentification(randomVolumes, volThreshORpos, randomVolumes.randVol);

% Cluster size distribution -> Nmin
expClusterSizes = histcounts(clusterIdx, numClusters);
ranClusterSizes = histcounts(randomClusterIdx, randomNumClusters);
Nmin = clusterSizeAnalyse(expClusterSizes, ranClusterSizes);
```

(Confirmed signatures: `voronoiVolumeAnalysis(clusterPos, pos, vis, bins, Vmax)`; `clusterIdentification(clusterPos, volThreshORpos, volVoronoi, posDelaunay)`; `clusterSizeAnalyse(expClusterSizes, ranClusterSizes)`.)

### 4. Cluster composition

```matlab
detEff = 0.52;   % MUST match the instrument — confirm with the user
concClusters = posCalculateConcentrationSimple(clusteredAtoms, detEff, {'unranged'}, 'clusters', 'mode', 'atomic');
disp(concClusters);
```

Exclude `{'unranged'}` by default. Confirm the detection efficiency (LEAP4000XHR 0.37, LEAP5000XR 0.52, LEAP5000XS 0.80, EIKOS 0.50).

### 5. (Optional) Visualise the clustered atoms

```matlab
load colorScheme.mat   % if not already in the workspace
scatterPlotPosData(clusteredAtoms, colorScheme);
```

## Report

Tell the user:
- Whether the K-S test passed (clusters statistically significant vs. random)
- Cluster percentage, `Nmin`, and cluster cutoff volume from `clusterParameter`
- Number of clusters / clustered atoms and their composition
- That results depend on the chosen species and the random reference

Follow-ups:
- Export: `writetable(concClusters, 'clusterComposition.csv')`
- For the simpler, GPU-accelerated density-based alternative, use `/apt-cluster` (DBSCAN: `clusterDBSCAN(double(clusterPos), epsilon, minPts)`).
