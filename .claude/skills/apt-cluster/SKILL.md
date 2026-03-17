---
name: apt-cluster
description: Run DBSCAN cluster analysis on atom probe data. Use when the user wants to find clusters, precipitates, or perform cluster analysis.
argument-hint: "[epsilon minPts]"
---

Run DBSCAN cluster analysis on the atom probe data. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data with ions allocated (has `ion` column)

If ions haven't been allocated, tell the user to run `/apt-composition` first.

## Parameters

Parse `$ARGUMENTS` for epsilon and minPts values. If not provided, use defaults:
- `epsilon` = 0.5 nm (neighborhood radius)
- `minPts` = 10 (minimum atoms in cluster)

If the user hasn't specified parameters, ask what they expect:
- For fine precipitates: epsilon ~0.3-0.5 nm, minPts ~5-15
- For larger precipitates: epsilon ~0.5-1.0 nm, minPts ~10-30
- These depend on the detection efficiency and atomic density

## Execution

### 1. Run DBSCAN

```matlab
[clusterIdx, clusterTable] = clusterDBSCAN(pos, epsilon, minPts);
```

Where `epsilon` and `minPts` come from the arguments or defaults.

### 2. Analyze results

```matlab
clusterSizeAnalyse;
```

### 3. Report

Tell the user:
- Number of clusters found
- Size distribution (if available)
- Suggest adjusting epsilon/minPts if results don't look right
- Remind them that cluster analysis results depend heavily on parameter choice
- Suggest using `voronoiVolumeAnalysis` as an alternative approach

## Notes

- DBSCAN is GPU-accelerated when available
- For large datasets (>1M atoms), consider using an ROI first to reduce computation time
- The user may want to filter by specific ion species before clustering
