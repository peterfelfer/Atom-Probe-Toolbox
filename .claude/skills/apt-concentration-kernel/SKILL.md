---
name: apt-concentration-kernel
description: Per-vertex concentration mapping onto a mesh/interface — compute atom-to-mesh distances, clip a slab around the interface, voxelise by nearest vertex, and apply a concentration kernel to colour each vertex. Use when the user wants a concentration map painted onto an isosurface or modelled grain boundary mesh.
---

Map the local concentration onto every vertex of a triangulated mesh (e.g. a modelled grain boundary or isosurface). Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — a DECOMPOSED, ranged pos file (allocated atoms; from `/apt-composition`)
- a mesh patch struct with `.faces` and `.vertices` — typically a modelled grain boundary `GB` (from the grain-boundary modelling workflow / isosurface). Below it is referred to as `GB6` for distances and `GB4` for display, matching the source workflow; substitute the user's actual variable name.
- `colorScheme` and instrument detection efficiency (see CLAUDE.md presets).

If `pos` is not decomposed/allocated, run `/apt-composition` first. If there is no mesh, the user must create one (isosurface or GB model) before running this.

## Steps

### 1. Distance of each atom to the mesh

```matlab
loc = posDistanceToMesh(pos, GB6);
```

`loc` is a table with each atom's signed `distance` to the mesh and its `closestVertex` index (distances are measured along vertex normals).

### 2. Clip to a slab around the interface

```matlab
clipDist = 2;   % nm — max distance from the mesh to keep
posGB = pos(abs(loc.distance) < clipDist, :);   % atoms within the slab
loc   = loc(abs(loc.distance) < clipDist, :);   % matching distances / closest-vertex indices
```

Pick `clipDist` to suit the feature width; 2 nm is the workflow default. (Note: `clipDist` must be numeric — do not quote it.)

### 3. Voxelise by nearest vertex

```matlab
vox = posNdBin(posGB, loc.closestVertex, {((1:length(GB6.vertices)+1) - .5)'});
```

`posNdBin` builds one sub-pos per vertex by binning atoms on their `closestVertex` index.

### 4. Apply the concentration kernel per vertex

```matlab
numVerts = length(vox);
ck = @(posGB) posCalculateConcentrationSimple(posGB, 0.37, {}, 'map');
conc = binApplyConcentrationKernel(vox, {1:numVerts}, ck, 'idx');
```

- The detection efficiency in `ck` (here `0.37`) MUST match the instrument — confirm with the user (LEAP4000XHR 0.37, LEAP5000XR 0.52, LEAP5000XS 0.80, EIKOS 0.50).
- The kernel computes one concentration per voxel/vertex. (Surface this default and ask if the user wants `{'unranged'}` excluded — the source workflow passes `{}`, i.e. no exclusions, for a per-vertex map.)

### 5. Display the concentration on the mesh

```matlab
figure
Cmap = conc(conc.format == 'concentration', :).C;
trisurf(GB4.faces, GB4.vertices(:,2), GB4.vertices(:,1), GB4.vertices(:,3), Cmap);
axisSpatialAptify
```

This colours each face/vertex of the mesh by its local concentration.

## Report

Tell the user: number of vertices mapped, the clip distance and detection efficiency used, and the concentration range across the mesh.

Follow-ups:
- Export the per-vertex concentration table: `writetable(conc, 'vertexConcentration.csv')`
- Choose a different species/mode by adjusting the `ck` handle (e.g. `'mode','atomic'`).
- For a profile across the interface instead of a surface map, see `/apt-proxigram`.
