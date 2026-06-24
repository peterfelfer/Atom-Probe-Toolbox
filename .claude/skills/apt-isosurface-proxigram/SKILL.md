---
name: apt-isosurface-proxigram
description: Full pipeline from voxelisation through proxigram — voxelise, build an isoconcentration surface, then compute single- and multi-ion proxigrams across it in one guided sequence. Use when the user wants a proxigram but does not yet have an interface.
argument-hint: "[isovalue]"
---

Run the complete interface-analysis pipeline: voxelise the reconstruction, build an isosurface, and compute proxigrams across it. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — ranged/allocated atom probe data (has an `ion` column, from `/apt-composition`)
- `colorScheme` — color scheme (from `/apt-setup`), for the multi-ion plot

If ions have not been allocated yet, tell the user to run `/apt-composition` first. If the user already has an interface `fv` in the workspace, they can skip straight to `/apt-proxigram` instead.

**Data-type note:** if `conc`/`concMap` or the `pos` columns are single precision and `isosurface` or a proxigram function throws a type error, cast with `double()`.

## Steps

### 1. Voxelisation

Build a 3D binning grid and voxelise all atoms. `voxBin` is the voxel size in nm; mode is always `'distance'` for 3D voxelisation.

```matlab
dist   = [pos.x pos.y pos.z];
voxBin = [2 2 2];          % voxel size in nm (1 = fine, 2 = general, 3 = fast)
mode   = 'distance';
[binCenters, binEdges] = binVectorsFromDistance(dist, voxBin, mode);
gridVec = binCenters;
vox = posToVoxel(pos, gridVec);   % counts of all atoms per voxel
```

### 2. Concentration map

Voxelise only the species of interest and divide by the total to get the per-voxel concentration (a fraction, 0..1).

```matlab
species = {'Al'};                          % species defining the isosurface; ask the user
voxIon  = posToVoxel(pos, gridVec, species);
concMap = voxIon ./ vox;                    % concentration field (fraction)
```

### 3. Create the isosurface

The `isovalue` is in at.% (e.g. 8 = 8 at.%). It is common to try several values to best capture the feature.

```matlab
isovalue = 8;   % at.% threshold (typical exploration range 1-40)
fv = isosurface(gridVec{2}, gridVec{1}, gridVec{3}, concMap, isovalue);
figure;
p = patch(fv, 'FaceColor', [0.8 0.8 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
axis equal; rotate3d on; light;
axisSpatialAptify;
title(sprintf('Isosurface at %g at.%%', isovalue));
```

Use `$ARGUMENTS` for the isovalue if the user gave one; otherwise start at `8` and suggest trying others.

**Critical axis-ordering note:** MATLAB's `isosurface` treats the first array dimension as Y and the second as X, so grid vectors are passed as `gridVec{2}, gridVec{1}, gridVec{3}` (X and Y swapped). Do not "fix" this to `{1},{2},{3}` — if the surface looks mirrored, this swap is what keeps it aligned with the data.

### 4. Single-ion proxigram

`patchCreateProxigram(posSpecies, pos, fv, bin)` returns `[proxi, binVector]`: `proxi` is the concentration (fraction) and `binVector` the distance from the interface in nm. `bin` is the distance bin width in nm.

```matlab
proxiSpecies = pos(pos.ion == 'Al', :);   % species to profile
proxiBin     = 0.5;                        % bin width in nm (0.25 = fine, 0.5 = default)
[proxi, binVector] = patchCreateProxigram(proxiSpecies, pos, fv, proxiBin);

figure;
plot(binVector, proxi * 100, 'LineWidth', 2);
xlabel('distance from interface [nm]');
ylabel('concentration [%]');
title('Single-ion proxigram');
grid on;
```

### 5. Multi-ion proxigram

The `binVector` is the same for every species, so all are collected into one table.

```matlab
ionList  = {'Al', 'Ni', 'Cr', 'Co'};   % ions to include; ask the user
proxiBin = 0.5;
proxiAll = table();

for i = 1:length(ionList)
    posSpecies = pos(pos.ion == ionList{i}, :);
    [proxi, binVector] = patchCreateProxigram(posSpecies, pos, fv, proxiBin);
    proxi = proxi * 100;   % convert to at.%
    if i == 1
        proxiAll = addvars(proxiAll, binVector', proxi', 'NewVariableNames', {'binVector', ionList{i}});
    else
        proxiAll = addvars(proxiAll, proxi', 'NewVariableNames', ionList(i));
    end
end
```

### 6. Plot the multi-ion proxigram

```matlab
figure;
lineWidth = 2;
for k = 1:length(ionList)
    idxElement = find(string(proxiAll.Properties.VariableNames) == ionList{k});
    plot(proxiAll.binVector, table2array(proxiAll(:, idxElement)), ...
        'DisplayName', ionList{k}, 'LineWidth', lineWidth, ...
        'Color', colorScheme.color(colorScheme.ion == ionList{k}, :));
    hold on
end
hold off
legend('Location', 'best');
xlabel('distance from interface [nm]');
ylabel('concentration [%]');
title('Multi-ion proxigram');
grid on;
```

## Report

Tell the user:
- The voxel size, isosurface species + isovalue, and proxigram bin width used.
- That `fv` (interface), `concMap` (concentration field), and `proxiAll` (one column per ion in at.%, plus `binVector`) are in the workspace.
- Suggest trying other isovalues; the interface `fv` persists, so they can re-run only the proxigram steps (or `/apt-proxigram`) without re-voxelising.

Suggested follow-ups: `/apt-interfacial-excess` (quantitative segregation across this interface) or `patchToObj(fv, 'isosurface', 'fileName.obj')` to export the surface. Export profiles with `writetable(proxiAll, 'proxigram.csv');`.
