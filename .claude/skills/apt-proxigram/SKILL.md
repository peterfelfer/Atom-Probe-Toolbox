---
name: apt-proxigram
description: Compute a proximity histogram (proxigram) — concentration vs. distance from an existing interface (isosurface, line, or point cloud) — for one or many species. Use when the user has an interface and wants a concentration profile across it.
argument-hint: "[species]"
---

Compute a proxigram (concentration vs. distance from an interface) for an interface that already exists in the workspace. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — ranged/allocated atom probe data (has an `ion` column, from `/apt-composition`)
- `fv` — an interface struct with `faces` and `vertices` fields (an isosurface, from `/apt-isosurface`)
- `colorScheme` — color scheme (from `/apt-setup`), for the multi-ion plot

If there is no interface yet, tell the user to run `/apt-isosurface` first (or `/apt-isosurface-proxigram` for the full voxelise → isosurface → proxigram pipeline).

**Data-type note:** if the `pos` columns are single precision and a proxigram function throws a type error, cast the coordinates with `double()` (e.g. build `posSpecies`/`pos` coordinate columns as `double(...)`).

## Steps

### 1. Single-ion proxigram (from an isosurface)

`patchCreateProxigram(posSpecies, pos, interface, bin)` returns `[proxi, binVector]`: `proxi` is the concentration (fraction, 0..1) and `binVector` is the distance from the interface in nm. `posSpecies` is the subset of `pos` for the species of interest, `pos` is the full parent table, `interface` is the `fv` struct, and `bin` is the distance bin width in nm.

```matlab
interface  = fv;                         % interface struct (faces + vertices)
proxiBin   = 0.5;                        % bin width in nm (0.25 = fine near interfaces, 0.5 = default)
posSpecies = pos(pos.ion == 'Al', :);    % the species of interest
[proxi, binVector] = patchCreateProxigram(posSpecies, pos, interface, proxiBin);

figure;
plot(binVector, proxi * 100, 'LineWidth', 2);
xlabel('distance from interface [nm]');
ylabel('concentration [%]');
title('Single-ion proxigram');
grid on;
```

Use `$ARGUMENTS` as the species if the user gave one; otherwise ask which species to profile. Note `proxi` is a fraction, so multiply by 100 for at.%.

### 2. Multi-ion proxigram

The `binVector` is identical for every species, so all concentrations can be collected into one table for comparison.

```matlab
species  = {'Cr', 'Ni', 'Al', 'Co'};   % ions to include; ask the user
proxiBin = 0.5;
proxiAll = table();

for i = 1:length(species)
    posSpecies = pos(pos.ion == species{1,i}, :);
    [proxi, binVector] = patchCreateProxigram(posSpecies, pos, interface, proxiBin);
    proxi = proxi * 100;   % convert to at.%
    if i == 1
        proxiAll = addvars(proxiAll, binVector', proxi', 'NewVariableNames', {'binVector', species{1,i}});
    else
        proxiAll = addvars(proxiAll, proxi', 'NewVariableNames', species{1,i});
    end
end
```

### 3. Plot the multi-ion proxigram

```matlab
figure;
lineWidth = 2;
for k = 1:length(species)
    idxElement = find(string(proxiAll.Properties.VariableNames) == species{1,k});
    plot(proxiAll.binVector, table2array(proxiAll(:, idxElement)), ...
        'DisplayName', species{1,k}, 'LineWidth', lineWidth, ...
        'Color', colorScheme.color(colorScheme.ion == species{1,k}, :));
    hold on
end
hold off
legend('Location', 'best');
xlabel('distance [nm]');
ylabel('concentration [%]');
grid on;
```

### 4. (Alternative) Proxigram from a line interface

If the interface is a line/edge object (e.g. a traced grain boundary) rather than an isosurface, use `lineCreateProxigram(posSpecies, pos, line, bin)`. The `line` argument is a struct with fields `vertices` (Nx3) and `edges` (Mx2); distances are measured perpendicular to the line segments. Same outputs as above.

```matlab
[proxi, binVector] = lineCreateProxigram(posSpecies, pos, line, proxiBin);
```

### 5. (Alternative) Proxigram from a point cloud / vertex set

If the reference is a set of vertices (no faces), use `pointCreateProxigram(posSpecies, pos, vertices, bin)`. `vertices` is an m-by-3 (or wider) matrix with m > 3; distances are measured to the nearest vertex. Same outputs as above.

```matlab
[proxi, binVector] = pointCreateProxigram(posSpecies, pos, vertices, proxiBin);
```

## Report

Tell the user:
- The bin width, species, and interface used.
- That `proxiAll` holds all profiles (column `binVector` plus one column per ion, in at.%).
- A proxigram is a function of distance, not a single number — for a quantitative segregation value across the interface use `/apt-interfacial-excess`.

Offer to export: `writetable(proxiAll, 'proxigram.csv');`
