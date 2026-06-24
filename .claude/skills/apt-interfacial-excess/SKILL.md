---
name: apt-interfacial-excess
description: Quantify interfacial excess (Gibbsian excess of solute, in at/nm^2) across an interface — either a single value via the Krakauer-Seidman integral method, or a 2D excess map over the surface. Use when the user wants the amount of an element segregated to an interface or grain boundary.
argument-hint: "[species]"
---

Quantify the interfacial excess of a species at an interface — the excess number of solute atoms associated with the interface, in at/nm^2. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — ranged/allocated atom probe data with `atom`/`ion` columns (from `/apt-composition`)
- `fv` — an interface struct with `faces` and `vertices` fields (an isosurface or grain boundary, from `/apt-isosurface`)

If there is no interface yet, tell the user to run `/apt-isosurface` first. For best statistics the parent `pos` should exclude unranged atoms. The integral method is well suited to quantifying a few hundred atoms segregated at the interface, so reasonable counting statistics matter.

**Data-type note:** if the `pos` columns are single precision and a function throws a type error, cast the coordinates with `double()`.

## Steps

### 1. Single interfacial-excess value (interactive)

`patchCreateExcessValue(posSpecies, parentPos, interface, vertices)` builds the Krakauer-Seidman integral curve and reports a single excess value. It has **no return value** — the result is displayed directly on the figure (below the integral curve). `posSpecies` is the subset of `pos` for the segregating species, `parentPos` is the full parent table (ideally unranged excluded), `interface` is the `fv` struct, and `vertices` optionally restricts the calculation to a sub-region.

```matlab
posSpecies = pos(pos.atom == 'B', :);   % segregating species (atomic symbol)
interface  = fv;                         % interface struct (faces + vertices)
vertices   = [];                         % [] = use the whole interface; else an Nx3 vertex subset
patchCreateExcessValue(posSpecies, pos, interface, vertices);
```

Use `$ARGUMENTS` as the species if the user gave one; otherwise ask which element segregates.

**Manual GUI step (do not automate):** the figure has interactive controls. Tell the user to use the *change Limits* button to drag/define the two sections used for the linear regression (the bulk slopes on either side of the interface), then the *fig/save* button to store the result as a figure. The displayed excess updates as they adjust the limits.

To restrict the calculation to a region away from the dataset edges (better statistics), pass a subset of `fv.vertices` as the `vertices` argument instead of `[]`.

### 2. (Alternative) 2D interfacial-excess map over the surface

`patchCreateInterfacialExcessMap(pos, interface, lim, deloc)` returns `IEmap` (excess per vertex, in at/nm^2) and renders it as a colored `trisurf` over the interface; it also prompts to save a `.ply`. Important: here `pos` is a **raw Nx3 coordinate matrix of the species atoms** (the function uses `pos(:,1:3)`), not a table. `lim` is the distance window in nm to count atoms near the interface, and `deloc` is an optional number of delocalisation (smoothing) passes.

```matlab
posSpecies = pos(pos.atom == 'B', :);
coords     = [posSpecies.x posSpecies.y posSpecies.z];   % Nx3 matrix, not a table
lim        = 1;     % count atoms within this distance (nm) of the interface
deloc      = 1;     % optional smoothing passes (omit for none)
IEmap = patchCreateInterfacialExcessMap(coords, fv, lim, deloc);
```

A `uiputfile` dialog will appear asking where to save the normalized map as `.ply`; tell the user they can cancel it to skip the export. The rendered map uses the same X/Y-swapped axis convention as the isosurface (it plots `vertices(:,2)` against `vertices(:,1)`).

## Report

Tell the user:
- The single excess value (read from the figure, in at/nm^2) and the species + interface used; remind them it depends on where they placed the regression limits.
- For the map: the `lim` and `deloc` used, that `IEmap` (at/nm^2 per vertex) is in the workspace, and the min/max range printed by the function.

Suggested follow-ups: `/apt-proxigram` for the full concentration-vs-distance profile across the same interface, or `patchToPly`/the built-in `.ply` export to view the excess map externally.
