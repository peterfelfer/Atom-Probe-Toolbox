---
name: apt-colorscheme
description: Create a new ion color scheme or add ions/colors to the existing one. Produces/updates the `colorScheme` workspace variable used by ranging and 3D visualization everywhere else. Use when the user wants custom ion colors, a new color scheme, or to add a missing ion color.
argument-hint: "[ions]"
---

Create a new `colorScheme` or extend the existing one with new ions and colors. The `colorScheme` table maps each ion/molecule to an RGB color and is consumed by `ionAdd`, `rangeAddAll`, and the 3D scatter plots. Use the MATLAB MCP server.

## Prerequisites

- None required to create a fresh scheme.
- To extend the current scheme, `colorScheme` should already be in the workspace (from `/apt-setup`, or `load('colorScheme.mat')`).
- To harvest ion names from a spectrum (step 3), the workspace must contain `spec` (from `/apt-spectrum`).

This skill creates or updates the `colorScheme` workspace variable. Generated colors start in HSV (S=0.8, V=1.0, equidistant hue) and are converted to RGB; because human color perception is non-linear, manual tweaks may still be desirable.

## Steps

### 1. Build an ionTable, then create a new colorScheme

Provide the ions as a cell array (`$ARGUMENTS`, e.g. `Fe C H N O`). Either list them by hand or pull every ion already on a spectrum.

```matlab
ions = {'Fe'; 'C'; 'H'; 'N'; 'O'};      % from $ARGUMENTS
ionTable = table(categorical(ions));
ionTable.Properties.VariableNames{1} = 'ion';

colorScheme_new = colorSchemeCreate(ionTable)
```

Alternative source for `ionTable` (uses an existing ranged spectrum):
```matlab
ionTable = ionsExtractFromMassSpec(spec);
colorScheme_new = colorSchemeCreate(ionTable);
```

If the user wants this to become the active scheme, assign it: `colorScheme = colorSchemeCreate(ionTable);`

### 2. Add a single new ion to the colorScheme

```matlab
newIon = '1865Da';         % name of the ion/molecule to add
selection = 'create';      % 'create' = farthest-from-existing color; 'select' = pick interactively; omit = random
colorScheme = colorSchemeIonAdd(colorScheme, newIon, selection)
```

- `'create'` chooses the color most distinct from all existing ones (recommended for batch additions).
- `'select'` lets the user pick the color in a GUI color picker — explain this is a manual action, do not automate it.
- Omitting the third argument generates a random color.
- If the ion already exists, the function prints "ion already exists in colorScheme" and returns the scheme unchanged — this is safe to ignore.

### 3. Add every ion from a spectrum at once

Useful after adding many new ions to a spectrum that aren't yet in the scheme:

```matlab
ionTable = ionsExtractFromMassSpec(spec);
newName = char(ionTable.ionName);
for i = 1:size(newName)
    colorScheme = colorSchemeIonAdd(colorScheme, newName(i,:));
end
```

"ion already exists in colorScheme" messages during the loop are expected and harmless.

## Report

Tell the user the `colorScheme` variable was created/updated and how many ions it now contains (`height(colorScheme)`). Note that this variable is now picked up automatically by `ionAdd`, `rangeAddAll`, and the scatter-plot widgets.

Suggest follow-ups:
- Persist it for reuse: `save('colorScheme.mat', 'colorScheme')`
- Re-run `/apt-visualize` or `/apt-spectrum` to see the new colors applied.
- If two ions look too similar, re-run step 2 with `selection = 'select'` to hand-pick a color.

## Toolbox functions used

`colorSchemeCreate`, `colorSchemeIonAdd`, `ionsExtractFromMassSpec`.
