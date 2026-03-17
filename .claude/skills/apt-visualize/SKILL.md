---
name: apt-visualize
description: Create a 3D scatter plot visualization of atom probe data. Use when the user wants to see a 3D view, atom map, or point cloud of their APT data.
argument-hint: "[ions to show]"
---

Create a 3D scatter plot of the atom probe data. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data, ideally with ions allocated (has `ion` column)
- `colorScheme` — color scheme for ion coloring

## Visualization

If `pos` has an `ion` column (ranges have been allocated):

```matlab
scatterPlotPosWidget(pos, colorScheme);
```

This opens an interactive widget where the user can toggle ion visibility, rotate, zoom, and adjust transparency.

If `pos` does NOT have an `ion` column yet (no ranges allocated), use the basic scatter plot:

```matlab
scatterPlotPosData(pos, colorScheme);
```

## Notes

Tell the user:
- The 3D scatter plot is displayed in a MATLAB figure window
- They can interact with it directly in MATLAB (rotate, zoom, toggle ions)
- This is a GUI-based tool — interaction happens in MATLAB, not through Claude
- If they want to filter to specific ions, they can use the widget controls in the MATLAB figure
