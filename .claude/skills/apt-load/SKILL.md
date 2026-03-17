---
name: apt-load
description: Load an atom probe data file (.pos or .epos) and display summary statistics (atom count, coordinate ranges, mass-to-charge range). Use when the user wants to load or open APT data.
argument-hint: "<file.pos|file.epos>"
---

Load an atom probe data file and show its summary. Use the MATLAB MCP server.

The file to load is: `$ARGUMENTS`

If no file is provided, ask the user for the file path. The file must be a `.pos` or `.epos` file.

Execute this MATLAB code (adapt the file path from the argument):

```matlab
pos = posLoad('$ARGUMENTS');

% Summary statistics
fprintf('Loaded: %s\n', '$ARGUMENTS');
fprintf('Atoms: %d\n', height(pos));
fprintf('Columns: %s\n', strjoin(pos.Properties.VariableNames, ', '));
fprintf('\nCoordinate ranges:\n');
fprintf('  x: [%.2f, %.2f] nm\n', min(pos.x), max(pos.x));
fprintf('  y: [%.2f, %.2f] nm\n', min(pos.y), max(pos.y));
fprintf('  z: [%.2f, %.2f] nm\n', min(pos.z), max(pos.z));
fprintf('  m/c: [%.2f, %.2f] Da\n', min(pos.mc), max(pos.mc));
```

Store the result in the workspace variable `pos` so subsequent commands can use it.

Report:
1. Number of atoms
2. File type (.pos vs .epos) and available columns
3. Spatial extent (x, y, z ranges in nm)
4. Mass-to-charge range
