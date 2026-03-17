---
name: apt-setup
description: Initialize a MATLAB session for atom probe analysis. Run this at the start of every APT session to set up the toolbox, load isotope tables, and load the default color scheme.
argument-hint: "[instrument]"
---

Initialize the MATLAB session for atom probe analysis. Use the MATLAB MCP server to execute the following setup code:

```matlab
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme_default.mat');
cfg = APTConfig.getInstance();
```

If the user provided an instrument name as `$ARGUMENTS`, also apply the instrument preset. Map common names to preset names:
- "LEAP 4000", "4000", "LEAP4000" → `cfg.applyInstrumentPreset('LEAP4000XHR')`
- "LEAP 5000 XR", "5000XR", "LEAP5000XR" → `cfg.applyInstrumentPreset('LEAP5000XR')`
- "LEAP 5000 XS", "5000XS", "LEAP5000XS" → `cfg.applyInstrumentPreset('LEAP5000XS')`
- "EIKOS", "eikos" → `cfg.applyInstrumentPreset('EIKOS')`

If no instrument is provided, ask the user which instrument their data was acquired on.

After setup, display the current configuration:
```matlab
cfg.display();
```

Report back to the user what was initialized and the active instrument settings (detection efficiency, flight path, t0).
