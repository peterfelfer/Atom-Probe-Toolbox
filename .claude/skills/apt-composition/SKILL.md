---
name: apt-composition
description: Calculate the composition (concentration) of atom probe data. Handles range extraction, ion allocation, and concentration calculation. Use when the user asks about composition, concentration, or chemical analysis.
argument-hint: "[instrument]"
---

Calculate the bulk composition from the current mass spectrum and pos data. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data (from `/apt-load`)
- `spec` — mass spectrum figure with ions and ranges defined (from `/apt-spectrum` + `rangeAddAll`)

If ranges haven't been defined yet, tell the user they need to define ranges first:
```matlab
rangeAddAll(spec, colorScheme, 0.04, true);
```

## Steps

### 1. Extract ranges from the spectrum

```matlab
rangeTable = rangesExtractFromMassSpec(spec);
disp(rangeTable);
```

### 2. Allocate ions

```matlab
pos = posAllocateRange(pos, rangeTable, 'decomposed');
```

Use `'decomposed'` mode to break molecular ions into constituent atoms. If the user specifically wants molecular ion identities, use `'direct'` instead.

### 3. Get detection efficiency

If the user provided an instrument as `$ARGUMENTS`, map it to detection efficiency:
- LEAP4000XHR → 0.37
- LEAP5000XR → 0.52
- LEAP5000XS → 0.80
- EIKOS → 0.50

Otherwise, read it from the config:
```matlab
cfg = APTConfig.getInstance();
detEff = cfg.reconstruction.detectionEfficiency;
```

If the detection efficiency is still the default (0.37) and the user hasn't confirmed their instrument, ask which instrument was used.

### 4. Calculate concentration

```matlab
conc = posCalculateConcentrationSimple(pos, detEff, {'unranged'}, '', 'mode', 'atomic');
disp(conc);
```

Exclude `'unranged'` atoms by default. If the user wants to include them or exclude other species, adjust the exclude list.

### 5. Calculate uncertainty

```matlab
unc = concentrationUncertainty(conc.counts);
```

### 6. Report

Present results as a formatted table showing:
- Each element/ion
- Concentration (at.%)
- Uncertainty (+/- at.%)
- Number of counts

Ask if the user wants to:
- Export results: `writetable(conc, 'composition.csv')`
- See isotopic or ionic concentration instead of atomic
- Exclude additional species
