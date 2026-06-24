---
name: apt-data-quality
description: Run quality-control checks on an atom probe dataset before analysis — basic stats, voltage curve, detector hit map (FDM), mass spectrum sanity, coordinate ranges, multi-hit fraction, and automated quality metrics. Use when the user wants to assess data quality or sanity-check a reconstruction.
---

Assess the quality and integrity of an APT dataset before committing to a full analysis. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data (from `/apt-load`)

For the full set of checks an **.epos** file is needed (it carries `VDC`, `detx`, `dety`, `multi`). A plain .pos file still supports the statistics, mass spectrum, coordinate, and metrics checks. Each .epos-only step below guards on the column being present.

If nothing is loaded, run:
```matlab
pos = posLoad;   % load .epos for full QC; .pos for the subset
```

## Steps

### 1. Basic statistics

```matlab
numAtoms = height(pos);
fprintf('Total atoms: %d\n', numAtoms);
fprintf('X range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.x), max(pos.x), max(pos.x)-min(pos.x));
fprintf('Y range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.y), max(pos.y), max(pos.y)-min(pos.y));
fprintf('Z range: %.1f to %.1f nm (span: %.1f nm)\n', min(pos.z), max(pos.z), max(pos.z)-min(pos.z));
fprintf('Mass-to-charge range: %.2f to %.2f Da\n', min(pos.mc), max(pos.mc));
```

### 2. Voltage curve (.epos)

VDC over the hit sequence shows measurement stability. A smooth, monotonically rising curve is good; jumps or drops suggest specimen fracture or laser/voltage drift.
```matlab
if ismember('VDC', pos.Properties.VariableNames)
    figure; plot(pos.VDC, 'LineWidth', 0.5);
    xlabel('ion hit index'); ylabel('V_{DC} [V]'); title('Voltage curve'); grid on;
else
    warning('VDC column not found. Load an .epos file for voltage curve analysis.');
end
```

### 3. Detector hit map / FDM (.epos)

The detector hit map reveals crystallographic poles, dead spots, and detector artefacts. Uniform coverage with visible poles is good; large dead regions or asymmetry are warning signs.
```matlab
if ismember('detx', pos.Properties.VariableNames)
    figure;
    FDM = hist3([pos.detx, pos.dety], [256 256]);
    imagesc(FDM'); axis equal; axis tight; colorbar;
    title('Detector hit map (FDM)'); xlabel('det x'); ylabel('det y');
else
    warning('Detector coordinates not found. Load an .epos file for FDM analysis.');
end
```

### 4. Mass spectrum overview

Coarse (0.1 Da) bins for overall shape and fine (0.01 Da) bins for peak detail, both on a log scale. Look for well-defined peaks, a reasonable background level, and no odd artefacts.
```matlab
figure;
subplot(2,1,1);
histogram(pos.mc, 'BinWidth', 0.1, 'EdgeColor', 'none'); set(gca, 'YScale', 'log');
xlabel('mass-to-charge [Da]'); ylabel('counts'); title('Mass spectrum — coarse (0.1 Da bins)');
xlim([0 min(200, max(pos.mc))]);
subplot(2,1,2);
histogram(pos.mc, 'BinWidth', 0.01, 'EdgeColor', 'none'); set(gca, 'YScale', 'log');
xlabel('mass-to-charge [Da]'); ylabel('counts'); title('Mass spectrum — fine (0.01 Da bins)');
xlim([0 min(100, max(pos.mc))]);
```

### 5. Coordinate outlier detection

Per-axis histograms expose atoms far outside the main point cloud (reconstruction artefacts / detector noise) that may need excluding.
```matlab
figure;
subplot(1,3,1); histogram(pos.x, 100); xlabel('x [nm]'); title('X distribution');
subplot(1,3,2); histogram(pos.y, 100); xlabel('y [nm]'); title('Y distribution');
subplot(1,3,3); histogram(pos.z, 100); xlabel('z [nm]'); title('Z distribution');
```

### 6. Multi-hit analysis (.epos)

A high multi-hit fraction signals detector pile-up that can bias composition.
```matlab
if ismember('multi', pos.Properties.VariableNames)
    multiHits = sum(pos.multi > 1);
    multiHitFraction = multiHits / height(pos) * 100;
    fprintf('Multi-hit events: %d (%.1f%% of total)\n', multiHits, multiHitFraction);
    figure; histogram(pos.multi, 'BinMethod', 'integers');
    xlabel('multiplicity'); ylabel('counts');
    title(sprintf('Multi-hit distribution (%.1f%% multi-hits)', multiHitFraction));
else
    warning('Multi column not found. Load an .epos file for multi-hit analysis.');
end
```

### 7. Automated quality metrics

`dataQualityMetrics` computes spatial-resolution estimates, density variation, detection efficiency, and other indicators in one call. `'showPlots', true` adds diagnostic plots.
```matlab
metrics = dataQualityMetrics(pos, 'showPlots', true);
disp(metrics.summary);
```

## Report

Summarise for the user:
- Atom count and coordinate spans.
- Voltage-curve stability (smooth vs. jumps), FDM coverage (poles vs. dead spots).
- Mass spectrum peak/background quality, any coordinate outliers.
- Multi-hit fraction (flag if high) and the `metrics.summary` verdict.
- If only a .pos file was loaded, note which checks (voltage, FDM, multi-hit) were skipped and that an .epos file enables them.

If the data looks clean, suggest proceeding to `/apt-spectrum` and `/apt-composition`.
