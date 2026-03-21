# aptDataCalibration

Calibration functions for raw atom probe data from the [pyccapt](https://github.com/mmonajem/pyccapt) instrument control system. These functions bridge pyccapt's raw HDF5 data format to the Atom Probe Toolbox, providing time-of-flight conversion, voltage correction, bowl correction, and flight-path estimation.

After calibration, the resulting pos table integrates seamlessly with all standard toolbox functions (`massSpecPlot`, `ionAdd`, `rangeAdd`, `posReconstruct3DGeiser`, `scatterPlotPosWidget`, etc.).

---

## Contents

- [Quick Start](#quick-start)
- [Functions](#functions)
  - [Data Loading](#data-loading)
  - [TOF Conversion](#tof-conversion)
  - [Calibration](#calibration)
  - [Interactive Tools](#interactive-tools)
- [Workflows](#workflows)
- [pyccapt HDF5 Format](#pyccapt-hdf5-format)
- [Integration with the Toolbox](#integration-with-the-toolbox)

---

## Quick Start

```matlab
% Setup
setupToolbox;
load('isotopeTable_naturalAbundances.mat');
load('colorScheme.mat');

% Load pyccapt raw data (auto-detected via posLoad)
pos = posLoad('measurement.h5');

% Convert TOF to mass-to-charge
pos.mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, 110, 38, ...
    'mode', 'voltage', 'VP', pos.VP);

% Calibrate
pos = voltageCorrection(pos, [26.5 27.5]);
pos = bowlCorrection(pos, [26.5 27.5]);

% Continue with standard toolbox workflow
spec = massSpecPlot(pos.mc, 0.01, 'normalised');
ionAdd(spec, 'Al', 1, isotopeTable, colorScheme);
```

---

## Functions

### Data Loading

| Function | Description |
|----------|-------------|
| `posLoadPyccapt(fileName)` | Load pyccapt raw HDF5 (`/dld/*` datasets) into a standard pos table. Detector coordinates are converted from cm to mm. Called automatically by `posLoad` when a pyccapt HDF5 file is detected. |

```matlab
% Direct loading
pos = posLoadPyccapt('measurement.h5');

% Or via the canonical loader (auto-detects pyccapt format)
pos = posLoad('measurement.h5');

% With filtering
pos = posLoadPyccapt('measurement.h5', 'maxTof', 5000, 'minTof', 100);
```

**Output table columns:**

| Column | Unit | Source |
|--------|------|--------|
| `ionIdx` | 1 | Sequential index |
| `tof` | ns | `/dld/t` |
| `VDC` | V | `/dld/high_voltage` |
| `VP` | V | `/dld/pulse` (or legacy keys) |
| `detx` | mm | `/dld/x` × 10 |
| `dety` | mm | `/dld/y` × 10 |
| `laserIntensity` | pJ | `/dld/laser_intensity` (if present) |
| `startCounter` | 1 | `/dld/start_counter` (if present) |

---

### TOF Conversion

| Function | Description |
|----------|-------------|
| `tofToMassToCharge(tof, VDC, detx, dety, L, t0)` | Convert time-of-flight to mass-to-charge ratio (Da). Supports laser and voltage pulse modes. |
| `massToChargeToTof(mc, VDC, detx, dety, L)` | Inverse conversion: mass-to-charge to time-of-flight (ns). |

**Laser mode** (default):

$$mc = \frac{2 \cdot V_{DC} \cdot e}{amu} \left(\frac{t - t_0}{L}\right)^2$$

**Voltage mode** (with pulse correction):

$$mc = \frac{2 \cdot \alpha \cdot (V_{DC} + \beta \cdot V_P) \cdot e}{amu} \left(\frac{t - t_0}{L}\right)^2$$

where α = 1.015 (amplification factor) and β = 0.7 (averaging factor).

```matlab
% Laser mode
mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, 110, 38);

% Voltage mode
mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, 110, 38, ...
    'mode', 'voltage', 'VP', pos.VP);

% Custom alpha/beta
mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, 110, 38, ...
    'mode', 'voltage', 'VP', pos.VP, 'alpha', 1.02, 'beta', 0.65);
```

---

### Calibration

| Function | Description |
|----------|-------------|
| `voltageCorrection(pos, peakRange)` | Fit quadratic *f(V) = a + bV + cV²* to normalised peak position vs. voltage, then correct: *mc = mc / √f(V)*. |
| `bowlCorrection(pos, peakRange)` | Fit 2-D quadratic surface *f(x,y) = a + bx + cy + dx² + exy + fy²* to normalised peak position across the detector, then correct: *mc = mc / f(x,y)*. |

```matlab
% Voltage correction with default options
[pos, vFit] = voltageCorrection(pos, [26.5 27.5]);

% Bowl correction with finer grid
[pos, bFit] = bowlCorrection(pos, [26.5 27.5], 'gridSize', 3);

% Iterative calibration (recommended)
for iter = 1:10
    pos = voltageCorrection(pos, [26.5 27.5], 'showPlot', false);
    pos = bowlCorrection(pos, [26.5 27.5], 'showPlot', false);
end
```

**Voltage correction options:**

| Option | Default | Description |
|--------|---------|-------------|
| `'sampleSize'` | 1000 | Ions per group (ionSequence mode) or voltage step (voltage mode) |
| `'mode'` | `'ionSequence'` | Grouping mode: `'ionSequence'` or `'voltage'` |
| `'binWidth'` | 0.01 | Histogram bin width (Da) for peak finding |
| `'sampleMethod'` | `'histogram'` | Peak location method: `'histogram'`, `'mean'`, or `'median'` |
| `'showPlot'` | true | Show diagnostic fit plot |

**Bowl correction options:**

| Option | Default | Description |
|--------|---------|-------------|
| `'gridSize'` | 5 | Detector grid cell size (mm) |
| `'binWidth'` | 0.01 | Histogram bin width (Da) |
| `'sampleMethod'` | `'histogram'` | Peak location method |
| `'showPlot'` | true | Show diagnostic surface plot |

---

### Interactive Tools

| Function | Description |
|----------|-------------|
| `fineTuneT0(pos, L, t0)` | Interactive slider to adjust t0 while watching the mass spectrum update in real time. |
| `determineLAndT0(pos, peaks)` | Estimate flight-path length *L* and propagation delay *t0* via linear regression on known isotope peaks. |
| `plotExperimentHistory(pos)` | 2-D histogram (ion index vs. TOF) with voltage overlay. Draw a rectangle to crop temporally. |
| `plotFieldDesorptionMap(pos)` | 2-D histogram of detector coordinates. Draw a circle to crop spatially. |

```matlab
% Fine-tune t0 interactively
[pos, t0] = fineTuneT0(pos, 110, 38, 'mode', 'voltage', 'VP', pos.VP);

% Determine L and t0 from known peaks
selectedPeaks = table([26.7; 13.2], [27.3; 13.8], [26.9815; 13.4908], ...
    'VariableNames', {'mcLow', 'mcHigh', 'idealMc'});
[L, t0, fitInfo] = determineLAndT0(pos, selectedPeaks);

% Temporal crop
pos = plotExperimentHistory(pos);

% Spatial crop
pos = plotFieldDesorptionMap(pos);
```

---

## Workflows

Two example workflow scripts are included:

| Script | Description |
|--------|-------------|
| `Workflow_pyccapt_DataProcessing.m` | Full end-to-end workflow: load → crop → calibrate → identify ions → reconstruct → visualise. Uses calibration functions from this folder, then hands off to standard toolbox functions. |
| `Workflow_pyccapt_LAndT0Determination.m` | Estimate propagation delay and flight-path length from known isotope peaks using linear regression. |

Open in the MATLAB editor and run section-by-section (`Ctrl+Enter`).

---

## pyccapt HDF5 Format

pyccapt raw HDF5 files store detector-level data under a `/dld` group:

```
/dld/
├── t                 → time-of-flight (ns)
├── high_voltage      → DC specimen voltage (V)
├── pulse             → pulse voltage (V)  [or /dld/voltage_pulse, /dld/pulse_voltage]
├── laser_intensity   → laser pulse energy (pJ, optional)
├── x                 → detector x position (cm)
├── y                 → detector y position (cm)
└── start_counter     → pulse counter (optional)
```

`posLoadPyccapt` handles legacy key variations automatically and converts detector coordinates from cm to the toolbox standard of mm.

---

## Integration with the Toolbox

After loading and calibrating pyccapt data, the resulting pos table is a standard toolbox table. All downstream functions work without modification:

```matlab
% Load and calibrate
pos = posLoad('measurement.h5');
pos.mc = tofToMassToCharge(pos.tof, pos.VDC, pos.detx, pos.dety, 110, 38);
pos = voltageCorrection(pos, [26.5 27.5], 'showPlot', false);
pos = bowlCorrection(pos, [26.5 27.5], 'showPlot', false);

% Standard toolbox workflow from here
spec = massSpecPlot(pos.mc, 0.01, 'normalised');
ionAdd(spec, 'Al', 1, isotopeTable, colorScheme);
rangeAddAll(spec, colorScheme, 3, true);
rangeTable = rangesExtractFromMassSpec(spec);
pos = posAllocateRange(pos, rangeTable, 'decompose');
conc = posCalculateConcentrationSimple(pos, 0.57, {'unranged'});
scatterPlotPosWidget(pos, colorScheme);
```

---

*Part of the [Atom Probe Toolbox](../README.md). (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg*
