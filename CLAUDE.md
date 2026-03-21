# Atom Probe Toolbox — Claude Code Guide

## Project Overview

The Atom Probe Toolbox is a MATLAB-based open-source toolkit for atom probe tomography (APT) data analysis. It covers the full APT pipeline: data import, mass spectrum analysis, ion identification, range definition, composition calculation, 3D visualization, spatial statistics, proxigram analysis, cluster analysis, and crystallography.

- **Language**: MATLAB R2025b
- **Required toolboxes**: Statistics and Machine Learning, Parallel Computing (optional), Image Processing (optional)
- **Author**: Prof. Peter Felfer Group, FAU Erlangen-Nurnberg
- **Users**: Materials scientists and atom probe researchers

## Setup

Always run at the start of a new MATLAB session:

```matlab
setupToolbox;
```

Most workflows also need:

```matlab
load('isotopeTable_naturalAbundances.mat');  % provides: isotopeTable
load('colorScheme.mat');                      % provides: colorScheme
```

## Instrument Presets

| Instrument     | Detection Eff. | Flight Path (mm) | t0 (ns) |
|----------------|---------------|-------------------|----------|
| LEAP 4000X HR  | 0.37          | 110               | 55.7     |
| LEAP 5000 XR   | 0.52          | 120               | 50.0     |
| LEAP 5000 XS   | 0.80          | 120               | 50.0     |
| EIKOS          | 0.50          | 100               | 45.0     |

Apply via APTConfig:

```matlab
cfg = APTConfig.getInstance();
cfg.applyInstrumentPreset('LEAP5000XR');
```

Or load from YAML config files in `configs/`:

```matlab
config = configYamlImport('configs/LEAP5000XR.yaml');
```

### Material Presets (atomic volumes in nm^3)

| Material | Atomic Volume |
|----------|--------------|
| Fe       | 0.01178      |
| Al       | 0.01660      |
| Cu       | 0.01182      |
| Ni       | 0.01094      |
| Ti       | 0.01766      |
| W        | 0.01583      |
| Si       | 0.02000      |
| Mg       | 0.02305      |

```matlab
cfg.applyMaterialPreset('Fe');
```

## Core Data Model

All data is stored in **MATLAB tables**.

### pos table columns

Basic (.pos): `ionIdx`, `x`, `y`, `z`, `mc`

Extended (.epos): adds `tof`, `VDC`, `VP`, `detx`, `dety`, `deltaP`, `multi`

After range allocation: adds `ion`, `charge`, `atom`

### Units

- Position: **nm**
- Mass-to-charge: **Da** (Dalton)
- Flight path: **mm**
- Time-of-flight: **ns**
- Atomic volume: **nm^3**

## Standard Analysis Workflow

### 1. Load data

```matlab
pos = posLoad('file.pos');       % or .epos
```

### 2. Mass spectrum

```matlab
spec = massSpecPlot(pos, 0.01, 'normalised');
```

### 3. Ion identification

```matlab
ionAdd(spec, 'Fe', 2, isotopeTable, colorScheme);   % adds Fe2+ markers
ionAdd(spec, 'C', 1, isotopeTable, colorScheme);     % adds C1+ markers
```

### 4. Range definition

```matlab
rangeAddAll(spec, colorScheme, 0.04, true);
```

Arguments: `(spec, colorScheme, margin, useMin)` — `margin` is the range margin in Da, `useMin` uses local minimum as range boundary.

### 5. Extract range table

```matlab
rangeTable = rangesExtractFromMassSpec(spec);
```

### 6. Allocate ions to atoms

```matlab
pos = posAllocateRange(pos, rangeTable, 'decomposed');
```

Mode `'decomposed'` decomposes molecular ions into constituent atoms. Mode `'direct'` keeps molecular ion identities.

### 7. Composition

```matlab
conc = posCalculateConcentrationSimple(pos, detEff, {}, '', 'mode', 'atomic');
```

Arguments: `(pos, detEff, excludeList, volumeName, 'mode', mode)`
- `detEff`: detection efficiency (from instrument preset)
- `excludeList`: cell array of ion names to exclude, e.g. `{'unranged'}`
- `mode`: `'atomic'`, `'isotopic'`, or `'ionic'`

### 8. Uncertainty

```matlab
unc = concentrationUncertainty(conc.counts);
```

### 9. 3D visualization

```matlab
scatterPlotPosWidget(pos, colorScheme);
```

### 10. ROI creation

```matlab
roiCreateBox;
roiCreateCylinder;
roiCreateSphere;
```

### 11. Proxigram

```matlab
lineCreateProxigram;
patchCreateProxigram;
```

### 12. Cluster analysis

```matlab
clusterDBSCAN(pos, epsilon, minPts);
```

## Function Reference

### Data I/O

- `posLoad(fileName)` — load .pos/.epos/.apt/.h5 file into table (auto-detects pyccapt HDF5)
- `posLoadPyccapt(fileName)` — load pyccapt raw HDF5 (`/dld/*` datasets) into standard pos table
- `posExport(pos, fileName)` — export pos table to file
- `posTableFromHDF5(fileName, datasetName)` — load from HDF5
- `posTableAddToHDF5(pos, fileName, datasetName)` — save to HDF5
- `rangesExtractFromFile(fileName)` — parse .rrng/.rng range file
- `stlToPatch(fileName)` — import STL mesh
- `patchToStl(patch, fileName)` — export STL mesh
- `patchToObj(patch, fileName)` — export OBJ mesh
- `patchToPly(patch, fileName)` — export PLY mesh

### Mass Spectrum

- `massSpecPlot(pos, binWidth, normMode)` — create mass spectrum plot; `normMode`: `'normalised'` or `'counts'`
- `massSpecAnalysis(pos)` — interactive mass spectrum analysis
- `massSpecComparison(specList)` — compare multiple spectra
- `backgroundEstimate(spec, method)` — fit background model
- `chargeStateRatioCalculate(spec)` — compute charge state ratios

### Ion & Range Management

- `ionAdd(spec, element, chargeState, isotopeTable, colorScheme)` — add ion markers to spectrum
- `ionConvertName(ionName)` — parse ion name string
- `ionConvertMode(ionTable, mode)` — convert between ionic/atomic modes
- `ionsCreateIsotopeList(element, isotopeTable)` — list isotopes for element
- `ionsCreateComplex(elements, isotopeTable)` — create complex ion isotope pattern
- `rangeAdd(spec, colorScheme)` — manually define a range on spectrum
- `rangeAddAll(spec, colorScheme, margin, useMin)` — add ranges for all identified ions
- `rangesExtractFromMassSpec(spec)` — extract range table from spectrum figure
- `rangesFromPos(pos)` — extract ranges from allocated pos table

### Allocation & Composition

- `posAllocateRange(pos, rangeTable, mode)` — allocate ions; `mode`: `'decomposed'` or `'direct'`
- `posCalculateConcentrationSimple(pos, detEff, excludeList, volumeName, 'mode', mode)` — bulk composition
- `posCalculateConcentrationBackgroundRemoved(...)` — composition with background subtraction
- `posCalculateConcentrationDeconvolved(...)` — composition with peak overlap deconvolution
- `concentrationUncertainty(counts)` — statistical counting uncertainty
- `binApplyConcentrationKernel(pos, kernelSize)` — local concentration kernel
- `calculate1Dprofile(pos, ...)` — 1D concentration profile
- `calculate2Dprofile(pos, ...)` — 2D concentration map

### ROI & Spatial Analysis

- `roiCreateBox` — interactive box ROI
- `roiCreateCylinder` — interactive cylinder ROI
- `roiCreateSphere` — interactive sphere ROI
- `roiCreatePlane` — interactive plane ROI
- `roiFromObj(fileName)` — ROI from OBJ mesh file
- `posInConvexHull(pos, hull)` — test points inside convex hull

### Proxigram

- `pointCreateProxigram` — proxigram from point
- `lineCreateProxigram` — proxigram from line/interface
- `patchCreateProxigram` — proxigram from isosurface
- `patchCreateInterfacialExcessMap` — interfacial excess on isosurface
- `patchCreateExcessValue` — total interfacial excess

### Cluster Analysis

- `clusterDBSCAN(pos, epsilon, minPts)` — DBSCAN clustering (GPU-accelerated)
- `clusterDetermination` — cluster analysis workflow
- `clusterIdentification` — identify cluster atoms
- `clusterSizeAnalyse` — cluster size distribution
- `voronoiVolumeAnalysis` — Voronoi cell volume analysis

### Spatial Statistics

- `spatialStatistics` — spatial distribution functions
- `spatialDistributionMap` — SDM computation
- `dataQualityMetrics` — reconstruction quality assessment

### Crystallography

- `crystalOrientation` — determine crystal orientation
- `ipfColor` — inverse pole figure coloring
- `ipfHistogram` — IPF histogram
- `ipfMesh` — IPF mesh
- `misorientation` — misorientation analysis
- `boundaryCharacter` — grain boundary character
- `stereoProj` — stereographic projection
- `plotOrientationCube` — orientation cube visualization

### Raw Data Calibration (aptDataCalibration/)

- `tofToMassToCharge(tof, VDC, detx, dety, L, t0, 'mode', 'voltage', 'VP', VP)` — TOF to mc conversion
- `massToChargeToTof(mc, VDC, detx, dety, L)` — mc to TOF (inverse)
- `voltageCorrection(pos, [lo hi], 'sampleSize', 1000)` — quadratic voltage-dependent mc correction
- `bowlCorrection(pos, [lo hi], 'gridSize', 5)` — 2-D quadratic spatial mc correction
- `fineTuneT0(pos, L, t0, 'mode', 'voltage', 'VP', VP)` — interactive t0 slider
- `determineLAndT0(pos, selectedPeaks, 'detectorWindow', 5)` — fit L and t0 via linear regression
- `plotExperimentHistory(pos)` — 2-D histogram (ion index vs TOF) with interactive temporal crop
- `plotFieldDesorptionMap(pos)` — detector FDM with interactive circular spatial crop

### Reconstruction

- `posReconstruct3DGeiser(pos, config)` — Geiser reconstruction
- `posReflectronCorrection(pos)` — reflectron time correction

### Visualization

- `scatterPlotPosData(pos, colorScheme)` — static 3D scatter plot
- `scatterPlotPosWidget(pos, colorScheme)` — interactive 3D scatter plot with GUI
- `colorSchemeCreate()` — create new color scheme
- `colorSchemeIonAdd(colorScheme, ionName, color)` — add ion to color scheme
- `movieCreateTurntableAnimation(...)` — create rotating 3D animation

### Configuration

- `APTConfig.getInstance()` — get singleton config object
- `APTConfig.get(key)` — get config value, e.g. `APTConfig.get('reconstruction.detectionEfficiency')`
- `APTConfig.set(key, value)` — set config value
- `cfg.applyInstrumentPreset('LEAP5000XR')` — apply instrument settings
- `cfg.applyMaterialPreset('Fe')` — apply material settings (atomic volume)
- `configYamlImport(fileName)` — load config struct from YAML
- `configYamlExport(config, fileName)` — save config struct to YAML
- `configExport(config, fileName)` — export config to JSON

### Utilities

- `setupToolbox` — add toolbox paths to MATLAB path
- `checkDependencies` — verify required toolboxes
- `batchProcess(fileList, processFcn)` — batch process multiple datasets
- `ProgressTracker` — progress bar utility
- `runTests` — run toolbox test suite

## Agentic Workflow Guidelines

When acting as an analysis agent:

1. **Always ask which instrument** the user's data was acquired on (or check if they've told you). The detection efficiency is critical for correct composition.
2. **Run `setupToolbox` first** in every new MATLAB session before calling any toolbox functions.
3. **Load lookup tables** (`isotopeTable_naturalAbundances.mat` and `colorScheme_default.mat`) before ion identification.
4. **For scripted analysis**: build a complete `.m` file, save it to the working directory, then execute via MATLAB MCP's `run_matlab_file` or use `evaluate_matlab_code` for short snippets.
5. **For interactive/GUI steps** (range picking with mouse, ROI manipulation, scatter plot interaction): explain to the user what to do in the MATLAB GUI — do not try to automate mouse interactions.
6. **Concentration calculation**: always ask about excluded species. Typically `{'unranged'}` should be excluded.
7. **Return results as tables**: use `disp(table)` or `writetable(table, 'results.csv')` to show output.
8. **Error handling**: if a function errors, read the error message, check function signatures in `functionSignatures.json`, and adjust arguments.
9. **Large datasets**: for datasets > 1M atoms, consider using chunked processing or ROIs to avoid memory issues.
