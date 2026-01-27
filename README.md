# Atom Probe Toolbox

**A comprehensive MATLAB toolbox for Atom Probe Tomography (APT) data analysis**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/Version-1.0-orange.svg)](https://github.com/peterfelfer/Atom-Probe-Toolbox)

Developed by the [Felfer Group](https://www.ww1.tf.fau.de/) at Friedrich-Alexander-Universität Erlangen-Nürnberg.

---

## Table of Contents
- [Versions](#versions)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [Function Reference](#function-reference)
- [Workflows](#workflows)
- [Contributing](#contributing)
- [Citation](#citation)
- [Authors](#authors)
- [License](#license)

---
## Versions
The top-level repository represents the current nightly builds. 
Older versions are in their respective subfolder. 

The status as of the publication of the following papers is contained in the Release 1.0

Heller et al., A MATLAB Toolbox for Findable, Accessible, Interoperable, and Reusable Atom Probe Data Science, Microscopy and Microanalysis (2025): [https://doi.org/10.1093/mam/ozae031](https://doi.org/10.1093/mam/ozae031)
Basic paper about the toolbox and functionality.

Heller et al., Compensating Image Distortions in a Commercial Reflectron-Type Atom Probe, Microscopy and Microanalysis (2025): [https://doi.org/10.1093/mam/ozae052](https://doi.org/10.1093/mam/ozae052)
This is relevant if you use the Reflectron distortion corrections for APT reconstruction. 

If you use this tool box for your atom probe tomography analysis, please consider citing the relevant papers. 

**Contributions and suggestions for improvement are always welcome! Happy Atom Probing**

## Features

### Data Import & Export
- Read `.pos`, `.epos`, and `.APT` (CAMECA) file formats
- HDF5 database creation with standardized metadata schema
- Export to `.obj`, `.ply`, `.h5` formats
- RRNG range file import

### Mass Spectrum Analysis
- Interactive mass spectrum plotting and ranging
- Automated peak detection and background fitting
- Isotope pattern matching
- Peak deconvolution
- Charge state ratio analysis

### 3D Visualization & ROI
- Interactive 3D point cloud visualization
- Region of Interest (ROI) selection: box, cylinder, sphere, plane
- ROI manipulation GUI (`roiManipulate.mlapp`)
- Custom color schemes per ion species

### Composition Analysis
- Local and global concentration calculation
- 1D and 2D concentration profiles
- Proxigram analysis (point, line, and surface-based)
- Concentration uncertainty quantification
- Interfacial excess calculations

### Cluster Analysis
- DBSCAN clustering with GPU acceleration (CPU fallback)
- Voronoi volume analysis
- Cluster size distribution analysis
- Spatial distribution functions (RDF, nearest neighbor)

### Crystallography
- Crystal orientation handling
- Inverse Pole Figure (IPF) coloring
- Grain boundary characterization
- Misorientation angle calculations
- CSL boundary identification (Sigma values)
- Stereographic projections

### Reconstruction & Correction
- Geiser reconstruction algorithm
- Reflectron distortion correction
- Field desorption map analysis

### Correlative Microscopy
- Volume resampling and rotation
- Isosurface extraction
- 3D data registration tools

---

## Requirements

### Required
- **MATLAB R2019b or later**
- **Statistics and Machine Learning Toolbox**
- **Image Processing Toolbox**

### Optional (for enhanced functionality)
- Curve Fitting Toolbox - for advanced peak fitting
- Optimization Toolbox - for parameter optimization
- Parallel Computing Toolbox - for accelerated batch processing

### Check Your Installation
```matlab
checkDependencies();
```

---

## Installation

### Option 1: Clone from GitHub (Recommended)
```bash
git clone https://github.com/peterfelfer/Atom-Probe-Toolbox.git
```

Then in MATLAB:
```matlab
cd('/path/to/Atom-Probe-Toolbox')
setupToolbox();
```

### Option 2: Download ZIP
1. Download the repository from [GitHub](https://github.com/peterfelfer/Atom-Probe-Toolbox)
2. Extract to your desired location
3. Run `setupToolbox()` in MATLAB

### Option 3: MATLAB Toolbox Installer
Download the `.mltbx` file from the Releases page and double-click to install.

### Make Installation Permanent
```matlab
setupToolbox('permanent', true);
```

---

## Quick Start

### Initialize the Toolbox
```matlab
setupToolbox();
```

### Load and Visualize Data
```matlab
% Load a .pos file
pos = posLoad('mydata.pos');

% Create a mass spectrum
massSpecPlot(pos.mc);

% 3D visualization
scatterPlotPosData(pos);
```

### Define Ions and Ranges
```matlab
% Create ion definitions
ions = ionAdd('Fe');
ions = ionAdd('Cr', ions);
ions = ionAdd('C', ions);

% Define mass-to-charge ranges
ranges = rangeAdd([55.8, 56.1], 'Fe', ions);
ranges = rangeAdd([51.9, 52.2], 'Cr', ions, ranges);
```

### Calculate Composition
```matlab
% Simple concentration
conc = posCalculateConcentrationSimple(pos, ranges);

% With uncertainty
unc = concentrationUncertainty(conc.counts);
```

### Cluster Analysis
```matlab
% Find clusters using DBSCAN
[clusterIdx, info] = clusterDBSCAN(pos, 0.5, 10);
fprintf('Found %d clusters\n', info.nClusters);
```

### Check Data Quality
```matlab
metrics = dataQualityMetrics(pos);
disp(metrics.summary);
```

---

## Documentation

### Getting Started Guide
Open the interactive getting started guide:
```matlab
open('doc/GettingStarted.mlx')
```

### Video Tutorials
Watch our [YouTube tutorial series](https://www.youtube.com/watch?v=Eid6iTGl48M&list=PLLr-VNeczShDzrm-wdIyz9ORjFL5KDeYT).

### Workflow Examples
The toolbox includes 30+ workflow files (`.mlx`) demonstrating common analysis tasks:
```matlab
% Open a workflow
open('Workflow_FirstSteps.mlx')
open('Workflow_ClusterDetermination.mlx')
open('Workflow_Proxigram.mlx')
```

---

## Function Reference

### Data I/O
| Function | Description |
|----------|-------------|
| `posLoad` | Load .pos/.epos files |
| `posExport` | Export position data |
| `rangesExtractFromFile` | Import RRNG range files |
| `hdf5FileCreateFromMetaDataList` | Create HDF5 database |
| `posTableAddToHDF5` | Save data to HDF5 |
| `posTableFromHDF5` | Load data from HDF5 |

### Ion & Range Management
| Function | Description |
|----------|-------------|
| `ionAdd` | Add ion to ion list |
| `ionConvertName` | Convert ion name formats |
| `ionsCreateIsotopeList` | Generate isotope patterns |
| `rangeAdd` | Add mass-to-charge range |
| `rangesFromPos` | Auto-generate ranges |

### Visualization
| Function | Description |
|----------|-------------|
| `massSpecPlot` | Plot mass spectrum |
| `scatterPlotPosData` | 3D atom visualization |
| `colorSchemeCreate` | Create color scheme |
| `colorSchemeIonAdd` | Add ion colors |

### ROI Operations
| Function | Description |
|----------|-------------|
| `roiCreateBox` | Create box ROI |
| `roiCreateCylinder` | Create cylindrical ROI |
| `roiCreateSphere` | Create spherical ROI |
| `roiCreatePlane` | Create planar ROI |
| `roiManipulate.mlapp` | Interactive ROI GUI |

### Analysis
| Function | Description |
|----------|-------------|
| `posCalculateConcentrationSimple` | Calculate composition |
| `concentrationUncertainty` | Uncertainty quantification |
| `pointCreateProxigram` | Point-based proxigram |
| `lineCreateProxigram` | Line-based proxigram |
| `patchCreateProxigram` | Surface proxigram |
| `spatialStatistics` | RDF, nearest neighbor analysis |
| `dataQualityMetrics` | Data quality assessment |
| `massSpecAnalysis` | Automated peak analysis |

### Cluster Analysis
| Function | Description |
|----------|-------------|
| `clusterDBSCAN` | DBSCAN clustering (GPU/CPU) |
| `clusterDetermination` | Cluster identification |
| `voronoiVolumeAnalysis` | Voronoi-based analysis |
| `clusterSizeAnalyse` | Cluster statistics |

### Crystallography
| Function | Description |
|----------|-------------|
| `crystalOrientation` | Orientation object |
| `ipfColor` | IPF coloring |
| `ipfHistogram` | IPF distribution |
| `misorientation` | Misorientation calculation |
| `boundaryCharacter` | GB characterization |
| `stereoProj` | Stereographic projection |

### Utilities
| Function | Description |
|----------|-------------|
| `setupToolbox` | Initialize toolbox |
| `checkDependencies` | Check requirements |
| `APTConfig` | Configuration management |
| `batchProcess` | Batch processing |
| `ProgressTracker` | Progress indication |
| `validateInputs` | Input validation |

---

## Workflows

| Workflow | Description |
|----------|-------------|
| `Workflow_FirstSteps.mlx` | Introduction and basics |
| `Workflow_1DConcentrationProfile.mlx` | 1D composition profiles |
| `Workflow_2DConcentrationMap.mlx` | 2D composition mapping |
| `Workflow_3D_Visualisation_of_APT_data.mlx` | 3D visualization |
| `Workflow_ClusterDetermination.mlx` | Cluster analysis |
| `Workflow_Proxigram.mlx` | Proxigram analysis |
| `Workflow_interfacialExcess.mlx` | Interfacial excess |
| `Workflow_Isosurface.mlx` | Isosurface extraction |
| `Workflow_ModellingGB.mlx` | Grain boundary modeling |
| `Workflow_HDF5_IO.mlx` | HDF5 database operations |
| `Workflow_reconstruction.mlx` | 3D reconstruction |
| `Workflow_ReflectronCorrection.mlx` | Reflectron correction |

---

## Configuration

### Instrument Presets
```matlab
cfg = APTConfig.getInstance();
cfg.applyInstrumentPreset('LEAP5000XR');  % or 'LEAP4000XHR', 'LEAP5000XS', 'EIKOS'
```

### Material Presets
```matlab
cfg.applyMaterialPreset('Al');  % Sets atomic volume for aluminum
% Available: Fe, Al, Cu, Ni, Ti, W, Si, Mg
```

### Custom Configuration
```matlab
cfg.reconstruction.detectionEfficiency = 0.52;
cfg.analysis.clusterEpsilon = 0.6;
cfg.save();  % Persist settings
```

---

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Ways to Contribute
- Report bugs and issues
- Suggest new features
- Submit pull requests
- Improve documentation
- Share workflow examples

---

## Citation

If you use this toolbox in your research, please cite:

```bibtex
@software{AtomProbeToolbox,
  author = {Felfer, Peter and Heller, Martina and Ott, Benedict and Dalbauer, Valentin},
  title = {Atom Probe Toolbox},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/peterfelfer/Atom-Probe-Toolbox}
}
```

---

## Authors

**Prof. Dr. Peter Felfer** - Principal Investigator

**Martina Heller** - Developer

**Benedict Ott** - Developer

**Dr. Valentin Dalbauer** - Developer

[Chair of General Materials Properties](https://www.wtm.tf.fau.de/)
Friedrich-Alexander-Universität Erlangen-Nürnberg

---

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- The atom probe community for feedback and feature requests
- [CAMECA](https://www.cameca.com/) for instrument support
- All contributors and users of this toolbox

---

**(c) 2024 Prof. Peter Felfer Group @ FAU Erlangen-Nürnberg**
