# Changelog

All notable changes to the Atom Probe Toolbox will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [1.1.0] - 2024-XX-XX

### Added

#### New Analysis Functions
- `concentrationUncertainty` - Proper uncertainty quantification using binomial, Poisson, or Clopper-Pearson methods
- `spatialStatistics` - Comprehensive spatial statistics including RDF g(r), nearest neighbor distributions, and Ripley's K function
- `dataQualityMetrics` - Data quality assessment with resolution estimation, artifact detection, and density analysis
- `massSpecAnalysis` - Automated mass spectrum analysis with peak detection, background fitting, and resolution estimation

#### Cluster Analysis
- `clusterDBSCAN` - DBSCAN clustering with automatic GPU acceleration and CPU fallback
  - Supports solute-specific clustering
  - Automatic memory management for large datasets
  - Returns comprehensive cluster statistics

#### Crystallography
- `misorientation` - Calculate misorientation angles between crystal orientations
  - Supports cubic, hexagonal, tetragonal, orthorhombic, monoclinic, triclinic crystal systems
  - Automatic CSL boundary identification (Sigma values)
  - Multiple input formats (Euler angles, quaternions, rotation matrices)

#### Infrastructure
- `setupToolbox` - One-command toolbox initialization
  - Adds all required paths
  - Optional dependency checking
  - Permanent path saving option
- `checkDependencies` - Comprehensive dependency verification
  - Checks MATLAB version and toolboxes
  - Verifies internal modules
  - Reports GPU and parallel computing availability
- `APTConfig` - Centralized configuration management
  - Instrument presets (LEAP 4000X HR, 5000XR, 5000XS, EIKOS)
  - Material presets (Fe, Al, Cu, Ni, Ti, W, Si, Mg)
  - Persistent user settings
- `batchProcess` - Batch processing framework
  - Sequential and parallel processing
  - Error handling and logging
  - Progress tracking
- `runTests` - Test runner framework for automated testing

#### Utilities
- `ProgressTracker` - Progress indication for long-running operations
  - Command-line and GUI modes
  - Automatic ETA calculation
  - Cancellation support
- `validateInputs` - Unified input validation
  - Supports posTable, ranges, ions, patches, quaternions, etc.
  - Informative error messages

### Changed
- Updated README.md with comprehensive documentation
- Reorganized `helptoc.xml` with proper function categories
- Standardized function documentation format

### Documentation
- New comprehensive README with installation guide, quick start, and API reference
- API_Reference.md - Quick reference guide for all functions
- CONTRIBUTING.md - Contribution guidelines
- TROUBLESHOOTING.md - Common issues and solutions
- Proper helptoc.xml structure for MATLAB help browser

---

## [1.0.0] - 2024-09-28

### Initial Release

#### Core Functionality
- Mass spectrum analysis and ion ranging
- 3D point cloud visualization
- Region of interest (ROI) selection and manipulation
- Concentration calculation
- Proxigram analysis
- Interfacial excess calculations
- Isosurface analysis
- HDF5 database support with metadata

#### Data Import/Export
- POS and EPOS file support
- APT file support (CAMECA)
- HDF5 database creation and querying
- OBJ and PLY mesh export

#### Analysis Methods
- Simple and kernelized concentration
- 1D and 2D concentration profiles
- Point, line, and surface-based proxigrams
- Voronoi volume analysis
- Cluster determination

#### Crystallography
- Crystal orientation handling
- IPF coloring and histograms
- Grain boundary characterization
- Stereographic projections

#### Reconstruction
- Geiser reconstruction algorithm
- Reflectron distortion correction

#### GUI Applications
- `ionFind.mlapp` - Ion identification
- `roiManipulate.mlapp` - ROI manipulation
- `cropping.mlapp` - Data cropping
- `reconExplorer.mlapp` - Reconstruction explorer
- `isosurfaceWidget.mlapp` - Isosurface visualization
- `volumeIsosurface.mlapp` - Volume isosurface extraction

#### Workflows
- 30+ workflow examples (.mlx files)
- Getting started guide
- Video tutorial series

---

## Future Plans

### Planned for v1.2.0
- [ ] Machine learning integration for peak identification
- [ ] Automated phase identification
- [ ] Enhanced correlative microscopy tools
- [ ] EBSD data integration

### Under Consideration
- [ ] Web-based visualization export
- [ ] Python interoperability
- [ ] Cloud processing support

---

## Version History

| Version | Date | MATLAB | Notes |
|---------|------|--------|-------|
| 1.1.0 | TBD | R2019b+ | New analysis functions, improved infrastructure |
| 1.0.0 | 2024-09-28 | R2019b+ | Initial public release |

---

## Contributors

- Prof. Dr. Peter Felfer
- Martina Heller
- Benedict Ott
- Dr. Valentin Dalbauer

---

**(c) Prof. Peter Felfer Group @ FAU Erlangen-Nürnberg**
