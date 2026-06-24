---
name: apt-reconstruct
description: Reconstruct 3D atom positions (x, y, z) from detector hit coordinates and the voltage curve using the Gault et al. straight-flight-path algorithm. Use when the user wants to reconstruct, build, or recompute the 3D point cloud from raw detector data.
argument-hint: "[element]"
---

Reconstruct the 3D atom positions from detector hit coordinates and the voltage curve, following Gault et al., Ultramicroscopy 111 (2011) 448-457. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — atom probe data **with detector coordinates and the voltage curve**: columns `detx`, `dety` (detector hit position in nm) and `VDC` (standing voltage). These exist only in extended formats (`.epos` / `.apt` / pyccapt HDF5), not in plain `.pos`. If only a `.pos` is loaded, tell the user to reload an `.epos`/`.apt`/HDF5 file via `/apt-load` first.

This skill overwrites `pos.x`, `pos.y`, `pos.z` with the newly computed reconstruction.

## Steps

### 1. Set reconstruction parameters

Instrument geometry and detection efficiency come from the instrument the data was acquired on. If `cfg` is available you can read defaults from it (`cfg = APTConfig.getInstance(); cfg.reconstruction.detectionEfficiency`, `cfg.reconstruction.flightPath`), otherwise ask the user which instrument was used and use the toolbox presets (LEAP4000XHR detEff 0.37 / flight 110, LEAP5000XR 0.52 / 120, LEAP5000XS 0.80 / 120, EIKOS 0.50 / 100).

```matlab
reconPara.detEff      = 0.36;   % detector efficiency (0..1)
reconPara.flightLength = 110;   % flight path in mm
reconPara.evapF       = 32.49;  % evaporation field in V/nm
reconPara.kf          = 3.3;    % field factor
reconPara.ICF         = 1.4;    % image compression factor
reconData = table;
```

`evapF`, `kf` and `ICF` are material/calibration dependent; the values above are the workflow defaults. Ask the user if they have calibrated values, otherwise keep these.

### 2. Average atomic density

Either set it directly (atoms/nm^3) or derive it from the main element of the sample. `$ARGUMENTS` is the element symbol (default `Fe`):

```matlab
load isotopeTable_naturalAbundances.mat;
element = "Fe";   % use $ARGUMENTS if provided
reconPara.avgDens = mean(isotopeTable.atomDensity(isotopeTable.element == element));
clearvars element isotopeTable
```

If the user knows the density directly, skip this and set `reconPara.avgDens = 85;` (example value).

### 3. Detector polar coordinates and radius evolution

```matlab
[reconData.ang, reconData.rad] = cart2pol(pos.detx, pos.dety);
reconPara.areaDet = (max(reconData.rad))^2 * pi();      % effective detector area
reconData.radEvol = pos.VDC/(reconPara.kf * reconPara.evapF);   % specimen radius (nm)
```

### 4. Calculate x and y coordinates

```matlab
reconData.thetaP = atan(reconData.rad / reconPara.flightLength);
reconData.theta  = reconData.thetaP + asin((reconPara.ICF - 1) * sin(reconData.thetaP));
[reconData.zP, reconData.distance] = pol2cart(reconData.theta, reconData.radEvol);
[pos.x, pos.y] = pol2cart(reconData.ang, reconData.distance);   % nm
```

### 5. Calculate z coordinate

```matlab
reconData.zP   = reconData.radEvol - reconData.zP;            % z shift from cap top
reconPara.omega = 1 ./ reconPara.avgDens;
reconData.dz   = reconPara.omega * reconPara.flightLength^2 * reconPara.kf^2 * reconPara.evapF^2 ...
                 / (reconPara.detEff * reconPara.areaDet * reconPara.ICF^2) * pos.VDC.^-2;
reconData.cumZ = cumsum(double(reconData.dz));               % accumulative depth
pos.z = reconData.cumZ + reconData.zP;                       % nm
```

## Report

Tell the user the reconstruction is complete and report the new coordinate ranges (`min`/`max` of `pos.x`, `pos.y`, `pos.z`) and the parameters used (detEff, flightLength, evapF, kf, ICF, avgDens). The intermediate values live in the `reconData` table and `reconPara` struct if they want to inspect them.

Suggest follow-ups:
- Visualize the result: `/apt-visualize`
- Build the mass spectrum: `/apt-spectrum`
- Export the reconstructed data: `posExport(pos, 'reconstructed.pos')`
- Adjust `kf` / `ICF` and re-run if the aspect ratio or feature spacing looks wrong (these are the main tuning knobs).

## Toolbox functions used

`cart2pol`, `pol2cart`, `cumsum` (MATLAB built-ins); reconstruction is implemented inline per the Gault algorithm. Optional: `APTConfig.getInstance`, `posExport`, `posLoad`.
