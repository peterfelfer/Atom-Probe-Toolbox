---
name: apt-roi
description: Define a region of interest (box, cylinder, sphere, or plane) on a 3D atom probe reconstruction, crop the data inside it, and plot a 1D concentration profile along the ROI axis. Use when the user wants to select, crop, or analyze a local sub-volume of their APT data.
argument-hint: "[shape]"
---

Create a region of interest on the current 3D reconstruction, crop the enclosed data, and optionally extract a 1D concentration profile along the ROI axis. Use the MATLAB MCP server.

## Prerequisites

The MATLAB workspace must contain:
- `pos` — loaded atom probe data with ions allocated (has `ion` column, from `/apt-composition`)
- `ax` — handle to the current 3D axis of a scatter plot (from `/apt-visualize`)
- `colorScheme` — color scheme (from `/apt-setup`)
- `ionTable` — table of ions in the dataset (only needed for the profile exclude-list step)

If there is no 3D figure yet, tell the user to run `/apt-visualize` first.
If the user opened a previously saved MATLAB figure instead, tell them to make it the current figure and run this once before anything else:

```matlab
ax = gca;
```

## Steps

### 1. Create the ROI

Pick the shape from `$ARGUMENTS` (box, cylinder, sphere, or plane). If the user did not specify, ask. All position/dimension inputs are optional — the defaults below match the workflow. Locations are `[x y z]` in nm.

Box (length, width, height in x/y/z):
```matlab
dimensions = [10 10 10];   % length, width, height (x,y,z) in nm
location   = [0 0 40];     % origin coordinates
roi_box = roiCreateBox(dimensions, location, ax);
```

Cylinder:
```matlab
radius      = 5;           % nm
height      = 10;          % nm
location    = [0 0 40];
numSegments = 32;          % number of facets approximating the cylinder
roi_cylinder = roiCreateCylinder(radius, height, location, numSegments, ax);
```

Sphere:
```matlab
radius       = 15;         % nm
subDivisions = 15;         % number of subdivisions
location     = [0 0 40];
roi_sphere = roiCreateSphere(radius, subDivisions, location, ax);
```

Plane (no enclosed volume — for orientation/cutting reference only):
```matlab
dimensions = [70 70];      % length (x) and width (y)
spacing    = [35 35];      % subdivision spacing in x and y
location   = [0 0 40];     % plane center
roi_plane = roiCreatePlane(dimensions, spacing, location, ax);
```

### 2. Manipulate the ROI (interactive, optional)

Tell the user this is a manual GUI step — do not automate it. Launch the app, then have the user move/rotate/resize the ROI by hand:

```matlab
roiManipulate;
```

In the app the user types the ROI variable name (e.g. `roi_box`), clicks **get object**, then edits size and location fields or drags the ROI. There is also an *execute after move* field for re-running a command on each change.

### 3. Crop the data inside the ROI

Only 3D ROIs with convex bounding surfaces (box, cylinder, sphere) can be used here — not the plane. Select the ROI created above:

```matlab
roi = roi_box;   % choose the ROI you created (roi_box / roi_cylinder / roi_sphere)
pos_roi = pos(posInConvexHull(pos, roi.Vertices), :);   % keep only ions inside the ROI
```

`pos_roi` is a standard pos table you can pass to any other skill (composition, cluster, etc.).

### 4. (Optional) 1D concentration profile along the ROI axis

A 1D profile only makes sense for a box or cylinder. Stop here if the user only wanted cropped data.

Transform the cropped data into the ROI's own coordinate system (needed if the ROI was rotated relative to the dataset axes):
```matlab
pos_roi_trans = coordSystemTransform(roi, pos_roi);
```

Create grid vectors along the profile axis (here z). Use mode `'distance'` for a fixed slice thickness, or `'count'` for a fixed number of atoms per slice:
```matlab
dist = pos_roi_trans.z;    % profile direction (z is the ROI long axis)
bin  = 1;                  % slice thickness in nm (mode 'distance')
mode = 'distance';
[binCenters, binEdges] = binVectorsFromDistance(dist, bin, mode);
```
For the count-based alternative use `bin = 5000; mode = 'count';` (atoms per slice).

Voxelise (slice) the cropped data:
```matlab
vox = posNdBin(pos_roi_trans, dist, binEdges);
```

Build the concentration kernel. `detEff` may be a decimal or a percentage; `excludeListConc` and `volumeName` are optional:
```matlab
detEff          = 37;                          % detection efficiency (% or decimal); ask which instrument
excludeListConc = {'unranged', 'unknown'};
volumeName      = 'AP06';
concentrationKernel = @(p) posCalculateConcentrationSimple(p, detEff, excludeListConc, volumeName);
```
Ask the user which instrument was used so `detEff` is correct (LEAP4000XHR 37, LEAP5000XR 52, LEAP5000XS 80, EIKOS 50).

Choose which species to plot (exclude everything else):
```matlab
listAllIons     = (cellstr(ionTable.ionName))';
desIons         = {'C', 'N', 'Fe N', 'N C', 'C2'};   % species the user wants to see
excludeListPlot = setdiff(listAllIons, desIons);
```

Apply the kernel and plot:
```matlab
conc1D = binApplyConcentrationKernel(vox, binCenters, concentrationKernel, {'nm'});
[p, ax, f] = concentrationProfilePlot(conc1D(conc1D.format=='concentration', :), excludeListPlot, colorScheme);
```

## Report

Tell the user:
- Which ROI shape/dimensions were created and how many atoms fell inside (`height(pos_roi)`).
- That `pos_roi` is now in the workspace for further analysis.
- If a profile was made, summarise the trend and offer to export: `writetable(conc1D, 'profile.csv')`.

Suggested follow-ups: `/apt-composition` (composition of just the ROI via `pos_roi`), `/apt-cluster` (clustering inside the ROI), or `/apt-isosurface` for a 3D concentration view.
