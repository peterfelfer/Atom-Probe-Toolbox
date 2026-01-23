# Troubleshooting Guide

Solutions to common issues when using the Atom Probe Toolbox.

---

## Table of Contents

- [Installation Issues](#installation-issues)
- [Data Import Problems](#data-import-problems)
- [Memory Issues](#memory-issues)
- [Performance Problems](#performance-problems)
- [Analysis Errors](#analysis-errors)
- [Visualization Issues](#visualization-issues)
- [HDF5 Database Issues](#hdf5-database-issues)
- [FAQ](#faq)

---

## Installation Issues

### "Function not found" errors

**Problem:** MATLAB cannot find toolbox functions.

**Solution:**
```matlab
% Run the setup script
setupToolbox();

% Or manually add paths
addpath(genpath('/path/to/Atom-Probe-Toolbox'));
```

### Missing toolbox dependencies

**Problem:** Error about missing Statistics or Image Processing Toolbox.

**Solution:**
```matlab
% Check what's installed
checkDependencies();

% Install missing toolboxes via MATLAB Add-On Explorer
% Or contact your MATLAB administrator
```

### MATLAB version incompatibility

**Problem:** Syntax errors or missing functions.

**Solution:**
The toolbox requires MATLAB R2019b or later. Key features requiring R2019b:
- `arguments` blocks for input validation
- Some string functions

If using an older version, you may need to modify functions that use `arguments` blocks.

---

## Data Import Problems

### Cannot read .pos file

**Problem:** Error when loading .pos or .epos files.

**Solutions:**
```matlab
% Check file exists
if ~isfile('mydata.pos')
    error('File not found');
end

% Try specifying the format explicitly
pos = posLoad('mydata.pos', 'format', 'pos');

% Check file permissions
fileattrib('mydata.pos')
```

### Corrupted or truncated files

**Problem:** File loads partially or with errors.

**Solution:**
```matlab
% Check file size - typical .pos files are ~16 bytes per atom
fileInfo = dir('mydata.pos');
estimatedAtoms = fileInfo.bytes / 16;
fprintf('Estimated atoms: %d\n', estimatedAtoms);

% Try loading with error tolerance
pos = posLoad('mydata.pos', 'skipErrors', true);
```

### Wrong coordinate system

**Problem:** Atoms appear in wrong positions or orientation.

**Solution:**
```matlab
% Check data range
fprintf('X: [%.2f, %.2f]\n', min(pos.x), max(pos.x));
fprintf('Y: [%.2f, %.2f]\n', min(pos.y), max(pos.y));
fprintf('Z: [%.2f, %.2f]\n', min(pos.z), max(pos.z));

% Data might need transformation
% Swap axes if needed
pos_corrected = pos;
pos_corrected.x = pos.y;
pos_corrected.y = pos.x;
```

---

## Memory Issues

### "Out of memory" errors

**Problem:** MATLAB runs out of memory with large datasets.

**Solutions:**

1. **Check available memory:**
```matlab
memory  % Windows only
% On Mac/Linux, check Activity Monitor
```

2. **Reduce dataset size:**
```matlab
% Random sampling
sampleIdx = randperm(height(pos), 1000000);
posSample = pos(sampleIdx, :);

% Spatial cropping
mask = pos.x > 0 & pos.x < 50;
posCropped = pos(mask, :);
```

3. **Use chunked processing:**
```matlab
% Process in chunks
cfg = APTConfig.getInstance();
cfg.performance.chunkSize = 500000;
```

4. **Clear unused variables:**
```matlab
clear largeVariable
pack  % Defragment memory (use sparingly)
```

### Slow HDF5 operations

**Problem:** Reading/writing HDF5 is very slow.

**Solution:**
```matlab
% Use chunked reading for large datasets
pos = posTableFromHDF5(h5file, dataset, 'chunkSize', 100000);

% Enable compression
cfg = APTConfig.getInstance();
cfg.io.hdf5Compression = 'deflate';
cfg.io.hdf5CompressionLevel = 6;
```

---

## Performance Problems

### Slow cluster analysis

**Problem:** DBSCAN or Voronoi analysis takes too long.

**Solutions:**

1. **Enable parallel processing:**
```matlab
% Start parallel pool
parpool('local');

% Enable parallel in config
cfg = APTConfig.getInstance();
cfg.performance.useParallel = true;
```

2. **Reduce dataset size:**
```matlab
% Analyze only solute atoms
soluteIdx = strcmp(pos.ion, 'Cu');
[idx, info] = clusterDBSCAN(pos, 0.5, 10, 'soluteIdx', soluteIdx);
```

3. **Adjust parameters:**
```matlab
% Larger epsilon = faster but less precise
[idx, info] = clusterDBSCAN(pos, 0.8, 10);  % Instead of 0.5
```

### Slow 3D visualization

**Problem:** Scatter plots are slow or unresponsive.

**Solutions:**

1. **Reduce point count:**
```matlab
% Display subset
sampleIdx = 1:10:height(pos);  % Every 10th atom
scatterPlotPosData(pos(sampleIdx, :));
```

2. **Use simpler markers:**
```matlab
scatter3(pos.x, pos.y, pos.z, 1, '.');  % Small dots
```

3. **Disable unnecessary features:**
```matlab
ax = gca;
ax.SortMethod = 'childorder';  % Faster rendering
```

---

## Analysis Errors

### "Index exceeds array bounds"

**Problem:** Array indexing error during analysis.

**Common causes and solutions:**

1. **Empty data after filtering:**
```matlab
% Check before using
mask = pos.mc > 55 & pos.mc < 56;
if ~any(mask)
    warning('No atoms in range');
    return;
end
posFiltered = pos(mask, :);
```

2. **Mismatched array sizes:**
```matlab
% Ensure arrays match
assert(height(pos) == length(ionLabels), 'Size mismatch');
```

### NaN or Inf in results

**Problem:** Analysis produces NaN or Inf values.

**Solutions:**

1. **Check for division by zero:**
```matlab
concentration = counts ./ totalCounts;
concentration(totalCounts == 0) = 0;  % Handle zero division
```

2. **Check input data:**
```matlab
% Find problematic values
hasNaN = any(isnan(pos{:, {'x','y','z'}}), 2);
hasInf = any(isinf(pos{:, {'x','y','z'}}), 2);
fprintf('NaN rows: %d, Inf rows: %d\n', sum(hasNaN), sum(hasInf));

% Remove them
pos = pos(~hasNaN & ~hasInf, :);
```

### Unexpected concentration values

**Problem:** Concentrations don't sum to 1 or seem wrong.

**Solutions:**

1. **Check range definitions:**
```matlab
% Ensure ranges don't overlap
for i = 1:height(ranges)
    for j = i+1:height(ranges)
        if ranges.mcbegin(i) < ranges.mcend(j) && ranges.mcend(i) > ranges.mcbegin(j)
            warning('Ranges %d and %d overlap', i, j);
        end
    end
end
```

2. **Check for unranged atoms:**
```matlab
% How many atoms are not assigned to any ion?
ranged = false(height(pos), 1);
for i = 1:height(ranges)
    ranged = ranged | (pos.mc >= ranges.mcbegin(i) & pos.mc <= ranges.mcend(i));
end
fprintf('Unranged atoms: %.1f%%\n', 100 * sum(~ranged) / height(pos));
```

---

## Visualization Issues

### Blank or black figure

**Problem:** Plot appears empty or all black.

**Solutions:**

1. **Check data range:**
```matlab
% Ensure data is not all zeros or identical
if range(pos.x) == 0
    warning('All x-coordinates are identical');
end
```

2. **Reset axes:**
```matlab
axis auto
view(3)
```

3. **Check alpha values:**
```matlab
% Fully transparent objects won't show
scatter3(pos.x, pos.y, pos.z, 1, 'filled', 'MarkerFaceAlpha', 0.5);
```

### Colors not showing correctly

**Problem:** Ion colors don't match expected.

**Solution:**
```matlab
% Check color scheme
colorScheme = colorSchemeCreate();
colorScheme = colorSchemeIonAdd('Fe', [1 0 0], colorScheme);  % Red
colorScheme = colorSchemeIonAdd('Cu', [0 0 1], colorScheme);  % Blue

% Apply explicitly
scatterPlotPosData(pos, 'colorScheme', colorScheme);
```

### Figure window freezes

**Problem:** MATLAB becomes unresponsive during plotting.

**Solution:**
```matlab
% Use drawnow to force updates
for i = 1:nFrames
    updatePlot(i);
    drawnow limitrate;  % Limits update frequency
end

% Or disable figure updates temporarily
set(gcf, 'Visible', 'off');
% ... do plotting ...
set(gcf, 'Visible', 'on');
drawnow;
```

---

## HDF5 Database Issues

### Cannot create HDF5 file

**Problem:** Error when creating new HDF5 database.

**Solutions:**

1. **Check write permissions:**
```matlab
[status, msg] = fileattrib(folder);
if ~status || ~msg.UserWrite
    error('No write permission');
end
```

2. **Close existing file handles:**
```matlab
% If file is locked
fclose('all');
% Or restart MATLAB
```

### Cannot read dataset

**Problem:** Dataset exists but cannot be read.

**Solution:**
```matlab
% Check dataset info
info = h5info(h5file, dataset);
disp(info);

% Check datatype
fprintf('Datatype: %s\n', info.Datatype.Class);

% Try reading with specific type
data = h5read(h5file, dataset, 'TextEncoding', 'UTF-8');
```

### Metadata not saving

**Problem:** Attributes not appearing in HDF5 file.

**Solution:**
```matlab
% Ensure attribute values are correct type
h5writeatt(h5file, dataset, 'description', 'My data');  % String
h5writeatt(h5file, dataset, 'count', int32(100));       % Integer
h5writeatt(h5file, dataset, 'voltage', 5.5);            % Double
```

---

## FAQ

### Q: Which MATLAB version do I need?

**A:** MATLAB R2019b or later is required. R2021a or later is recommended for best performance.

### Q: Can I use the toolbox without the Statistics Toolbox?

**A:** Some functions will work, but many analysis features require the Statistics and Machine Learning Toolbox. Run `checkDependencies()` to see what's available.

### Q: How do I cite this toolbox?

**A:** See the CITATION section in README.md for the recommended citation format.

### Q: My data is in a different format. Can I still use the toolbox?

**A:** Yes, if you can load your data as a table with columns `x`, `y`, `z`, and `mc` (mass-to-charge), you can use most functions. See `posTableFromHDF5` for the expected format.

### Q: How do I report a bug?

**A:** Open an issue on GitHub with:
- MATLAB version
- Operating system
- Minimal code to reproduce
- Full error message

### Q: Can I contribute new functions?

**A:** Yes! See CONTRIBUTING.md for guidelines.

---

## Still Having Issues?

1. Check the [GitHub Issues](https://github.com/peterfelfer/Atom-Probe-Toolbox/issues)
2. Watch the [Video Tutorials](https://www.youtube.com/playlist?list=PLLr-VNeczShDzrm-wdIyz9ORjFL5KDeYT)
3. Open a new issue with detailed information

---

**(c) Prof. Peter Felfer Group @ FAU Erlangen-Nürnberg**
