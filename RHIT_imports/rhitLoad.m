function [hits, histograms, metadata] = rhitLoad(fileName)
% RHITLOAD Load raw detector data from a Cameca RHIT file.
%
% [hits, histograms, metadata] = rhitLoad(fileName)
%
% The RHIT format is a CERN ROOT-based proprietary format used by Cameca
% LEAP atom probes. It contains raw (pre-reconstruction) detector hit data.
%
% This function requires Python 3 with the 'uproot', 'h5py', and 'numpy'
% packages installed. Install them via:
%   pip3 install uproot h5py numpy
%
% INPUTS:
%   fileName - Path to the .RHIT file. If empty, opens a file dialog.
%
% OUTPUTS:
%   hits        - Table with per-event detector data in epos-compatible
%                 format. Standard columns come first:
%                   ionIdx     - sequential event index (1-based)
%                   detx, dety - detector position (mm, converted from LSB)
%                   mc         - mass-to-charge (Da, voltage-corrected;
%                                dominant peak ~0.3%, others ~2-5% without
%                                bowl correction)
%                   tof        - time-of-flight (ns)
%                   VDC        - specimen DC voltage (V)
%                   Vref       - reference voltage (V, = VDC * kf)
%                 Additional columns:
%                   detxRaw, detyRaw - raw detector position (LSB)
%                   pulse      - pulse voltage/energy (arb)
%                   freq       - pulse frequency (Hz)
%                   chi2       - fit quality
%                   hreg       - hit region
%                   tElapsed   - elapsed time (s)
%                   erate      - evaporation rate
%                   Pres       - chamber pressure (mbar)
%                   Temp       - specimen temperature (K)
%                   AmbTemp    - ambient temperature (C)
%                   VMcpGain, VAnodeAccel - MCP/anode voltages (V)
%                   Noise      - noise level
%                   and more instrument parameters
%
%   histograms  - Struct with stored histograms:
%                   .massSpectrum   - mass-to-charge spectrum
%                   .tofRaw         - raw time-of-flight spectrum
%                   .voltageHistory - voltage vs event number
%                   .erateHistory   - evaporation rate history
%                   .detectorXY     - 2D detector hit map
%                 Each histogram has .values, .edges, and .title fields.
%
%   metadata    - Struct with file metadata and summary info
%
% NOTE: The RHIT file contains RAW detector data, not reconstructed 3D
% positions. To obtain a standard pos table (x,y,z,mc), you need to:
%   1. Convert detector x,y from LSB to mm using detector calibration
%   2. Apply mass calibration: mc = f(tof, V, flightPath)
%   3. Perform 3D reconstruction (e.g., posReconstruct3DGeiser)
%
% EXAMPLE:
%   [hits, hist] = rhitLoad('R56_03622.RHIT');
%   plot(hist.massSpectrum.edges(1:end-1), hist.massSpectrum.values);
%   xlabel('Mass-to-charge (Da)'); ylabel('Counts');
%
% See also: posLoad, posReconstruct3DGeiser
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
end

if isempty(fileName)
    [file, path] = uigetfile({'*.RHIT;*.rhit', 'RHIT files (*.RHIT)'}, ...
        'Select RHIT file');
    if isequal(file, 0)
        hits = table();
        histograms = struct();
        metadata = struct();
        return;
    end
    fileName = fullfile(path, file);
end

% Locate the Python extraction script
scriptDir = fileparts(mfilename('fullpath'));
pyScript = fullfile(scriptDir, 'rhitExtract.py');

if ~isfile(pyScript)
    error('rhitLoad:missingScript', ...
        'Cannot find rhitExtract.py in %s', scriptDir);
end

% Create temporary HDF5 file for data transfer
h5File = [tempname, '.h5'];
cleanupObj = onCleanup(@() deleteIfExists(h5File));

% Run Python extraction
cmd = sprintf('python3 "%s" "%s" "%s"', pyScript, fileName, h5File);
[status, output] = system(cmd);

if status ~= 0
    error('rhitLoad:pythonError', ...
        'Python extraction failed:\n%s\nEnsure uproot, h5py, numpy are installed:\n  pip3 install uproot h5py numpy', ...
        output);
end

if ~isfile(h5File)
    error('rhitLoad:noOutput', 'Python script did not produce output file.');
end

fprintf('%s\n', strtrim(output));

% Read the HDF5 data into MATLAB
rawHits = readHits(h5File);
histograms = readHistograms(h5File);
metadata = readMetadata(h5File, fileName);

% Convert to pos-table-compatible format
hits = convertToPosFormat(rawHits, metadata);

fprintf('Loaded %d raw detector hits from %s\n', height(hits), fileName);

end


function hits = readHits(h5File)
% Read the per-event hit data from the HDF5 file

info = h5info(h5File, '/hits');

% Get all scalar datasets
hits = table();
for i = 1:numel(info.Datasets)
    dsName = info.Datasets(i).Name;
    data = h5read(h5File, ['/hits/' dsName]);
    hits.(dsName) = data;
end

% Read struct sub-groups (Thresh, UpperThresh, Walk)
for i = 1:numel(info.Groups)
    grpName = info.Groups(i).Name;
    [~, shortName] = fileparts(grpName);
    for j = 1:numel(info.Groups(i).Datasets)
        dsName = info.Groups(i).Datasets(j).Name;
        fieldName = [shortName '_' dsName];
        data = h5read(h5File, [grpName '/' dsName]);
        hits.(fieldName) = data;
    end
end
end


function histograms = readHistograms(h5File)
% Read stored histograms

histograms = struct();
info = h5info(h5File, '/histograms');

for i = 1:numel(info.Groups)
    grpPath = info.Groups(i).Name;
    [~, hName] = fileparts(grpPath);

    h = struct();

    % Read datasets
    for j = 1:numel(info.Groups(i).Datasets)
        dsName = info.Groups(i).Datasets(j).Name;
        h.(dsName) = h5read(h5File, [grpPath '/' dsName]);
    end

    % Read attributes
    for j = 1:numel(info.Groups(i).Attributes)
        attrName = info.Groups(i).Attributes(j).Name;
        h.(attrName) = info.Groups(i).Attributes(j).Value;
    end

    histograms.(hName) = h;
end
end


function metadata = readMetadata(h5File, fileName)
% Collect metadata

metadata = struct();
metadata.fileName = fileName;
metadata.format = 'RHIT (Cameca ROOT)';

% Read file-level attributes
info = h5info(h5File, '/');
for i = 1:numel(info.Attributes)
    attrName = info.Attributes(i).Name;
    attrName = matlab.lang.makeValidName(attrName);
    metadata.(attrName) = info.Attributes(i).Value;
end

% Hit count
hitsInfo = h5info(h5File, '/hits');
for i = 1:numel(hitsInfo.Attributes)
    if strcmp(hitsInfo.Attributes(i).Name, 'num_entries')
        metadata.numHits = double(hitsInfo.Attributes(i).Value);
    end
end

% pElf summary data if available
try
    elfInfo = h5info(h5File, '/pElf');
    elf = struct();
    for i = 1:numel(elfInfo.Datasets)
        dsName = elfInfo.Datasets(i).Name;
        elf.(dsName) = h5read(h5File, ['/pElf/' dsName]);
    end
    metadata.pElf = elf;
catch
    % pElf not available
end

% Instrument parameters from CRunHeader
try
    paramsInfo = h5info(h5File, '/instrumentParams');
    params = struct();
    for i = 1:numel(paramsInfo.Attributes)
        attrName = paramsInfo.Attributes(i).Name;
        attrName = matlab.lang.makeValidName(attrName);
        params.(attrName) = paramsInfo.Attributes(i).Value;
    end
    for i = 1:numel(paramsInfo.Datasets)
        dsName = paramsInfo.Datasets(i).Name;
        dsName = matlab.lang.makeValidName(dsName);
        params.(dsName) = h5read(h5File, ['/instrumentParams/' paramsInfo.Datasets(i).Name]);
    end
    metadata.instrumentParams = params;
catch
    % instrumentParams not available
end
end


function hits = convertToPosFormat(rawHits, metadata)
% Convert raw RHIT table to pos-table-compatible column names and units.
%
% Mapping:
%   RHIT 'x','y' (LSB)  →  'detxRaw','detyRaw' (LSB) + 'detx','dety' (mm)
%   RHIT 'v'             →  'VDC' (specimen DC voltage, V)
%   RHIT 'z' (seq idx)   →  'ionIdx' (1-based)
%   RHIT 'tof'           →  'tof' (ns)
%   RHIT 'Vref'          →  'Vref' (reference voltage, V)
%   All other columns are kept with their original names.

hits = rawHits;

% --- Detector coordinates: LSB → mm ---
if isfield(metadata, 'instrumentParams') && ...
        isfield(metadata.instrumentParams, 'lsb_to_mm')
    lsb2mm = metadata.instrumentParams.lsb_to_mm;
else
    lsb2mm = 18.5 / 750;  % default for LEAP 4000 Hamamatsu MCP
end

hits.detxRaw = hits.x;
hits.detyRaw = hits.y;
hits.detx = double(hits.x) * lsb2mm;
hits.dety = double(hits.y) * lsb2mm;
hits.x = [];
hits.y = [];

% --- ionIdx from sequence index (1-based) ---
hits.ionIdx = double(hits.z) + 1;
hits.z = [];

% --- Rename voltage column ---
hits.VDC = hits.v;
hits.v = [];

% --- Compute mass-to-charge ratio with voltage correction ---
% mc = (2*e / (u * L^2)) * VDC * (tof - t0)^2 * 1e-18 / sqrt(kf)
%
% Parameters from CRunHeader:
%   L  = flight path (mm)
%   t0 = time-of-flight offset (ns)
%   kf = voltage correction factor (Vref/VDC ratio, typically 1.03)
%
% The sqrt(kf) correction accounts for the reflectron geometry where kf
% affects both the effective voltage and the flight path. Validated against
% paired epos data: dominant peak accuracy ~0.3%.
%
% No bowl correction (position-dependent flight path) is applied, so
% secondary peaks may have ~2-5% offset. The fully calibrated Cameca
% spectrum is available in histograms.massSpectrum.
if isfield(metadata, 'instrumentParams')
    L  = metadata.instrumentParams.flight_path_mm / 1000; % m
    t0 = metadata.instrumentParams.t0_ns;                 % ns
    % kf = Vref / VDC (voltage correction factor)
    if ismember('Vref', hits.Properties.VariableNames) && hits.VDC(1) > 0
        kf = double(hits.Vref(1)) / double(hits.VDC(1));
    else
        kf = 1.03;
    end
else
    L  = 0.382;
    t0 = 45.0;
    kf = 1.03;
end

e_const = 1.602176634e-19;   % C
u_const = 1.66053906660e-27; % kg
C_calib = 2 * e_const / (u_const * L^2) * 1e-18;

hits.mc = C_calib .* double(hits.VDC) .* (double(hits.tof) - t0).^2 ...
    ./ sqrt(kf);

% --- Reorder: put standard epos-like columns first ---
stdCols = {'ionIdx', 'detx', 'dety', 'mc', 'tof', 'VDC'};
if ismember('Vref', hits.Properties.VariableNames)
    stdCols{end+1} = 'Vref';
end
otherCols = setdiff(hits.Properties.VariableNames, stdCols, 'stable');
hits = hits(:, [stdCols, otherCols]);

% --- Units ---
units = repmat({''}, 1, width(hits));
unitMap = struct( ...
    'ionIdx', '1', ...
    'mc', 'Da', ...
    'detx', 'mm', ...
    'dety', 'mm', ...
    'detxRaw', 'LSB', ...
    'detyRaw', 'LSB', ...
    'tof', 'ns', ...
    'VDC', 'V', ...
    'Vref', 'V', ...
    'pulse', 'arb', ...
    'freq', 'Hz', ...
    'tElapsed', 's', ...
    'erate', '1', ...
    'TargetErate', '1', ...
    'Pres', 'mbar', ...
    'Temp', 'K', ...
    'AmbTemp', 'C', ...
    'VMcpGain', 'V', ...
    'VAnodeAccel', 'V');

names = hits.Properties.VariableNames;
for i = 1:numel(names)
    if isfield(unitMap, names{i})
        units{i} = unitMap.(names{i});
    end
end
hits.Properties.VariableUnits = units;

end


function deleteIfExists(filePath)
if isfile(filePath)
    delete(filePath);
end
end
