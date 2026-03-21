function pos = posLoadPyccapt(fileName, options)
% POSLOADPYCCAPT Load a pyccapt raw HDF5 file into a standard pos table.
%
% pos = posLoadPyccapt()
% pos = posLoadPyccapt(fileName)
% pos = posLoadPyccapt(fileName, 'maxTof', 5000)
%
% Reads the /dld/* datasets written by pyccapt instrument control and
% returns a MATLAB table compatible with the Atom Probe Toolbox.
%
% INPUT
%   fileName - path to .h5 file (file dialog if omitted)
%
% OPTIONS
%   'maxTof'  - discard ions with tof > maxTof (ns, default Inf)
%   'minTof'  - discard ions with tof < minTof (ns, default 50)
%   'quiet'   - suppress output messages (default false)
%
% OUTPUT
%   pos - table with columns: ionIdx, tof, VDC, VP, detx, dety
%         Units:              1,      ns,  V,   V,  mm,   mm
%         Optional extra columns: laserIntensity (pJ), startCounter (1)
%
% NOTE: The table does not contain x, y, z, or mc columns. Use
%       tofToMassToCharge to compute mc from tof, then a reconstruction
%       function (e.g., posReconstruct3DGeiser) for x, y, z.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
    options.maxTof (1,1) double = Inf
    options.minTof (1,1) double = 50
    options.quiet (1,1) logical = false
end

% File dialog if no path given
if isempty(fileName)
    [file, path] = uigetfile({'*.h5;*.hdf5', 'HDF5 files (*.h5, *.hdf5)'}, ...
                              'Select pyccapt HDF5 file');
    if isequal(file, 0)
        pos = table();
        return;
    end
    fileName = fullfile(path, file);
end

% Verify /dld group exists
info = h5info(fileName);
groupNames = {info.Groups.Name};
if ~any(strcmp(groupNames, '/dld'))
    error('posLoadPyccapt:noDldGroup', ...
        'File does not contain a /dld group. Not a pyccapt raw HDF5 file.');
end

% --- Read datasets ---
tof = readDataset(fileName, '/dld/t');
VDC = readDataset(fileName, '/dld/high_voltage');
n   = numel(tof);

% Pulse voltage: try several legacy key names
VP = tryRead(fileName, {'/dld/pulse', '/dld/voltage_pulse', '/dld/pulse_voltage'});
if isempty(VP)
    VP = zeros(n, 1);
end

% Detector positions: pyccapt stores in cm, toolbox uses mm
detx = readDataset(fileName, '/dld/x') * 10;   % cm -> mm
dety = readDataset(fileName, '/dld/y') * 10;   % cm -> mm

% Optional columns
laserIntensity = tryRead(fileName, {'/dld/laser_intensity'});
startCounter   = tryRead(fileName, {'/dld/start_counter'});

% --- Build table ---
ionIdx = (1:n)';
pos = table(ionIdx, tof, VDC, VP, detx, dety);
pos.Properties.VariableUnits = {'1', 'ns', 'V', 'V', 'mm', 'mm'};

if ~isempty(laserIntensity)
    pos.laserIntensity = laserIntensity(:);
    pos.Properties.VariableUnits{end} = 'pJ';
end
if ~isempty(startCounter)
    pos.startCounter = startCounter(:);
    pos.Properties.VariableUnits{end} = '1';
end

% --- Filter invalid rows ---
bad = (pos.tof < options.minTof) | ...
      (pos.tof > options.maxTof) | ...
      (pos.VDC < 0) | ...
      (pos.tof == 0 & pos.detx == 0 & pos.dety == 0);

if any(bad)
    if ~options.quiet
        fprintf('Removed %d invalid rows (%.1f%%).\n', sum(bad), 100*sum(bad)/n);
    end
    pos(bad, :) = [];
    pos.ionIdx = (1:height(pos))';
end

if ~options.quiet
    fprintf('Loaded %d ions from pyccapt file: %s\n', height(pos), fileName);
end

end


% ===== Local helpers =====================================================

function data = readDataset(fileName, datasetPath)
    data = double(h5read(fileName, datasetPath));
    data = data(:);
end

function data = tryRead(fileName, paths)
    data = [];
    for k = 1:numel(paths)
        try
            data = double(h5read(fileName, paths{k}));
            data = data(:);
            return;
        catch
        end
    end
end
