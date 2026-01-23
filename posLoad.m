function pos = posLoad(fileName, options)
% POSLOAD Canonical loader for APT datasets (.pos, .epos, .apt, .h5)
%
% pos = posLoad(fileName)
% pos = posLoad(fileName, 'format', 'pos')
% pos = posLoad(fileName, 'skipErrors', true)
%
% OPTIONS:
%   'format'     - Explicit format override: 'pos', 'epos', 'apt', 'h5'
%   'skipErrors' - Return empty table on failure (default: false)
%   'quiet'      - Suppress informational output (default: false)
%   'dataset'    - Reserved for future HDF5 dataset override
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
    options.format (1,:) char = ''
    options.skipErrors (1,1) logical = false
    options.quiet (1,1) logical = false
    options.dataset (1,:) char = ''
end

if isempty(fileName)
    [file, path] = uigetfile({
        '*.pos;*.epos;*.apt;*.h5;*.hdf5', 'APT files (*.pos, *.epos, *.apt, *.h5, *.hdf5)';
        '*.*', 'All files (*.*)'
    }, 'Select APT data file');
    if isequal(file, 0)
        pos = table();
        return;
    end
    fileName = fullfile(path, file);
end

if isstring(fileName)
    fileName = char(fileName);
end

if istable(fileName)
    pos = fileName;
    return;
end

format = normalizeFormat(options.format, fileName);

try
    switch format
        case {'pos', 'epos'}
            pos = posToTable(fileName);
        case 'apt'
            pos = aptToTable(fileName);
        case {'h5', 'hdf5'}
            if ~isempty(options.dataset)
                error('posLoad:datasetUnsupported', ...
                    'HDF5 dataset override not supported.');
            end
            pos = posTableFromHDF5(fileName);
        otherwise
            error('posLoad:unknownFormat', 'Unknown file format: %s', format);
    end
catch ME
    if options.skipErrors
        warning('posLoad:loadFailed', 'Failed to load %s: %s', fileName, ME.message);
        pos = table();
        return;
    end
    rethrow(ME);
end

pos = normalizePosTable(pos);

if ~options.quiet
    if istable(pos)
        fprintf('Loaded %d ions from %s\n', height(pos), fileName);
    end
end

end

function format = normalizeFormat(formatOpt, fileName)
    format = lower(strtrim(formatOpt));
    format = regexprep(format, '^\.+', '');
    if isempty(format) || strcmp(format, 'auto')
        [~, ~, ext] = fileparts(fileName);
        format = lower(strrep(ext, '.', ''));
    end
end

function pos = normalizePosTable(pos)
    if isempty(pos) || ~istable(pos)
        return;
    end

    if ~ismember('ionIdx', pos.Properties.VariableNames)
        pos = addvars(pos, (1:height(pos))', 'Before', 1, 'NewVariableNames', 'ionIdx');
    end

    units = pos.Properties.VariableUnits;
    if isempty(units) || numel(units) ~= width(pos)
        units = repmat({''}, 1, width(pos));
    end

    unitMap = struct(...
        'ionIdx', '1', ...
        'x', 'nm', ...
        'y', 'nm', ...
        'z', 'nm', ...
        'mc', 'Da', ...
        'tof', 'ns', ...
        'VDC', 'V', ...
        'VP', 'V', ...
        'detx', 'mm', ...
        'dety', 'mm', ...
        'deltaP', '1', ...
        'multi', '1', ...
        'atomNum', '1');

    names = pos.Properties.VariableNames;
    for i = 1:numel(names)
        name = names{i};
        if isfield(unitMap, name)
            units{i} = unitMap.(name);
        elseif isempty(units{i})
            units{i} = '';
        end
    end

    pos.Properties.VariableUnits = units;
end
