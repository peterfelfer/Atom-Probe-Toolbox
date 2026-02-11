function config = configImport(fileName)
% CONFIGIMPORT Import toolbox configuration structs from MAT/JSON/YAML.
%
% config = configImport(fileName)

arguments
    fileName (1,:) char
end

if ~isfile(fileName)
    error('configImport:fileNotFound', ...
        'Configuration file not found: %s', fileName);
end

[~, ~, ext] = fileparts(fileName);
ext = lower(ext);

switch ext
    case '.mat'
        loaded = load(fileName);
        if isfield(loaded, 'config')
            config = loaded.config;
        else
            names = fieldnames(loaded);
            if isempty(names)
                error('configImport:emptyMat', ...
                    'MAT file ''%s'' does not contain any variables.', fileName);
            end
            config = loaded.(names{1});
        end

    case '.json'
        txt = fileread(fileName);
        cfg = jsondecode(txt);
        config = configDenormalize(cfg);

    case {'.yaml', '.yml'}
        config = configYamlImport(fileName);

    otherwise
        error('configImport:unsupportedExtension', ...
            'Unsupported extension ''%s''. Use .mat, .json, .yaml, or .yml.', ext);
end

if ~isstruct(config) || ~isscalar(config)
    error('configImport:invalidType', ...
        'Imported configuration must be a scalar struct.');
end

end
