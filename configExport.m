function configExport(config, fileName, options)
% CONFIGEXPORT Export toolbox configuration structs to MAT/JSON/YAML.
%
% configExport(config, fileName)
% configExport(config, fileName, 'prettyJson', true)

arguments
    config (1,1) struct
    fileName (1,:) char
    options.prettyJson (1,1) logical = true
end

configValidate(config);
[~, ~, ext] = fileparts(fileName);
ext = lower(ext);

switch ext
    case '.mat'
        save(fileName, 'config');

    case '.json'
        cfg = configNormalize(config);
        jsonText = jsonencode(cfg, 'PrettyPrint', options.prettyJson);
        writeTextFile(fileName, jsonText);

    case {'.yaml', '.yml'}
        configYamlExport(config, fileName);

    otherwise
        error('configExport:unsupportedExtension', ...
            'Unsupported extension ''%s''. Use .mat, .json, .yaml, or .yml.', ext);
end

end

function writeTextFile(fileName, textIn)
fid = fopen(fileName, 'w');
if fid < 0
    error('configExport:fileOpenFailed', ...
        'Could not open file for writing: %s', fileName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fwrite(fid, textIn, 'char');
end
