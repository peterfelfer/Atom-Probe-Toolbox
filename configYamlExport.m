function configYamlExport(config, fileName)
% CONFIGYAMLEXPORT Export a configuration struct to YAML (dependency-free subset).
%
% configYamlExport(config, fileName)

arguments
    config (1,1) struct
    fileName (1,:) char
end

configValidate(config);
cfg = configNormalize(config);

fid = fopen(fileName, 'w');
if fid < 0
    error('configYamlExport:fileOpenFailed', ...
        'Could not open file for writing: %s', fileName);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '# Atom Probe Toolbox configuration (YAML subset)\n');
writeNode(fid, cfg, 0);

end

function writeNode(fid, value, indent)
if isstruct(value)
    writeStruct(fid, value, indent);
elseif iscell(value)
    writeCell(fid, value, indent);
elseif isnumeric(value)
    writeNumeric(fid, value, indent);
elseif islogical(value)
    writeLogical(fid, value, indent);
elseif ischar(value)
    fprintf(fid, '%s%s\n', indentStr(indent), quoteString(value));
elseif isempty(value)
    fprintf(fid, '%snull\n', indentStr(indent));
else
    error('configYamlExport:unsupportedType', ...
        'Unsupported value type ''%s'' for YAML export.', class(value));
end
end

function writeStruct(fid, s, indent)
if isempty(s)
    fprintf(fid, '%s{}\n', indentStr(indent));
    return;
end

if ~isscalar(s)
    items = num2cell(s(:));
    writeCell(fid, items, indent);
    return;
end

fieldNames = fieldnames(s);
if isempty(fieldNames)
    fprintf(fid, '%s{}\n', indentStr(indent));
    return;
end

for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    value = s.(fieldName);
    key = formatKey(fieldName);
    if isInlineValue(value)
        fprintf(fid, '%s%s: %s\n', indentStr(indent), key, encodeInline(value));
    else
        fprintf(fid, '%s%s:\n', indentStr(indent), key);
        writeNode(fid, value, indent + 2);
    end
end
end

function writeCell(fid, c, indent)
if isempty(c)
    fprintf(fid, '%s[]\n', indentStr(indent));
    return;
end

for i = 1:numel(c)
    value = c{i};
    if isInlineValue(value)
        fprintf(fid, '%s- %s\n', indentStr(indent), encodeInline(value));
    else
        fprintf(fid, '%s-\n', indentStr(indent));
        writeNode(fid, value, indent + 2);
    end
end
end

function writeNumeric(fid, x, indent)
if isempty(x)
    fprintf(fid, '%s[]\n', indentStr(indent));
    return;
end

if isscalar(x) || isvector(x)
    fprintf(fid, '%s%s\n', indentStr(indent), encodeInline(x));
    return;
end

for i = 1:size(x, 1)
    fprintf(fid, '%s- %s\n', indentStr(indent), encodeInline(x(i, :)));
end
end

function writeLogical(fid, x, indent)
if isempty(x)
    fprintf(fid, '%s[]\n', indentStr(indent));
    return;
end

if isscalar(x) || isvector(x)
    fprintf(fid, '%s%s\n', indentStr(indent), encodeInline(x));
    return;
end

for i = 1:size(x, 1)
    fprintf(fid, '%s- %s\n', indentStr(indent), encodeInline(x(i, :)));
end
end

function tf = isInlineValue(value)
if isempty(value)
    tf = true;
    return;
end

if ischar(value)
    tf = true;
    return;
end

if isnumeric(value) || islogical(value)
    tf = isscalar(value) || isvector(value);
    return;
end

if iscell(value)
    tf = all(cellfun(@isInlineValue, value(:)));
    return;
end

if isstruct(value)
    tf = false;
    return;
end

tf = false;
end

function out = encodeInline(value)
if isempty(value)
    out = 'null';
    return;
end

if ischar(value)
    out = quoteString(value);
    return;
end

if islogical(value)
    if isscalar(value)
        out = ternary(value, 'true', 'false');
    else
        items = arrayfun(@(v) ternary(v, 'true', 'false'), value(:)', 'UniformOutput', false);
        out = ['[', strjoin(items, ', '), ']'];
    end
    return;
end

if isnumeric(value)
    if isscalar(value)
        out = numericToken(value);
    else
        items = arrayfun(@(v) numericToken(v), value(:)', 'UniformOutput', false);
        out = ['[', strjoin(items, ', '), ']'];
    end
    return;
end

if iscell(value)
    parts = cell(size(value));
    for i = 1:numel(value)
        parts{i} = encodeInline(value{i});
    end
    out = ['[', strjoin(parts, ', '), ']'];
    return;
end

error('configYamlExport:unsupportedInlineType', ...
    'Unsupported inline value type ''%s''.', class(value));
end

function token = numericToken(v)
if isnan(v)
    token = '.nan';
elseif isinf(v) && v > 0
    token = '.inf';
elseif isinf(v) && v < 0
    token = '-.inf';
elseif abs(v - round(v)) < eps(max(1, abs(v)))
    token = sprintf('%.0f', round(v));
else
    token = sprintf('%.16g', v);
end
end

function out = quoteString(str)
str = char(str);
str = strrep(str, '\\', '\\\\');
str = strrep(str, '"', '\\"');
out = ['"', str, '"'];
end

function out = indentStr(indent)
out = repmat(' ', 1, max(0, indent));
end

function key = formatKey(fieldName)
if all(isstrprop(fieldName, 'alphanum') | fieldName == '_')
    key = fieldName;
else
    key = quoteString(fieldName);
end
end

function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
