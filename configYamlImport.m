function config = configYamlImport(fileName)
% CONFIGYAMLIMPORT Import a configuration struct from YAML (dependency-free subset).
%
% config = configYamlImport(fileName)

arguments
    fileName (1,:) char
end

if ~isfile(fileName)
    error('configYamlImport:fileNotFound', ...
        'YAML file not found: %s', fileName);
end

text = fileread(fileName);
lines = splitlines(string(text));
lines = stripCarriageReturns(lines);

[idxStart, hasContent] = firstContentLine(lines, 1);
if ~hasContent
    config = struct();
    return;
end

[value, idxOut] = parseBlock(lines, idxStart, 0); %#ok<ASGLU>
idxAfter = skipIgnorable(lines, idxOut);
if idxAfter <= numel(lines)
    error('configYamlImport:trailingContent', ...
        'Unexpected trailing YAML content near line %d.', idxAfter);
end

config = configDenormalize(value);
if ~isstruct(config) || ~isscalar(config)
    error('configYamlImport:invalidTopLevel', ...
        'Top-level YAML node must decode to a scalar struct.');
end

end

function lines = stripCarriageReturns(lines)
for i = 1:numel(lines)
    lines(i) = replace(lines(i), sprintf('\r'), '');
end
end

function idx = skipIgnorable(lines, idx)
while idx <= numel(lines)
    line = char(lines(idx));
    line = stripInlineComment(line);
    if isempty(strtrim(line))
        idx = idx + 1;
    else
        break;
    end
end
end

function [idx, hasContent] = firstContentLine(lines, idx)
idx = skipIgnorable(lines, idx);
hasContent = idx <= numel(lines);
end

function [value, idxOut] = parseBlock(lines, idx, indent)
idx = skipIgnorable(lines, idx);
if idx > numel(lines)
    value = [];
    idxOut = idx;
    return;
end

line = char(lines(idx));
lineNoComment = stripInlineComment(line);
currentIndent = countLeadingSpaces(lineNoComment);
if currentIndent < indent
    value = [];
    idxOut = idx;
    return;
end
if currentIndent > indent
    error('configYamlImport:indentation', ...
        'Unexpected indentation at line %d.', idx);
end

trimmed = strtrim(lineNoComment);
    if startsWith(trimmed, '-')
        [value, idxOut] = parseSequence(lines, idx, indent);
    else
        [value, idxOut] = parseMapping(lines, idx, indent);
end
end

function [value, idxOut] = parseMapping(lines, idx, indent)
value = struct();
idxCurrent = idx;

while true
    idxCurrent = skipIgnorable(lines, idxCurrent);
    if idxCurrent > numel(lines)
        break;
    end

    line = char(lines(idxCurrent));
    lineNoComment = stripInlineComment(line);
    currentIndent = countLeadingSpaces(lineNoComment);
    if currentIndent < indent
        break;
    end
    if currentIndent > indent
        error('configYamlImport:indentation', ...
            'Unexpected indentation at line %d.', idxCurrent);
    end

    trimmed = strtrim(lineNoComment);
    if startsWith(trimmed, '-')
        break;
    end

    colonPos = find(trimmed == ':', 1, 'first');
    if isempty(colonPos)
        error('configYamlImport:missingColon', ...
            'Invalid mapping entry at line %d.', idxCurrent);
    end

    rawKey = strtrim(trimmed(1:colonPos-1));
    rawValue = strtrim(trimmed(colonPos+1:end));
    key = parseYamlKey(rawKey, idxCurrent);

    if isempty(rawValue)
        [child, nextIdx] = parseBlock(lines, idxCurrent + 1, indent + 2);
        value.(key) = child;
        idxCurrent = nextIdx;
    else
        value.(key) = parseScalarOrFlow(rawValue, idxCurrent);
        idxCurrent = idxCurrent + 1;
    end
end

idxOut = idxCurrent;
end

function [value, idxOut] = parseSequence(lines, idx, indent)
items = {};
idxCurrent = idx;

while true
    idxCurrent = skipIgnorable(lines, idxCurrent);
    if idxCurrent > numel(lines)
        break;
    end

    line = char(lines(idxCurrent));
    lineNoComment = stripInlineComment(line);
    currentIndent = countLeadingSpaces(lineNoComment);
    if currentIndent < indent
        break;
    end
    if currentIndent > indent
        error('configYamlImport:indentation', ...
            'Unexpected indentation at line %d.', idxCurrent);
    end

    trimmed = strtrim(lineNoComment);
    if ~startsWith(trimmed, '-')
        break;
    end

    if strcmp(trimmed, '-')
        itemToken = '';
    elseif startsWith(trimmed, '- ')
        itemToken = strtrim(trimmed(3:end));
    else
        error('configYamlImport:flowSyntax', ...
            'Invalid sequence item at line %d.', idxCurrent);
    end
    if isempty(itemToken)
        [child, nextIdx] = parseBlock(lines, idxCurrent + 1, indent + 2);
        items{end+1, 1} = child; %#ok<AGROW>
        idxCurrent = nextIdx;
    else
        items{end+1, 1} = parseScalarOrFlow(itemToken, idxCurrent); %#ok<AGROW>
        idxCurrent = idxCurrent + 1;
    end
end

value = items;
idxOut = idxCurrent;
end

function value = parseScalarOrFlow(token, lineNumber)
token = strtrim(token);
if isempty(token)
    value = [];
    return;
end

if startsWith(token, '[') && endsWith(token, ']')
    value = parseFlowSequence(token, lineNumber);
    return;
end

value = parseScalar(token);
end

function value = parseFlowSequence(token, lineNumber)
inner = strtrim(token(2:end-1));
if isempty(inner)
    value = [];
    return;
end

parts = splitFlowItems(inner, lineNumber);
values = cell(numel(parts), 1);
for i = 1:numel(parts)
    values{i} = parseScalarOrFlow(parts{i}, lineNumber);
end

if all(cellfun(@(x) isnumeric(x) && isscalar(x), values))
    value = cell2mat(values)';
elseif all(cellfun(@(x) islogical(x) && isscalar(x), values))
    value = logical(cell2mat(values)');
else
    value = values;
end
end

function parts = splitFlowItems(text, lineNumber)
parts = {};
startIdx = 1;
depth = 0;
inDouble = false;
inSingle = false;
i = 1;

while i <= strlength(text)
    ch = extractBetween(text, i, i);
    c = char(ch);

    if inDouble
        if c == '\\'
            i = i + 1;
        elseif c == '"'
            inDouble = false;
        end
    elseif inSingle
        if c == ''''
            inSingle = false;
        end
    else
        if c == '"'
            inDouble = true;
        elseif c == ''''
            inSingle = true;
        elseif c == '['
            depth = depth + 1;
        elseif c == ']'
            depth = depth - 1;
            if depth < 0
                error('configYamlImport:flowSyntax', ...
                    'Invalid flow sequence syntax at line %d.', lineNumber);
            end
        elseif c == ',' && depth == 0
            parts{end+1, 1} = strtrim(char(extractBetween(text, startIdx, i-1))); %#ok<AGROW>
            startIdx = i + 1;
        end
    end
    i = i + 1;
end

parts{end+1, 1} = strtrim(char(extractBetween(text, startIdx, strlength(text))));
end

function value = parseScalar(token)
if strcmp(token, 'null') || strcmp(token, '~')
    value = [];
    return;
end
if strcmp(token, 'true')
    value = true;
    return;
end
if strcmp(token, 'false')
    value = false;
    return;
end
if strcmp(token, '.nan') || strcmpi(token, 'nan')
    value = NaN;
    return;
end
if strcmp(token, '.inf') || strcmpi(token, 'inf')
    value = Inf;
    return;
end
if strcmp(token, '-.inf') || strcmpi(token, '-inf')
    value = -Inf;
    return;
end

if startsWith(token, '"') && endsWith(token, '"')
    raw = token(2:end-1);
    raw = strrep(raw, '\\"', '"');
    raw = strrep(raw, '\\\\', '\\');
    value = raw;
    return;
end

if startsWith(token, '''') && endsWith(token, '''')
    raw = token(2:end-1);
    raw = strrep(raw, '''''', '''');
    value = raw;
    return;
end

numericValue = str2double(token);
if ~isnan(numericValue)
    value = numericValue;
else
    value = token;
end
end

function key = parseYamlKey(rawKey, lineNumber)
if startsWith(rawKey, '"') && endsWith(rawKey, '"')
    key = rawKey(2:end-1);
elseif startsWith(rawKey, '''') && endsWith(rawKey, '''')
    key = rawKey(2:end-1);
else
    key = rawKey;
end

if isempty(key)
    error('configYamlImport:emptyKey', 'Empty key at line %d.', lineNumber);
end

if ~isvarname(key)
    error('configYamlImport:invalidKey', ...
        'YAML key ''%s'' at line %d is not a valid MATLAB struct field name.', ...
        key, lineNumber);
end
end

function n = countLeadingSpaces(line)
if isempty(line)
    n = 0;
    return;
end
n = 0;
while n < numel(line) && line(n+1) == ' '
    n = n + 1;
end
end

function out = stripInlineComment(line)
out = line;
inDouble = false;
inSingle = false;
i = 1;
while i <= numel(line)
    c = line(i);
    if inDouble
        if c == '\' && i < numel(line)
            i = i + 1;
        elseif c == '"'
            inDouble = false;
        end
    elseif inSingle
        if c == ''''
            inSingle = false;
        end
    else
        if c == '"'
            inDouble = true;
        elseif c == ''''
            inSingle = true;
        elseif c == '#'
            out = line(1:i-1);
            return;
        end
    end
    i = i + 1;
end
end
