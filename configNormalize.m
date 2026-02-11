function valueOut = configNormalize(valueIn)
% CONFIGNORMALIZE Convert MATLAB values into a text-serializable config subset.
%
% valueOut = configNormalize(valueIn)
%
% Supported output types are struct, cell, char, double, logical, and []
% with nested combinations. This function is used before JSON/YAML export.

valueOut = normalizeValue(valueIn, 'root');

end

function valueOut = normalizeValue(valueIn, path)
if nargin < 2
    path = 'root';
end

if isempty(valueIn)
    valueOut = [];
    return;
end

if istable(valueIn)
    error('configNormalize:unsupportedTable', ...
        'Unsupported table at %s. Convert tables to struct first.', path);
end

if isdatetime(valueIn) || isduration(valueIn) || isa(valueIn, 'calendarDuration')
    if isscalar(valueIn)
        valueOut = char(string(valueIn));
    else
        valueOut = cellstr(string(valueIn(:)));
    end
    return;
end

if iscategorical(valueIn)
    if isscalar(valueIn)
        valueOut = char(string(valueIn));
    else
        valueOut = cellstr(string(valueIn(:)));
    end
    return;
end

if isstring(valueIn)
    if isscalar(valueIn)
        valueOut = char(valueIn);
    else
        valueOut = cellstr(valueIn(:));
    end
    return;
end

if ischar(valueIn)
    valueOut = valueIn;
    return;
end

if isnumeric(valueIn) || islogical(valueIn)
    valueOut = valueIn;
    return;
end

if iscell(valueIn)
    valueOut = cell(size(valueIn));
    for i = 1:numel(valueIn)
        valueOut{i} = normalizeValue(valueIn{i}, sprintf('%s{%d}', path, i));
    end
    return;
end

if isstruct(valueIn)
    if isscalar(valueIn)
        valueOut = struct();
        fieldNames = fieldnames(valueIn);
        for i = 1:numel(fieldNames)
            fieldName = fieldNames{i};
            valueOut.(fieldName) = normalizeValue(valueIn.(fieldName), ...
                sprintf('%s.%s', path, fieldName));
        end
    else
        valueOut = cell(numel(valueIn), 1);
        for i = 1:numel(valueIn)
            valueOut{i} = normalizeValue(valueIn(i), sprintf('%s(%d)', path, i));
        end
    end
    return;
end

error('configNormalize:unsupportedType', ...
    'Unsupported value type ''%s'' at %s.', class(valueIn), path);
end
