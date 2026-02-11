function valueOut = configDenormalize(valueIn)
% CONFIGDENORMALIZE Convert imported config values into practical MATLAB types.
%
% valueOut = configDenormalize(valueIn)
%
% This function primarily converts YAML parser cell-array sequences back into
% numeric/logical vectors/matrices or struct arrays where possible.

valueOut = denormalizeValue(valueIn);

end

function valueOut = denormalizeValue(valueIn)
if isempty(valueIn)
    valueOut = valueIn;
    return;
end

if isnumeric(valueIn) || islogical(valueIn) || ischar(valueIn)
    valueOut = valueIn;
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

if isstruct(valueIn)
    if isscalar(valueIn)
        valueOut = struct();
        fieldNames = fieldnames(valueIn);
        for i = 1:numel(fieldNames)
            fieldName = fieldNames{i};
            valueOut.(fieldName) = denormalizeValue(valueIn.(fieldName));
        end
    else
        valueOut = valueIn;
        for i = 1:numel(valueIn)
            valueOut(i) = denormalizeValue(valueIn(i));
        end
    end
    return;
end

if iscell(valueIn)
    valueCell = cell(size(valueIn));
    for i = 1:numel(valueIn)
        valueCell{i} = denormalizeValue(valueIn{i});
    end

    % Collapse to numeric/logical vector when possible
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), valueCell(:)))
        valueOut = cell2mat(valueCell(:))';
        return;
    end
    if all(cellfun(@(x) islogical(x) && isscalar(x), valueCell(:)))
        valueOut = logical(cell2mat(valueCell(:))');
        return;
    end

    % Collapse nested numeric rows to matrix when possible
    if all(cellfun(@(x) isnumeric(x) && isvector(x), valueCell(:)))
        lengths = cellfun(@numel, valueCell(:));
        if ~isempty(lengths) && all(lengths == lengths(1))
            try
                valueOut = vertcat(valueCell{:});
                return;
            catch
            end
        end
    end

    % Collapse to struct array when all entries are scalar structs
    if all(cellfun(@(x) isstruct(x) && isscalar(x), valueCell(:)))
        try
            valueOut = vertcat(valueCell{:});
            return;
        catch
        end
    end

    valueOut = valueCell;
    return;
end

valueOut = valueIn;
end
