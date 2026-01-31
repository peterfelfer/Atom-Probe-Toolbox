function col = colorSchemeSync(col1, col2, varargin)
% colorSchemeSync synchronizes colors of col1 from col2.
%
% col = colorSchemeSync(col1, col2)
% col = colorSchemeSync(col1, col2, 'expand', true)
% col = colorSchemeSync(col1, col2, true)
%
% INPUT
% col1: table with columns ion and color (base scheme)
% col2: table with columns ion and color (source scheme)
%
% name-value options:
%   'expand'  true/false (default: false)
%             when true, append ions from col2 not present in col1
%
% OUTPUT
% col: color scheme table containing all ions in col1 with colors from col2
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if ~istable(col1) || ~istable(col2)
    error('colorSchemeSync:invalidInput', 'Inputs must be tables.');
end
if ~all(ismember({'ion','color'}, col1.Properties.VariableNames)) || ...
   ~all(ismember({'ion','color'}, col2.Properties.VariableNames))
    error('colorSchemeSync:invalidColorScheme', ...
        'Both color schemes must have columns ion and color.');
end

expand = false;
if ~isempty(varargin)
    if isscalar(varargin) && (islogical(varargin{1}) || isnumeric(varargin{1}))
        expand = logical(varargin{1});
    else
        if mod(numel(varargin), 2) ~= 0
            error('colorSchemeSync:invalidOptions', 'Options must be name-value pairs.');
        end
        for k = 1:2:numel(varargin)
            key = lower(string(varargin{k}));
            val = varargin{k + 1};
            switch key
                case "expand"
                    expand = logical(val);
                otherwise
                    error('colorSchemeSync:invalidOption', 'Unknown option "%s".', key);
            end
        end
    end
end

ion1 = col1.ion;
ion2 = col2.ion;

ion1Str = string(ion1);
ion2Str = string(ion2);

colors = col1.color;
if size(colors, 2) ~= 3
    error('colorSchemeSync:invalidColorScheme', 'col1.color must be Nx3 RGB values.');
end
if size(col2.color, 2) ~= 3
    error('colorSchemeSync:invalidColorScheme', 'col2.color must be Nx3 RGB values.');
end

[isMatch, idx] = ismember(ion1Str, ion2Str);
colors(isMatch, :) = col2.color(idx(isMatch), :);

col = col1;
col.color = colors;

if expand
    addMask = ~ismember(ion2Str, ion1Str);
    if any(addMask)
        addIons = col2.ion(addMask);
        addColors = col2.color(addMask, :);
        if iscategorical(ion1)
            col.ion = addcats(col.ion, categories(categorical(addIons)));
        end
        addTable = table(addIons, addColors, 'VariableNames', {'ion','color'});
        col = [col; addTable];
    end
end
end
