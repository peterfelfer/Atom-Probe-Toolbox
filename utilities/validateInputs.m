function validateInputs(varargin)
% VALIDATEINPUTS Unified input validation for Atom Probe Toolbox functions
%
% validateInputs(name1, value1, type1, name2, value2, type2, ...)
%
% Validates input arguments with informative error messages.
%
% INPUTS:
%   Arguments come in triplets: (name, value, type)
%   - name: string, the parameter name for error messages
%   - value: the value to validate
%   - type: string or cell array specifying validation rules
%
% SUPPORTED TYPES:
%   'posTable'      - Valid position table with required columns
%   'ranges'        - Valid ranges table
%   'ions'          - Valid ions table
%   'numeric'       - Numeric array (any size)
%   'scalar'        - Single numeric value
%   'positive'      - Positive numeric value(s)
%   'nonnegative'   - Non-negative numeric value(s)
%   'integer'       - Integer value(s)
%   'vector'        - 1D numeric array
%   'matrix'        - 2D numeric array
%   'string'        - Character array or string
%   'logical'       - Logical value
%   'table'         - MATLAB table
%   'struct'        - Structure
%   'cell'          - Cell array
%   'function'      - Function handle
%   'file'          - Existing file path
%   'directory'     - Existing directory path
%   'patch'         - Valid patch structure (vertices, faces)
%   'colorScheme'   - Valid color scheme table
%   'quaternion'    - Valid quaternion (4-element unit vector)
%   'rotationMatrix'- Valid 3x3 rotation matrix
%   {'option', {'a','b','c'}} - Value must be one of the listed options
%   {'size', [m,n]}          - Array must have specified size
%   {'minLength', n}         - Minimum length for arrays
%   {'columns', {'col1','col2'}} - Table must have specified columns
%
% EXAMPLES:
%   validateInputs('pos', pos, 'posTable', 'radius', r, 'positive');
%   validateInputs('method', method, {'option', {'linear','nearest'}});
%   validateInputs('data', data, {'size', [3, NaN]}); % 3 rows, any columns
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if mod(nargin, 3) ~= 0
    error('validateInputs:invalidArguments', ...
        'Arguments must come in triplets: (name, value, type)');
end

nArgs = nargin / 3;

for i = 1:nArgs
    idx = (i-1)*3 + 1;
    name = varargin{idx};
    value = varargin{idx + 1};
    typeSpec = varargin{idx + 2};

    validateSingleInput(name, value, typeSpec);
end

end

function validateSingleInput(name, value, typeSpec)
% Validate a single input argument

% Handle cell array type specifications (compound types)
if iscell(typeSpec)
    validateCompoundType(name, value, typeSpec);
    return;
end

% Simple type validation
switch lower(typeSpec)
    case 'postable'
        validatePosTable(name, value);

    case 'ranges'
        validateRanges(name, value);

    case 'ions'
        validateIons(name, value);

    case 'numeric'
        if ~isnumeric(value)
            error('validateInputs:invalidType', ...
                '%s must be numeric.', name);
        end

    case 'scalar'
        if ~isnumeric(value) || ~isscalar(value)
            error('validateInputs:invalidType', ...
                '%s must be a numeric scalar.', name);
        end

    case 'positive'
        if ~isnumeric(value) || any(value(:) <= 0)
            error('validateInputs:invalidValue', ...
                '%s must be positive.', name);
        end

    case 'nonnegative'
        if ~isnumeric(value) || any(value(:) < 0)
            error('validateInputs:invalidValue', ...
                '%s must be non-negative.', name);
        end

    case 'integer'
        if ~isnumeric(value) || any(mod(value(:), 1) ~= 0)
            error('validateInputs:invalidValue', ...
                '%s must be integer-valued.', name);
        end

    case 'vector'
        if ~isnumeric(value) || (~isvector(value) && ~isempty(value))
            error('validateInputs:invalidType', ...
                '%s must be a numeric vector.', name);
        end

    case 'matrix'
        if ~isnumeric(value) || ndims(value) > 2
            error('validateInputs:invalidType', ...
                '%s must be a 2D numeric matrix.', name);
        end

    case 'string'
        if ~ischar(value) && ~isstring(value)
            error('validateInputs:invalidType', ...
                '%s must be a character array or string.', name);
        end

    case 'logical'
        if ~islogical(value) && ~(isnumeric(value) && all(value(:) == 0 | value(:) == 1))
            error('validateInputs:invalidType', ...
                '%s must be logical.', name);
        end

    case 'table'
        if ~istable(value)
            error('validateInputs:invalidType', ...
                '%s must be a table.', name);
        end

    case 'struct'
        if ~isstruct(value)
            error('validateInputs:invalidType', ...
                '%s must be a structure.', name);
        end

    case 'cell'
        if ~iscell(value)
            error('validateInputs:invalidType', ...
                '%s must be a cell array.', name);
        end

    case 'function'
        if ~isa(value, 'function_handle')
            error('validateInputs:invalidType', ...
                '%s must be a function handle.', name);
        end

    case 'file'
        if ~ischar(value) && ~isstring(value)
            error('validateInputs:invalidType', ...
                '%s must be a file path string.', name);
        end
        if ~isfile(value)
            error('validateInputs:fileNotFound', ...
                'File not found: %s', value);
        end

    case 'directory'
        if ~ischar(value) && ~isstring(value)
            error('validateInputs:invalidType', ...
                '%s must be a directory path string.', name);
        end
        if ~isfolder(value)
            error('validateInputs:directoryNotFound', ...
                'Directory not found: %s', value);
        end

    case 'patch'
        validatePatch(name, value);

    case 'colorscheme'
        validateColorScheme(name, value);

    case 'quaternion'
        if ~isnumeric(value) || numel(value) ~= 4
            error('validateInputs:invalidQuaternion', ...
                '%s must be a 4-element quaternion.', name);
        end
        normVal = norm(value);
        if abs(normVal - 1) > 1e-6
            error('validateInputs:invalidQuaternion', ...
                '%s must be a unit quaternion (norm = 1, got %.6f).', name, normVal);
        end

    case 'rotationmatrix'
        if ~isnumeric(value) || ~isequal(size(value), [3, 3])
            error('validateInputs:invalidRotationMatrix', ...
                '%s must be a 3x3 matrix.', name);
        end
        % Check orthogonality: R * R' should be identity
        RRt = value * value';
        if norm(RRt - eye(3), 'fro') > 1e-6
            error('validateInputs:invalidRotationMatrix', ...
                '%s must be an orthogonal matrix (R * R'' = I).', name);
        end
        % Check determinant = 1 (proper rotation, not reflection)
        if abs(det(value) - 1) > 1e-6
            error('validateInputs:invalidRotationMatrix', ...
                '%s must have determinant = 1 (got %.6f).', name, det(value));
        end

    otherwise
        error('validateInputs:unknownType', ...
            'Unknown validation type: %s', typeSpec);
end

end

function validateCompoundType(name, value, typeSpec)
% Handle compound type specifications like {'option', {'a','b'}}

if isempty(typeSpec)
    return;
end

specType = typeSpec{1};

switch lower(specType)
    case 'option'
        options = typeSpec{2};
        if ~any(strcmpi(value, options))
            optStr = strjoin(options, ', ');
            error('validateInputs:invalidOption', ...
                '%s must be one of: %s (got: %s)', name, optStr, char(value));
        end

    case 'size'
        expectedSize = typeSpec{2};
        actualSize = size(value);

        % NaN in expectedSize means "any" for that dimension
        for d = 1:length(expectedSize)
            if ~isnan(expectedSize(d)) && ...
               (d > length(actualSize) || actualSize(d) ~= expectedSize(d))
                error('validateInputs:invalidSize', ...
                    '%s has wrong size. Expected [%s], got [%s].', ...
                    name, num2str(expectedSize), num2str(actualSize));
            end
        end

    case 'minlength'
        minLen = typeSpec{2};
        if length(value) < minLen
            error('validateInputs:invalidLength', ...
                '%s must have at least %d elements (got %d).', ...
                name, minLen, length(value));
        end

    case 'columns'
        requiredCols = typeSpec{2};
        if ~istable(value)
            error('validateInputs:invalidType', ...
                '%s must be a table.', name);
        end
        missingCols = setdiff(requiredCols, value.Properties.VariableNames);
        if ~isempty(missingCols)
            error('validateInputs:missingColumns', ...
                '%s is missing required columns: %s', name, strjoin(missingCols, ', '));
        end

    case 'range'
        bounds = typeSpec{2};
        if ~isnumeric(value) || any(value(:) < bounds(1)) || any(value(:) > bounds(2))
            error('validateInputs:outOfRange', ...
                '%s must be in range [%g, %g].', name, bounds(1), bounds(2));
        end

    otherwise
        error('validateInputs:unknownCompoundType', ...
            'Unknown compound validation type: %s', specType);
end

end

function validatePosTable(name, value)
% Validate a position table

if ~istable(value)
    error('validateInputs:invalidPosTable', ...
        '%s must be a table.', name);
end

requiredCols = {'x', 'y', 'z'};
missingCols = setdiff(requiredCols, lower(value.Properties.VariableNames));
if ~isempty(missingCols)
    error('validateInputs:invalidPosTable', ...
        '%s must contain columns: x, y, z. Missing: %s', name, strjoin(missingCols, ', '));
end

end

function validateRanges(name, value)
% Validate a ranges table

if ~istable(value)
    error('validateInputs:invalidRanges', ...
        '%s must be a table.', name);
end

requiredCols = {'mcbegin', 'mcend'};
actualCols = lower(value.Properties.VariableNames);
missingCols = setdiff(lower(requiredCols), actualCols);
if ~isempty(missingCols)
    error('validateInputs:invalidRanges', ...
        '%s must contain columns: mcbegin, mcend. Missing: %s', name, strjoin(missingCols, ', '));
end

end

function validateIons(name, value)
% Validate an ions table

if ~istable(value)
    error('validateInputs:invalidIons', ...
        '%s must be a table.', name);
end

requiredCols = {'ion', 'chargeState'};
actualCols = value.Properties.VariableNames;
missingCols = setdiff(requiredCols, actualCols);
if ~isempty(missingCols)
    error('validateInputs:invalidIons', ...
        '%s must contain columns: ion, chargeState. Missing: %s', name, strjoin(missingCols, ', '));
end

end

function validatePatch(name, value)
% Validate a patch structure

if ~isstruct(value)
    error('validateInputs:invalidPatch', ...
        '%s must be a structure.', name);
end

if ~isfield(value, 'vertices') || ~isfield(value, 'faces')
    error('validateInputs:invalidPatch', ...
        '%s must have ''vertices'' and ''faces'' fields.', name);
end

if size(value.vertices, 2) ~= 3
    error('validateInputs:invalidPatch', ...
        '%s.vertices must be Nx3.', name);
end

if size(value.faces, 2) < 3
    error('validateInputs:invalidPatch', ...
        '%s.faces must have at least 3 columns.', name);
end

end

function validateColorScheme(name, value)
% Validate a color scheme table

if ~istable(value)
    error('validateInputs:invalidColorScheme', ...
        '%s must be a table.', name);
end

requiredCols = {'ion', 'color'};
missingCols = setdiff(requiredCols, value.Properties.VariableNames);
if ~isempty(missingCols)
    error('validateInputs:invalidColorScheme', ...
        '%s must contain columns: ion, color. Missing: %s', name, strjoin(missingCols, ', '));
end

end
