function outCat = ionConvertMode(inCat, targetMode)
% ionConvertMode converts ion categorical names between concentration modes.
%
% outCat = ionConvertMode(inCat, targetMode)
%
% Converts ion names from ionic form to either isotopic or atomic form by
% stripping charge state and/or isotope information.
%
% INPUT
% inCat:        categorical array of ion names (e.g., '56Fe ++', '16O2 +')
% targetMode:   'ionic' | 'isotopic' | 'atomic'
%               - 'ionic': no change (returns input as-is)
%               - 'isotopic': strips charge state (e.g., '56Fe ++' -> '56Fe')
%               - 'atomic': strips isotope and charge (e.g., '56Fe ++' -> 'Fe')
%
% OUTPUT
% outCat:       categorical array with converted names
%
% EXAMPLE
% ions = categorical({'56Fe ++', '54Fe ++', '56Fe +', '16O +'});
%
% % Isotopic mode - groups by isotope (ignores charge state)
% isotopic = ionConvertMode(ions, 'isotopic');
% % Result: {'56Fe', '54Fe', '56Fe', '16O'}
%
% % Atomic mode - groups by element (ignores isotope and charge)
% atomic = ionConvertMode(ions, 'atomic');
% % Result: {'Fe', 'Fe', 'Fe', 'O'}
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 2
    targetMode = 'ionic';
end

targetMode = lower(string(targetMode));

% If ionic mode, return as-is
if targetMode == "ionic"
    outCat = inCat;
    return;
end

% Get unique categories and convert each
cats = categories(inCat);
newCats = cell(size(cats));

for i = 1:numel(cats)
    ionName = cats{i};

    % Handle special cases
    if strcmpi(ionName, 'unranged') || strcmpi(ionName, 'unknown') || isempty(ionName)
        newCats{i} = ionName;
        continue;
    end

    switch targetMode
        case "isotopic"
            % Strip charge state only, keep isotopes
            newCats{i} = stripChargeState(ionName);

        case "atomic"
            % Strip both isotopes and charge state
            newCats{i} = stripIsotopesAndCharge(ionName);

        otherwise
            error('ionConvertMode:invalidMode', ...
                'Mode must be ''ionic'', ''isotopic'', or ''atomic''.');
    end
end

% Create mapping from old categories to new
% Multiple old categories may map to the same new category
outCat = inCat;
for i = 1:numel(cats)
    outCat(inCat == cats{i}) = newCats{i};
end

% Remove unused categories
outCat = removecats(outCat);

end

function name = stripChargeState(ionName)
% Remove charge state symbols (+ and -) from ion name
% Handles both plain format ('56Fe ++') and LaTeX format ('^{56}Fe^{++}')

name = ionName;

% Handle LaTeX charge state format: '^{+}', '^{++}', '^{-}', etc. (with or without space)
name = regexprep(name, '\s*\^\{[\+\-]+\}\s*$', '');

% Handle plain charge state: ' ++', ' +', ' --', etc.
name = regexprep(name, '\s+[\+\-]+\s*$', '');

% Handle charge state directly attached (no space): '++', '+' at end
name = regexprep(name, '[\+\-]+$', '');

name = strtrim(name);
end

function name = stripIsotopesAndCharge(ionName)
% Remove isotope numbers and charge state from ion name
% Handles both plain and LaTeX formats:
% Plain: '56Fe2 16O3 ++' -> 'Fe2 O3'
% LaTeX compact: '^{56}Fe_{2}^{16}O_{3}^{++}' -> 'Fe_{2}O_{3}'
% LaTeX spaced: '^{56}Fe_{2} ^{16}O_{3} ^{++}' -> 'Fe_{2} O_{3}'

% First strip charge state
ionName = stripChargeState(ionName);

% Remove LaTeX isotope superscripts: ^{56} before element
ionName = regexprep(ionName, '\^\{\d+\}', '');

% Remove plain isotope numbers before element (e.g., '56Fe' -> 'Fe')
% Only if followed by uppercase letter (element symbol)
ionName = regexprep(ionName, '(\s|^)(\d+)([A-Z])', '$1$3');

% Keep LaTeX subscript format as-is for LaTeX output: _{2} stays _{2}
% For plain format, the subscript doesn't exist

% Handle remaining leading digits (plain format without spaces)
ionName = regexprep(ionName, '^(\d+)([A-Z])', '$2');

name = strtrim(ionName);
end
