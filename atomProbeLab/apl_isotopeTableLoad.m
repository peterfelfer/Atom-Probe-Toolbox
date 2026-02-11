
function isotopeTable = apl_isotopeTableLoad(isotopeTableIn)
% apl_isotopeTableLoad loads an isotope table from table or .mat file.
% Default is isotopeTable_naturalAbundances.mat in toolbox root.
%
% INPUT
%   isotopeTableIn: table | struct | char | string | []
%
% OUTPUT
%   isotopeTable: table
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

if nargin < 1 || isempty(isotopeTableIn)
    thisDir = fileparts(mfilename('fullpath'));
    toolboxRoot = fileparts(thisDir);
    isotopeTableIn = fullfile(toolboxRoot, 'isotopeTable_naturalAbundances.mat');
end

if istable(isotopeTableIn)
    isotopeTable = isotopeTableIn;
    return;
end

if isstruct(isotopeTableIn)
    f = fieldnames(isotopeTableIn);
    for k = 1:numel(f)
        if istable(isotopeTableIn.(f{k}))
            isotopeTable = isotopeTableIn.(f{k});
            return;
        end
    end
    error('apl_isotopeTableLoad:invalidStruct', 'No table found in struct input.');
end

if isstring(isotopeTableIn)
    isotopeTableIn = char(isotopeTableIn);
end

if ischar(isotopeTableIn)
    if ~exist(isotopeTableIn, 'file')
        error('apl_isotopeTableLoad:notFound', 'Isotope table file not found: %s', isotopeTableIn);
    end
    data = load(isotopeTableIn);
    f = fieldnames(data);
    for k = 1:numel(f)
        if istable(data.(f{k}))
            isotopeTable = data.(f{k});
            return;
        end
    end
    error('apl_isotopeTableLoad:noTable', 'No table found in %s', isotopeTableIn);
end

error('apl_isotopeTableLoad:invalidInput', 'Unsupported isotopeTable input.');
end
