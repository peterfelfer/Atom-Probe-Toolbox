function [excess, info] = posCalculateInterfacialExcess(pos, dist, species, varargin)
% posCalculateInterfacialExcess calculates interfacial excess from distance data
%
% [excess, info] = posCalculateInterfacialExcess(pos, dist, species)
% [excess, info] = posCalculateInterfacialExcess(pos, dist, species, 'Name', Value)
%
% This function calculates the interfacial excess of specified species
% relative to an interface defined by a distance array. The distance can
% be to any feature: mesh interface, plane, line, point, isosurface, etc.
%
% INPUT
% pos:      pos table with ranged atoms (must have 'atom' or 'ion' column)
%
% dist:     distance array (same length as pos) - distance to interface/feature
%           Positive and negative values indicate sides of the interface.
%           The interface location is where |dist| is minimum.
%
% species:  cell array of species to calculate excess for
%           Examples: {'Al', 'Ti'}      - atomic mode
%                     {'Al++', 'Ti+'}   - ionic mode
%                     {'27Al', '48Ti'}  - isotopic mode
%
% OPTIONS
% 'mode'        'atomic' | 'ionic' | 'isotopic' (default: auto-detect from species)
%               - atomic: group by element (ignores isotope and charge)
%               - ionic: each ion species separately
%               - isotopic: group by isotope (ignores charge)
%
% 'area'        interface area in nm^2 (optional)
%               If provided, excess is reported in at/nm^2
%
% 'interface'   patch struct with .vertices and .faces (optional)
%               If provided, area is calculated automatically from mesh
%
% 'side'        'both' | 'positive' | 'negative' (default: 'both')
%               - 'both': fit linear regions on both sides of interface
%               - 'positive': surface/cluster mode, fit only positive side
%               - 'negative': fit only negative side
%
% 'gui'         true | false (default: true)
%               If true, launches interactive GUI for adjusting fit limits
%
% 'threshold'   fraction of max derivative for auto fit limits (default: 0.2)
%               Only used when gui=false
%
% OUTPUT
% excess:   table with columns:
%           - species: species name
%           - excess_atoms: excess in number of atoms
%           - excess_per_nm2: excess in at/nm^2 (NaN if no area)
%           - partial1: partial excess on negative side
%           - partial2: partial excess on positive side
%           - fit_limits: [lim1 lim2 lim3 lim4] indices
%
% info:     struct with:
%           - cumulative: cumulative curves for each species
%           - dist_sorted: sorted distance array
%           - sortIdx: sort indices
%           - interfaceLoc: interface location index
%           - area: interface area (or NaN)
%           - fits: linear fit coefficients
%           - guiHandle: handle to GUI figure (if gui=true)
%
% EXAMPLE
% % Calculate Al and Ti excess relative to a grain boundary mesh
% dist = patchPointDistance(pos, interface);  % or any distance calculation
% [excess, info] = posCalculateInterfacialExcess(pos, dist, {'Al', 'Ti'}, ...
%     'interface', interface);
%
% % Calculate excess without GUI (automatic fitting)
% [excess, info] = posCalculateInterfacialExcess(pos, dist, {'Al'}, ...
%     'gui', false, 'area', 125.4);
%
% SEE ALSO
% patchCreateExcessValue, patchCreateProxigram, patchPointDistance
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Input validation
if ~istable(pos)
    error('posCalculateInterfacialExcess:invalidInput', ...
        'pos must be a table.');
end

if ~isnumeric(dist) || length(dist) ~= height(pos)
    error('posCalculateInterfacialExcess:invalidInput', ...
        'dist must be a numeric array with same length as pos.');
end

if ~iscell(species)
    if ischar(species) || isstring(species)
        species = {char(species)};
    else
        error('posCalculateInterfacialExcess:invalidInput', ...
            'species must be a cell array of strings.');
    end
end

%% Parse options
options = struct(...
    'mode', '', ...
    'area', NaN, ...
    'interface', [], ...
    'side', 'both', ...
    'gui', true, ...
    'threshold', 0.2);

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "mode"
            options.mode = lower(char(value));
        case "area"
            options.area = value;
        case "interface"
            options.interface = value;
        case "side"
            options.side = lower(char(value));
        case "gui"
            options.gui = logical(value);
        case "threshold"
            options.threshold = value;
        otherwise
            error('posCalculateInterfacialExcess:invalidOption', ...
                'Unknown option "%s".', name);
    end
end

%% Determine mode from species format if not specified
mode = options.mode;
if isempty(mode)
    mode = detectModeFromSpecies(species);
end

%% Calculate area from interface if provided
area = options.area;
if ~isnan(area)
    % Area provided directly, use it
elseif ~isempty(options.interface) && isstruct(options.interface)
    % Calculate area from mesh
    area = calculateMeshArea(options.interface);
else
    % No area available
    area = NaN;
end

%% Get species labels from pos table
[speciesLabels, speciesColumn] = getSpeciesLabels(pos, mode);

%% Sort by distance
dist = dist(:);
[dist_sorted, sortIdx] = sort(dist);
speciesLabels_sorted = speciesLabels(sortIdx);

numAtoms = length(dist);

%% Find interface location (where |dist| is minimum)
[~, interfaceLoc] = min(abs(dist_sorted));

%% Build cumulative curves for each species
numSpecies = length(species);
cumulative = zeros(numAtoms, numSpecies);
speciesFound = false(numSpecies, 1);

for s = 1:numSpecies
    speciesName = species{s};

    % Match species in the labels
    isMatch = matchSpecies(speciesLabels_sorted, speciesName, mode);

    if any(isMatch)
        speciesFound(s) = true;
        cumulative(:, s) = cumsum(double(isMatch));
    else
        warning('posCalculateInterfacialExcess:speciesNotFound', ...
            'Species "%s" not found in pos table.', speciesName);
    end
end

%% Determine initial fit limits
isSurface = strcmp(options.side, 'positive') || strcmp(options.side, 'negative');
limits = calculateInitialLimits(numAtoms, interfaceLoc, isSurface, options.side);

%% Build info struct
info = struct();
info.cumulative = cumulative;
info.dist_sorted = dist_sorted;
info.sortIdx = sortIdx;
info.interfaceLoc = interfaceLoc;
info.area = area;
info.species = species;
info.speciesFound = speciesFound;
info.mode = mode;
info.side = options.side;
info.limits = limits;
info.guiHandle = [];

%% Launch GUI or auto-calculate
if options.gui
    % Launch interactive GUI
    [excess, info] = interfacialExcessGUI(info);
else
    % Auto-calculate using threshold method
    [excess, info] = calculateExcessAutomatic(info, options.threshold);
end

end

%% ========== Helper Functions ==========

function mode = detectModeFromSpecies(species)
% Detect mode from species format
% - Contains + or -: ionic
% - Starts with number: isotopic
% - Otherwise: atomic

    mode = 'atomic';
    for s = 1:length(species)
        sp = species{s};
        if contains(sp, '+') || contains(sp, '-')
            mode = 'ionic';
            return;
        end
        if ~isempty(regexp(sp, '^\d', 'once'))
            mode = 'isotopic';
            return;
        end
    end
end

function area = calculateMeshArea(interface)
% Calculate total area of a triangular mesh

    if ~isfield(interface, 'vertices') || ~isfield(interface, 'faces')
        area = NaN;
        return;
    end

    v = interface.vertices;
    f = interface.faces;

    % Calculate area of each triangle
    v1 = v(f(:,1), :);
    v2 = v(f(:,2), :);
    v3 = v(f(:,3), :);

    % Cross product method
    crossProd = cross(v2 - v1, v3 - v1, 2);
    triangleAreas = 0.5 * sqrt(sum(crossProd.^2, 2));

    area = sum(triangleAreas);
end

function [labels, columnName] = getSpeciesLabels(pos, mode)
% Get species labels from pos table based on mode

    hasAtom = ismember('atom', pos.Properties.VariableNames);
    hasIon = ismember('ion', pos.Properties.VariableNames);

    switch mode
        case 'atomic'
            if hasAtom
                labels = pos.atom;
                columnName = 'atom';
            elseif hasIon
                labels = ionConvertMode(pos.ion, 'atomic');
                columnName = 'ion';
            else
                error('posCalculateInterfacialExcess:missingColumn', ...
                    'pos table must have ''atom'' or ''ion'' column for atomic mode.');
            end

        case 'ionic'
            if hasIon
                labels = pos.ion;
                columnName = 'ion';
            else
                error('posCalculateInterfacialExcess:missingColumn', ...
                    'pos table must have ''ion'' column for ionic mode.');
            end

        case 'isotopic'
            if hasIon
                labels = ionConvertMode(pos.ion, 'isotopic');
                columnName = 'ion';
            elseif hasAtom && ismember('isotope', pos.Properties.VariableNames)
                labels = buildIsotopeLabels(pos);
                columnName = 'atom';
            else
                error('posCalculateInterfacialExcess:missingColumn', ...
                    'pos table must have ''ion'' column or ''atom''+''isotope'' for isotopic mode.');
            end

        otherwise
            error('posCalculateInterfacialExcess:invalidMode', ...
                'Mode must be ''atomic'', ''ionic'', or ''isotopic''.');
    end

    % Ensure categorical
    if ~iscategorical(labels)
        labels = categorical(labels);
    end
end

function labels = buildIsotopeLabels(pos)
% Build isotope labels from atom and isotope columns

    atomCol = pos.atom;
    isoCol = pos.isotope;

    if iscategorical(atomCol)
        atomStr = string(atomCol);
    else
        atomStr = string(atomCol);
    end

    if isnumeric(isoCol)
        isoStr = strings(size(isoCol));
        validIso = ~isnan(isoCol) & isoCol > 0;
        isoStr(validIso) = string(isoCol(validIso));
    else
        isoStr = string(isoCol);
    end

    labels = isoStr + atomStr;
    labels = categorical(labels);
end

function isMatch = matchSpecies(labels, speciesName, mode)
% Match species name against labels

    % Convert to string for comparison
    labelStr = string(labels);
    speciesStr = string(speciesName);

    % Direct match
    isMatch = labelStr == speciesStr;

    % If no match and mode is atomic, try case-insensitive
    if ~any(isMatch) && strcmp(mode, 'atomic')
        isMatch = strcmpi(labelStr, speciesStr);
    end

    % If still no match and species has charge, try matching without charge
    if ~any(isMatch) && (contains(speciesStr, '+') || contains(speciesStr, '-'))
        baseSpecies = regexprep(speciesStr, '[+-]+$', '');
        isMatch = labelStr == baseSpecies | strcmpi(labelStr, baseSpecies);
    end
end

function limits = calculateInitialLimits(numAtoms, interfaceLoc, isSurface, side)
% Calculate initial fit limit positions

    delta = round(numAtoms / 10);

    if isSurface
        if strcmp(side, 'positive')
            % Surface mode: fit only positive side
            limits = [1, 1, interfaceLoc + delta, numAtoms];
            if limits(3) > numAtoms
                limits(3) = round((interfaceLoc + numAtoms) / 2);
            end
        else
            % Negative side only
            limits = [1, interfaceLoc - delta, numAtoms, numAtoms];
            if limits(2) < 1
                limits(2) = round(interfaceLoc / 2);
            end
        end
    else
        % Both sides
        limL = max(1, interfaceLoc - delta);
        limU = min(numAtoms, interfaceLoc + delta);

        % Ensure limits are within reasonable range
        if limL < round(interfaceLoc / 2)
            limL = round(interfaceLoc / 2);
        end
        if (limU - interfaceLoc) > (numAtoms - interfaceLoc) / 2
            limU = interfaceLoc + round((numAtoms - interfaceLoc) / 2);
        end

        limits = [1, limL, limU, numAtoms];
    end

    limits = round(limits);
    limits = max(1, min(numAtoms, limits));
end

function [excess, info] = calculateExcessAutomatic(info, threshold)
% Calculate excess using automatic threshold-based fitting

    numAtoms = length(info.dist_sorted);
    numSpecies = length(info.species);
    limits = info.limits;
    area = info.area;
    interfaceLoc = info.interfaceLoc;
    cumulative = info.cumulative;
    isSurface = strcmp(info.side, 'positive') || strcmp(info.side, 'negative');

    % Prepare output
    speciesNames = info.species(:);
    excess_atoms = zeros(numSpecies, 1);
    excess_per_nm2 = zeros(numSpecies, 1);
    partial1 = zeros(numSpecies, 1);
    partial2 = zeros(numSpecies, 1);
    fit_limits = repmat(limits, numSpecies, 1);

    x = (1:numAtoms)';

    for s = 1:numSpecies
        if ~info.speciesFound(s)
            excess_atoms(s) = NaN;
            excess_per_nm2(s) = NaN;
            partial1(s) = NaN;
            partial2(s) = NaN;
            continue;
        end

        cum = cumulative(:, s);

        % Linear regression
        if isSurface
            % Surface mode: one fit only
            if strcmp(info.side, 'positive')
                a1 = 0; b1 = 0;
                fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
                a2 = fit2(1); b2 = fit2(2);
            else
                fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
                a1 = fit1(1); b1 = fit1(2);
                a2 = 0; b2 = cum(end);
            end
        else
            % Both sides
            fit1 = polyfit(x(limits(1):limits(2)), cum(limits(1):limits(2)), 1);
            a1 = fit1(1); b1 = fit1(2);
            fit2 = polyfit(x(limits(3):limits(4)), cum(limits(3):limits(4)), 1);
            a2 = fit2(1); b2 = fit2(2);
        end

        % Calculate intersections at interface
        intersect1 = a1 * interfaceLoc + b1;
        intersect2 = a2 * interfaceLoc + b2;

        % Excess
        excess_atoms(s) = intersect2 - intersect1;
        if ~isnan(area) && area > 0
            excess_per_nm2(s) = excess_atoms(s) / area;
        else
            excess_per_nm2(s) = NaN;
        end

        % Partial excesses
        cumAtInterface = cum(round(interfaceLoc));
        partial1(s) = cumAtInterface - intersect1;
        partial2(s) = intersect2 - cumAtInterface;

        if ~isnan(area) && area > 0
            partial1(s) = partial1(s) / area;
            partial2(s) = partial2(s) / area;
        end
    end

    % Build output table
    excess = table(speciesNames, excess_atoms, excess_per_nm2, partial1, partial2, fit_limits, ...
        'VariableNames', {'species', 'excess_atoms', 'excess_per_nm2', 'partial1', 'partial2', 'fit_limits'});

    info.limits = limits;
    info.excess = excess;
end
