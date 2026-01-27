function [colorScheme, info] = colorSchemeUpdate(colorScheme, source, varargin)
% colorSchemeUpdate updates a color scheme with new ions from plots or tables.
%
% colorScheme = colorSchemeUpdate(colorScheme, source)
% [colorScheme, info] = colorSchemeUpdate(colorScheme, source, 'overwrite', true)
%
% INPUT
% colorScheme: table with columns ion and color
% source:      graphics handle (mass spectrum) or pos table
%
% OPTIONS (name-value pairs)
% overwrite:       overwrite existing colors (default: false)
% candidateCount:  number of HSV samples used to generate new colors (default: 10000)
% hueMax:          upper hue limit for HSV sampling (default: 0.8)
% satMin:          lower saturation limit for HSV sampling (default: 0.2)
% value:           HSV value (default: 1)
%
% OUTPUT
% colorScheme: updated color scheme
% info:        struct with added/updated ions
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 2 || isempty(source)
    source = gcf;
end

opts = struct('overwrite', false, ...
    'candidateCount', 10000, ...
    'hueMax', 0.8, ...
    'satMin', 0.2, ...
    'value', 1);
if ~isempty(varargin)
    opts = parseOptions(opts, varargin{:});
end

if ~istable(colorScheme) || ~all(ismember({'ion','color'}, colorScheme.Properties.VariableNames))
    error('colorSchemeUpdate:invalidColorScheme', ...
        'colorScheme must be a table with columns ion and color.');
end

if isgraphics(source)
    [sourceIons, sourceColors] = extractFromGraphics(source);
    [colorScheme, infoPlot] = mergeColorsFromSource(colorScheme, sourceIons, sourceColors, opts.overwrite);
    info = infoPlot;
    return;
end

if istable(source)
    sourceIons = extractFromPos(source);
    sourceColors = [];
else
    error('colorSchemeUpdate:invalidSource', ...
        'source must be a graphics handle or a table.');
end

[colorScheme, infoPlot] = mergeColorsFromSource(colorScheme, string.empty(0,1), [], opts.overwrite);
missing = setdiff(string(sourceIons), string(colorScheme.ion), 'stable');
if isempty(missing)
    info = infoPlot;
    info.addedIons = string.empty(0,1);
    return;
end

existingColors = colorScheme.color;
newColors = generateColorsHSV(existingColors, numel(missing), opts);

newRows = table(categorical(missing), newColors, 'VariableNames', {'ion','color'});
colorScheme = [colorScheme; newRows];

info = infoPlot;
info.addedIons = missing;
info.addedColors = newColors;
end

function [colorScheme, info] = mergeColorsFromSource(colorScheme, ionNames, colors, overwrite)
    ionNames = string(ionNames);
    valid = ionNames ~= "";
    ionNames = ionNames(valid);
    colors = colors(valid, :);

    if isempty(ionNames)
        info = struct('addedIons', string.empty(0,1), 'updatedIons', string.empty(0,1), ...
            'addedColors', [], 'updatedColors', []);
        return;
    end

    existingNames = string(colorScheme.ion);
    added = string.empty(0,1);
    updated = string.empty(0,1);
    addedColors = [];
    updatedColors = [];

    for i = 1:numel(ionNames)
        name = ionNames(i);
        color = colors(i, :);
        idx = find(existingNames == name, 1, 'first');
        if isempty(idx)
            newRow = table(categorical(name), color, 'VariableNames', {'ion','color'});
            colorScheme = [colorScheme; newRow];
            added(end+1,1) = name;
            addedColors(end+1,:) = color;
            existingNames(end+1,1) = name;
        elseif overwrite
            colorScheme.color(idx, :) = color;
            updated(end+1,1) = name;
            updatedColors(end+1,:) = color;
        end
    end

    info = struct('addedIons', added, 'updatedIons', updated, ...
        'addedColors', addedColors, 'updatedColors', updatedColors);
end

function [ionNames, colors] = extractFromGraphics(source)
    ax = source;
    if ~isgraphics(ax, 'axes')
        ax = ancestor(source, 'axes');
    end
    if isempty(ax) || ~isgraphics(ax)
        error('colorSchemeUpdate:invalidGraphics', ...
            'Could not resolve an axes from the provided graphics handle.');
    end

    plots = ax.Children;
    ionNames = string.empty(0,1);
    colors = zeros(0, 3);

    for pl = 1:numel(plots)
        try
            type = plots(pl).UserData.plotType;
        catch
            type = "unknown";
        end
        if type == "ion"
            name = ionFromIonPlot(plots(pl));
            color = plotColor(plots(pl));
        elseif type == "range"
            name = ionFromRangePlot(plots(pl));
            color = plotFaceColor(plots(pl));
        elseif type == "background"
            name = "background";
            color = plotFaceColor(plots(pl));
        else
            continue;
        end

        if ~isempty(name) && ~isempty(color)
            ionNames(end+1,1) = string(name);
            colors(end+1,:) = color;
        end
    end
end

function ionNames = extractFromPos(pos)
    if ismember('atom', pos.Properties.VariableNames)
        ions = pos.atom;
    elseif ismember('ion', pos.Properties.VariableNames)
        ions = pos.ion;
    else
        error('colorSchemeUpdate:missingIonField', ...
            'pos must contain ion or atom column.');
    end
    if iscategorical(ions)
        ions = string(ions);
    else
        ions = string(ions);
    end
    ions = ions(~ismissing(ions) & ions ~= "");
    ionNames = unique(ions, 'stable');
end

function name = ionFromIonPlot(h)
    name = "";
    if isfield(h.UserData, 'ion') && ~isempty(h.UserData.ion)
        try
            ionTable = h.UserData.ion{1};
            name = ionConvertName(ionTable.element);
        catch
            name = string(h.DisplayName);
        end
    elseif isprop(h, 'DisplayName')
        name = string(h.DisplayName);
    end
end

function name = ionFromRangePlot(h)
    name = "";
    if isfield(h.UserData, 'ion')
        if istable(h.UserData.ion)
            name = ionConvertName(h.UserData.ion.element);
        else
            name = string(h.UserData.ion);
        end
    elseif isprop(h, 'DisplayName')
        name = string(h.DisplayName);
    end
end

function color = plotColor(h)
    color = [];
    if isprop(h, 'Color')
        color = h.Color;
    end
    if isempty(color) || numel(color) ~= 3
        color = [];
    end
end

function color = plotFaceColor(h)
    color = [];
    if isprop(h, 'FaceColor')
        color = h.FaceColor;
    end
    if isempty(color) || numel(color) ~= 3
        color = [];
    end
end

function colors = generateColorsHSV(existingColors, nAdd, opts)
    if nAdd == 0
        colors = zeros(0,3);
        return;
    end
    existingHSV = [];
    if ~isempty(existingColors)
        existingHSV = rgb2hsv(existingColors);
    end
    candidates = sampleHSV(opts.candidateCount, opts.hueMax, opts.satMin, opts.value);

    if isempty(existingHSV)
        idx = round(linspace(1, size(candidates, 1), nAdd));
        chosenHSV = candidates(idx, :);
    else
        chosenHSV = zeros(nAdd, 3);
        for k = 1:nAdd
            dist = pdist2(existingHSV, candidates, "euclidean");
            md = min(dist, [], 1);
            [~, idx] = max(md);
            chosenHSV(k, :) = candidates(idx, :);
            existingHSV = [existingHSV; candidates(idx, :)];
            candidates(idx, :) = [];
        end
    end
    colors = hsv2rgb(chosenHSV);
end

function candidates = sampleHSV(n, hueMax, satMin, value)
    nH = max(2, ceil(sqrt(n)));
    nS = max(2, ceil(n / nH));
    h = linspace(0, hueMax, nH);
    s = linspace(satMin, 1, nS);
    [H, S] = meshgrid(h, s);
    candidates = [H(:), S(:), repmat(value, numel(H), 1)];
    if size(candidates, 1) > n
        candidates = candidates(1:n, :);
    end
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('colorSchemeUpdate:invalidOptions', ...
            'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        value = varargin{k+1};
        switch name
            case "overwrite"
                opts.overwrite = logical(value);
            case "candidatecount"
                opts.candidateCount = max(10, round(value));
            case "huemax"
                opts.hueMax = max(0, min(1, value));
            case "satmin"
                opts.satMin = max(0, min(1, value));
            case "value"
                opts.value = max(0, min(1, value));
            otherwise
                error('colorSchemeUpdate:invalidOption', ...
                    'Unknown option "%s".', name);
        end
    end
end
