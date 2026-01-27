function [spec, info] = colorSchemeApplyToMassSpec(spec, colorScheme)
% colorSchemeApplyToMassSpec applies a color scheme to a mass spectrum plot.
%
% spec = colorSchemeApplyToMassSpec(spec, colorScheme)
% [spec, info] = colorSchemeApplyToMassSpec(spec, colorScheme)
%
% INPUT
% spec:        graphics handle (area plot, axes, or figure) of mass spectrum
% colorScheme: table with columns ion and color
%
% OUTPUT
% spec:        same handle as input
% info:        struct with applied and missing ions
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 1 || isempty(spec)
    spec = gcf;
end

if ~istable(colorScheme) || ~all(ismember({'ion','color'}, colorScheme.Properties.VariableNames))
    error('colorSchemeApplyToMassSpec:invalidColorScheme', ...
        'colorScheme must be a table with columns ion and color.');
end

ax = resolveAxes(spec);

csNames = string(colorScheme.ion);
csColors = colorScheme.color;
if size(csColors, 2) ~= 3
    error('colorSchemeApplyToMassSpec:invalidColorScheme', ...
        'colorScheme.color must be Nx3 RGB values.');
end

function ax = resolveAxes(spec)
    ax = gobjects(0, 1);

    if isgraphics(spec, 'axes')
        ax = spec;
    elseif isgraphics(spec, 'figure')
        axList = findobj(spec, 'Type', 'axes');
        if ~isempty(axList)
            ax = axList(1);
        end
    elseif isgraphics(spec)
        ax = ancestor(spec, 'axes');
    end

    if isempty(ax) || ~isgraphics(ax, 'axes')
        if isstruct(spec) && isfield(spec, 'Parent')
            ax = spec.Parent;
        elseif isobject(spec) && isprop(spec, 'Parent')
            ax = spec.Parent;
        end
    end

    if isempty(ax) || ~isgraphics(ax, 'axes')
        error('colorSchemeApplyToMassSpec:invalidHandle', ...
            'Could not resolve an axes from the provided handle.');
    end
end

plots = ax.Children;
missing = string.empty(0, 1);
applied = string.empty(0, 1);
appliedHandles = gobjects(0, 1);

for pl = 1:numel(plots)
    h = plots(pl);
    try
        plotType = h.UserData.plotType;
    catch
        plotType = "unknown";
    end

    ionName = "";
    colorTarget = [];
    if plotType == "ion"
        ionName = ionFromIonPlot(h);
        colorTarget = "line";
    elseif plotType == "range"
        ionName = ionFromRangePlot(h);
        colorTarget = "face";
    elseif plotType == "background"
        ionName = "background";
        colorTarget = "face";
    else
        continue;
    end

    if ionName == ""
        continue;
    end

    idx = find(csNames == ionName, 1, 'first');
    if isempty(idx)
        missing(end+1, 1) = ionName;
        continue;
    end

    color = csColors(idx, :);
    if colorTarget == "face" && isprop(h, 'FaceColor')
        h.FaceColor = color;
    elseif colorTarget == "line" && isprop(h, 'Color')
        h.Color = color;
    end

    applied(end+1, 1) = ionName;
    appliedHandles(end+1, 1) = h;
end

info = struct();
info.appliedIons = unique(applied, 'stable');
info.missingIons = unique(missing, 'stable');
info.handles = appliedHandles;
end

function name = ionFromIonPlot(h)
    name = "";
    if isfield(h.UserData, 'ion') && ~isempty(h.UserData.ion)
        try
            if iscell(h.UserData.ion)
                ionTable = h.UserData.ion{1};
            else
                ionTable = h.UserData.ion;
            end
            if istable(ionTable)
                name = ionConvertName(ionTable.element);
            else
                name = string(ionTable);
            end
        catch
            name = string(h.DisplayName);
        end
    elseif isprop(h, 'DisplayName')
        name = string(h.DisplayName);
    end
    name = strtrim(name);
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
    name = strtrim(name);
end
