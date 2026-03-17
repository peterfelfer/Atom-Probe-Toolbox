function profile = massSpecProfileFromSpec(spec, colorScheme)
% MASSSPECPROFILEFROMSPEC Capture ions/ranges/colors from a mass spectrum.
%
% profile = massSpecProfileFromSpec(spec)
% profile = massSpecProfileFromSpec(spec, colorScheme)
%
% INPUT
%   spec        mass spectrum area handle, axes, or figure
%   colorScheme optional color scheme table (ion, color)
%
% OUTPUT
%   profile     scalar struct compatible with configExport/configImport
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 2
    colorScheme = [];
end

spec = resolveMassSpecHandle(spec);
ax = ancestor(spec, 'axes');
fig = ancestor(spec, 'figure');

ionTable = ionsExtractFromMassSpec(spec);
rangeTable = rangesExtractFromMassSpec(spec);
colorScheme = filterColorScheme(colorScheme, ionTable, rangeTable);

profile = struct();
profile.version = 1;
profile.schema = 'massSpecProfile';
profile.meta = struct();
profile.meta.createdAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
profile.meta.sourceFigure = string(fig.Name);

% Keep the payload table-free so YAML/JSON export works via configExport.
profile.ions = tableToStructArray(ionTable);
profile.ranges = tableToStructArray(rangeTable);
profile.colorScheme = tableToStructArray(colorScheme);

plotSettings = struct();
plotSettings.xlim = ax.XLim;
plotSettings.ylim = ax.YLim;
plotSettings.yscale = string(ax.YScale);
plotSettings.binWidth = estimateBinWidth(spec);
profile.plotSettings = plotSettings;

end

function spec = resolveMassSpecHandle(specIn)
if isgraphics(specIn, 'matlab.graphics.chart.primitive.Area')
    spec = specIn;
    return;
end

if isgraphics(specIn, 'axes')
    ax = specIn;
elseif isgraphics(specIn, 'figure')
    axesList = findobj(specIn, 'Type', 'axes');
    if isempty(axesList)
        error('massSpecProfileFromSpec:noAxes', 'No axes found in provided figure.');
    end
    ax = axesList(1);
else
    error('massSpecProfileFromSpec:invalidHandle', ...
        'Input must be a mass spectrum area handle, axes, or figure.');
end

spec = findMassSpectrumInAxes(ax);
if isempty(spec)
    error('massSpecProfileFromSpec:noMassSpectrum', ...
        'Could not resolve a mass spectrum area plot.');
end
end

function spec = findMassSpectrumInAxes(ax)
spec = [];
plots = ax.Children;
for i = 1:numel(plots)
    try
        if isfield(plots(i).UserData, 'plotType') && plots(i).UserData.plotType == "massSpectrum"
            spec = plots(i);
            return;
        end
    catch
    end
end
for i = 1:numel(plots)
    if isa(plots(i), 'matlab.graphics.chart.primitive.Area')
        spec = plots(i);
        return;
    end
end
end

function out = tableToStructArray(tbl)
if isempty(tbl) || ~istable(tbl)
    out = struct.empty(0, 1);
    return;
end
out = table2struct(tbl);
end

function colorSchemeOut = filterColorScheme(colorSchemeIn, ionTable, rangeTable)
colorSchemeOut = colorSchemeIn;
if isempty(colorSchemeOut) || ~istable(colorSchemeOut) || ...
        ~all(ismember({'ion', 'color'}, colorSchemeOut.Properties.VariableNames))
    return;
end

used = strings(0, 1);

if ~isempty(ionTable) && istable(ionTable) && ismember('ionName', ionTable.Properties.VariableNames)
    ionNames = string(ionTable.ionName);
    used = [used; normalizeNamesForColorScheme(ionNames)]; %#ok<AGROW>
end

if ~isempty(rangeTable) && istable(rangeTable) && ismember('rangeName', rangeTable.Properties.VariableNames)
    rangeNames = string(rangeTable.rangeName);
    used = [used; normalizeNamesForColorScheme(rangeNames)]; %#ok<AGROW>
end

if any(string(colorSchemeOut.ion) == "background")
    used = [used; "background"]; %#ok<AGROW>
end

used = unique(used(used ~= ""), 'stable');
if isempty(used)
    return;
end

mask = ismember(string(colorSchemeOut.ion), used);
colorSchemeOut = colorSchemeOut(mask, :);
end

function namesOut = normalizeNamesForColorScheme(namesIn)
namesOut = strings(numel(namesIn), 1);
for i = 1:numel(namesIn)
    name = string(namesIn(i));
    if strlength(name) == 0
        namesOut(i) = "";
        continue;
    end
    try
        name = string(ionConvertMode(categorical(name), 'atomic'));
    catch
    end
    namesOut(i) = strtrim(name);
end
end

function bw = estimateBinWidth(spec)
if isempty(spec) || ~isgraphics(spec) || numel(spec.XData) < 2
    bw = NaN;
    return;
end
dx = diff(spec.XData(:));
dx = dx(isfinite(dx) & dx > 0);
if isempty(dx)
    bw = NaN;
else
    bw = median(dx);
end
end
