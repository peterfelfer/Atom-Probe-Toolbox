function [report, colorSchemeOut] = massSpecProfileApply(spec, profile, options)
% MASSSPECPROFILEAPPLY Apply a mass spectrum profile to a mass spectrum.
%
% [report, colorSchemeOut] = massSpecProfileApply(spec, profile)
% [report, colorSchemeOut] = massSpecProfileApply(spec, profile, 'clearExisting', true)
%
% INPUT
%   spec      mass spectrum area handle, axes, or figure
%   profile   struct created by massSpecProfileFromSpec/configImport
%
% OPTIONS
%   'colorScheme'   color scheme table override (ion, color)
%   'isotopeTable'  isotope table override
%   'clearExisting' clear current ion/range/text plots before apply
%   'applyIons'     apply ions from profile
%   'applyRanges'   apply ranges from profile
%
% OUTPUT
%   report         struct with counts and warnings
%   colorSchemeOut color scheme after merge/apply
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    spec
    profile (1,1) struct
    options.colorScheme = []
    options.isotopeTable = []
    options.clearExisting (1,1) logical = false
    options.applyIons (1,1) logical = true
    options.applyRanges (1,1) logical = true
end

spec = resolveMassSpecHandle(spec);
profile = migrateProfile(profile);

colorSchemeProfile = structArrayToColorScheme(profile.colorScheme);
if isempty(options.colorScheme)
    colorSchemeOut = colorSchemeProfile;
else
    colorSchemeOut = normalizeColorScheme(options.colorScheme);
    colorSchemeOut = mergeColorSchemes(colorSchemeOut, colorSchemeProfile);
end
colorSchemeOut = ensureBackgroundColor(colorSchemeOut);

if isempty(options.isotopeTable)
    isotopeTable = loadIsotopeTable();
else
    isotopeTable = options.isotopeTable;
end

if options.clearExisting
    clearMassSpecOverlayPlots(spec);
end

report = struct();
report.ionsAdded = 0;
report.rangesAdded = 0;
report.warnings = strings(0, 1);

if options.applyIons
    ionTable = structArrayToIonTable(profile.ions);
    for i = 1:height(ionTable)
        try
            ionName = normalizeIonName(ionTable, i);
            chargeState = normalizeChargeState(ionTable, i);
            isTracer = false;
            if ismember('isTracer', ionTable.Properties.VariableNames)
                isTracer = logical(ionTable.isTracer(i));
            end
            if isTracer && ~contains(lower(ionName), 'tracer')
                ionName = strtrim(ionName + " tracer");
            end

            colorSchemeOut = ensureIonColor(colorSchemeOut, ionName, ionTable, i);
            ionAdd(spec, char(ionName), chargeState, isotopeTable, colorSchemeOut, ...
                0, 0.01, 'most abundant', 0.1);
            report.ionsAdded = report.ionsAdded + 1;
        catch ME
            report.warnings(end+1, 1) = "Ion apply failed (" + string(i) + "): " + string(ME.message); %#ok<AGROW>
        end
    end
end

if options.applyRanges
    rangeTable = structArrayToRangeTable(profile.ranges);
    for i = 1:height(rangeTable)
        try
            if ~ismember('mcbegin', rangeTable.Properties.VariableNames) || ...
                    ~ismember('mcend', rangeTable.Properties.VariableNames)
                continue;
            end

            limits = [double(rangeTable.mcbegin(i)), double(rangeTable.mcend(i))];
            if any(~isfinite(limits)) || limits(2) <= limits(1)
                continue;
            end

            manualName = "";
            if ismember('rangeName', rangeTable.Properties.VariableNames)
                manualName = string(rangeTable.rangeName(i));
            end
            if strlength(manualName) == 0 || strcmpi(manualName, 'not assigned')
                rangeAdd(spec, colorSchemeOut, 'background', limits);
            else
                colorSchemeOut = ensureIonColor(colorSchemeOut, manualName, [], []);
                rangeAdd(spec, colorSchemeOut, char(manualName), limits);
            end
            report.rangesAdded = report.rangesAdded + 1;
        catch ME
            report.warnings(end+1, 1) = "Range apply failed (" + string(i) + "): " + string(ME.message); %#ok<AGROW>
        end
    end
end

try
    if isfield(profile, 'plotSettings') && isstruct(profile.plotSettings)
        ax = ancestor(spec, 'axes');
        s = profile.plotSettings;
        if isfield(s, 'xlim') && isnumeric(s.xlim) && numel(s.xlim) == 2
            ax.XLim = s.xlim;
        end
        if isfield(s, 'ylim') && isnumeric(s.ylim) && numel(s.ylim) == 2
            ax.YLim = s.ylim;
        end
        if isfield(s, 'yscale') && strlength(string(s.yscale)) > 0
            ax.YScale = char(string(s.yscale));
        end
    end
catch
end

try
    massSpecReorderPlot(spec);
catch
end

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
        error('massSpecProfileApply:noAxes', 'No axes found in provided figure.');
    end
    ax = axesList(1);
else
    error('massSpecProfileApply:invalidHandle', ...
        'Input must be a mass spectrum area handle, axes, or figure.');
end

spec = findMassSpectrumInAxes(ax);
if isempty(spec)
    error('massSpecProfileApply:noMassSpectrum', ...
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

function profileOut = migrateProfile(profileIn)
profileOut = profileIn;
if ~isfield(profileOut, 'schema')
    profileOut.schema = 'massSpecProfile';
end
if ~isfield(profileOut, 'version')
    profileOut.version = 1;
end
if ~isfield(profileOut, 'ions')
    if isfield(profileOut, 'ionTable')
        profileOut.ions = profileOut.ionTable;
    else
        profileOut.ions = struct.empty(0, 1);
    end
end
if ~isfield(profileOut, 'ranges')
    if isfield(profileOut, 'rangeTable')
        profileOut.ranges = profileOut.rangeTable;
    else
        profileOut.ranges = struct.empty(0, 1);
    end
end
if ~isfield(profileOut, 'colorScheme')
    profileOut.colorScheme = struct.empty(0, 1);
end
if ~isfield(profileOut, 'plotSettings')
    profileOut.plotSettings = struct();
end
end

function tbl = structArrayToIonTable(data)
if isempty(data)
    tbl = table();
    return;
end
if istable(data)
    tbl = data;
    return;
end
if isstruct(data)
    tbl = struct2table(data);
else
    tbl = table();
end
end

function tbl = structArrayToRangeTable(data)
if isempty(data)
    tbl = table();
    return;
end
if istable(data)
    tbl = data;
    return;
end
if isstruct(data)
    tbl = struct2table(data);
else
    tbl = table();
end
end

function ionName = normalizeIonName(ionTable, rowIdx)
ionName = "";
if ismember('ionName', ionTable.Properties.VariableNames)
    ionName = string(ionTable.ionName(rowIdx));
end
if strlength(ionName) == 0 && ismember('ion', ionTable.Properties.VariableNames)
    try
        ionValue = ionTable.ion(rowIdx);
        if iscell(ionValue)
            ionValue = ionValue{1};
        end
        if isstruct(ionValue) && isfield(ionValue, 'element')
            ionName = string(ionValue(1).element);
        else
            ionName = string(ionValue);
        end
    catch
    end
end
if strlength(ionName) == 0
    error('massSpecProfileApply:invalidIonName', ...
        'Could not resolve ion name from profile.');
end
end

function chargeState = normalizeChargeState(ionTable, rowIdx)
chargeState = 1;
if ismember('chargeState', ionTable.Properties.VariableNames)
    chargeState = ionTable.chargeState(rowIdx);
end
chargeState = round(double(chargeState));
if ~isfinite(chargeState) || chargeState <= 0
    chargeState = 1;
end
end

function colorScheme = ensureIonColor(colorScheme, ionName, ionTable, rowIdx)
colorScheme = normalizeColorScheme(colorScheme);
ionBase = normalizeIonForColorScheme(ionName);
if any(string(colorScheme.ion) == ionBase)
    return;
end

if nargin >= 3 && ~isempty(ionTable) && ismember('color', ionTable.Properties.VariableNames)
    try
        c = ionTable.color(rowIdx, :);
        if size(c, 2) == 3 && all(isfinite(c(:)))
            colorScheme(end+1, :) = {categorical(ionBase), c};
            return;
        end
    catch
    end
end

try
    colorScheme = colorSchemeIonAdd(colorScheme, char(ionBase), 'create');
catch
    % Last-resort fallback color.
    colorScheme(end+1, :) = {categorical(ionBase), rand(1, 3)};
end
end

function ionOut = normalizeIonForColorScheme(ionIn)
ionOut = string(ionIn);
try
    ionOut = string(ionConvertMode(categorical(ionOut), 'atomic'));
catch
end
ionOut = strtrim(ionOut);
if strlength(ionOut) == 0
    ionOut = "unknown";
end
end

function colorScheme = structArrayToColorScheme(data)
if isempty(data)
    colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
        'VariableNames', {'ion', 'color'});
    return;
end
if istable(data)
    colorScheme = normalizeColorScheme(data);
    return;
end
if isstruct(data)
    try
        colorScheme = normalizeColorScheme(struct2table(data));
    catch
        colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
            'VariableNames', {'ion', 'color'});
    end
else
    colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
        'VariableNames', {'ion', 'color'});
end
end

function colorScheme = normalizeColorScheme(colorSchemeIn)
if isempty(colorSchemeIn)
    colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
        'VariableNames', {'ion', 'color'});
    return;
end
if ~istable(colorSchemeIn)
    error('massSpecProfileApply:invalidColorScheme', ...
        'colorScheme must be a table with columns ion and color.');
end
if ~all(ismember({'ion', 'color'}, colorSchemeIn.Properties.VariableNames))
    error('massSpecProfileApply:invalidColorScheme', ...
        'colorScheme must contain columns ion and color.');
end
colorScheme = colorSchemeIn(:, {'ion', 'color'});
colorScheme.ion = categorical(string(colorScheme.ion));
colorVals = colorScheme.color;
if size(colorVals, 2) ~= 3
    error('massSpecProfileApply:invalidColorSchemeColor', ...
        'colorScheme.color must be Nx3.');
end
colorScheme.color = double(colorVals);
end

function colorSchemeOut = mergeColorSchemes(colorSchemeA, colorSchemeB)
colorSchemeOut = colorSchemeA;
if isempty(colorSchemeB)
    return;
end
ionsA = string(colorSchemeOut.ion);
for i = 1:height(colorSchemeB)
    ionName = string(colorSchemeB.ion(i));
    idx = find(ionsA == ionName, 1, 'first');
    if isempty(idx)
        colorSchemeOut(end+1, :) = colorSchemeB(i, :);
        ionsA(end+1, 1) = ionName; %#ok<AGROW>
    else
        colorSchemeOut.color(idx, :) = colorSchemeB.color(i, :);
    end
end
end

function colorScheme = ensureBackgroundColor(colorScheme)
if any(string(colorScheme.ion) == "background")
    return;
end
colorScheme(end+1, :) = {categorical("background"), [0.7 0.7 0.7]};
end

function isotopeTable = loadIsotopeTable()
data = load('isotopeTable_naturalAbundances.mat', 'isotopeTable');
isotopeTable = data.isotopeTable;
end

function clearMassSpecOverlayPlots(spec)
ax = ancestor(spec, 'axes');
plots = ax.Children;
for i = 1:numel(plots)
    h = plots(i);
    if ~isgraphics(h) || h == spec
        continue;
    end
    deleteThis = false;
    try
        if isfield(h.UserData, 'plotType')
            pt = h.UserData.plotType;
            deleteThis = any(pt == ["ion", "range", "text", "background"]);
        end
    catch
    end
    if deleteThis
        delete(h);
    end
end
end
