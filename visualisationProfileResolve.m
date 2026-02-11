function [resolvedProfile, info] = visualisationProfileResolve(pos, colorScheme, visualisationProfile, options)
% VISUALISATIONPROFILERESOLVE Resolve a profile against a new POS dataset.
%
% [resolvedProfile, info] = visualisationProfileResolve(pos, colorScheme, visualisationProfile)

arguments
    pos table
    colorScheme
    visualisationProfile (1,1) struct
    options.speciesPolicy (1,1) string = ""
end

profile = visualisationProfileMigrate(visualisationProfile);
visualisationProfileValidate(profile);

settings = profile.settings;
policies = profile.policies;

if strlength(options.speciesPolicy) > 0
    speciesPolicy = lower(options.speciesPolicy);
else
    speciesPolicy = lower(string(policies.speciesPolicy));
end
if strlength(speciesPolicy) == 0
    speciesPolicy = "expand";
end

mode = normalizeMode(settings.mode, pos);
splitIsotope = logical(settings.splitIsotope);
splitCharge = logical(settings.splitCharge);
showUnranged = logical(settings.showUnranged);

[groupNames, baseNames, displayNames, groupIndices, groupCounts] = computeGroupsLocal( ...
    pos, mode, splitIsotope, splitCharge, showUnranged);

colors = mapColorsLocal(baseNames, colorScheme);

markerSizeGlobal = double(settings.markerSizeGlobal);
if ~isfinite(markerSizeGlobal) || markerSizeGlobal <= 0
    markerSizeGlobal = 15;
end

sampleValue = double(settings.sampleValue);
if ~isfinite(sampleValue) || sampleValue < 0
    sampleValue = 1;
end
sampleFractions = initSampleFractions(groupCounts, sampleValue);

n = numel(groupNames);
visible = true(n, 1);
markerSize = repmat(markerSizeGlobal, n, 1);
markerSizeLinked = true(n, 1);

saved = profile.species;
savedNames = string(saved.name(:));

matchedMask = false(n, 1);
matchedSavedMask = false(numel(savedNames), 1);

for i = 1:n
    idxSaved = find(savedNames == displayNames(i), 1, 'first');
    if isempty(idxSaved)
        continue;
    end
    matchedMask(i) = true;
    matchedSavedMask(idxSaved) = true;

    visible(i) = logical(saved.visible(idxSaved));
    sampleFractions(i) = clamp01(double(saved.fraction(idxSaved)));
    markerSize(i) = max(0.1, double(saved.markerSize(idxSaved)));
    markerSizeLinked(i) = logical(saved.markerSizeLinked(idxSaved));

    if size(saved.color, 1) >= idxSaved && size(saved.color, 2) == 3
        c = double(saved.color(idxSaved, :));
        if all(isfinite(c))
            colors(i, :) = min(1, max(0, c));
        end
    end
end

newSpecies = displayNames(~matchedMask);
missingSpecies = savedNames(~matchedSavedMask);

switch speciesPolicy
    case "expand"
        % Keep defaults for new species.

    case "intersection"
        visible(~matchedMask) = false;

    case "strict"
        if ~isempty(newSpecies) || ~isempty(missingSpecies)
            error('visualisationProfileResolve:speciesMismatch', ...
                ['Strict species policy violated. New species: %s. Missing species: %s.'], ...
                strjoin(newSpecies, ', '), strjoin(missingSpecies, ', '));
        end

    otherwise
        error('visualisationProfileResolve:invalidSpeciesPolicy', ...
            'Unknown species policy ''%s''.', speciesPolicy);
end

% Linked species always follow the global marker size.
markerSize(markerSizeLinked) = markerSizeGlobal;

resolvedProfile = profile;
resolvedProfile.version = 1;
resolvedProfile.schema = 'visualisationProfile';
resolvedProfile.settings.mode = mode;
resolvedProfile.settings.splitIsotope = splitIsotope;
resolvedProfile.settings.splitCharge = splitCharge;
resolvedProfile.settings.showUnranged = showUnranged;
resolvedProfile.settings.markerSizeGlobal = markerSizeGlobal;
resolvedProfile.settings.sampleValue = sampleValue;

resolvedProfile.species = struct();
resolvedProfile.species.name = displayNames;
resolvedProfile.species.baseName = baseNames;
resolvedProfile.species.count = groupCounts;
resolvedProfile.species.visible = visible;
resolvedProfile.species.fraction = sampleFractions;
resolvedProfile.species.markerSize = markerSize;
resolvedProfile.species.markerSizeLinked = markerSizeLinked;
resolvedProfile.species.color = colors;

resolvedProfile.meta.resolvedAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

info = struct();
info.mode = mode;
info.speciesPolicy = speciesPolicy;
info.groupNames = groupNames;
info.displayNames = displayNames;
info.baseNames = baseNames;
info.groupCounts = groupCounts;
info.groupIndices = groupIndices;
info.matchedSpecies = displayNames(matchedMask);
info.newSpecies = newSpecies;
info.missingSpecies = missingSpecies;

end

function mode = normalizeMode(modeIn, pos)
mode = lower(string(modeIn));
if mode == "ion"
    mode = "ionic";
elseif mode == "atom"
    mode = "atomic";
end

if mode ~= "ionic" && mode ~= "atomic"
    if ismember('ion', pos.Properties.VariableNames)
        mode = "ionic";
    else
        mode = "atomic";
    end
end
end

function fractions = initSampleFractions(counts, sampleValue)
if sampleValue <= 1
    fractions = repmat(sampleValue, size(counts));
else
    fractions = min(1, sampleValue ./ max(counts, 1));
end
fractions = max(0, min(1, fractions));
end

function x = clamp01(x)
x = min(1, max(0, x));
end

function [groupNames, baseNames, displayNames, groupIndices, groupCounts] = computeGroupsLocal( ...
    pos, mode, splitIsotope, splitCharge, showUnranged)

if mode == "ionic"
    if ismember('ion', pos.Properties.VariableNames)
        base = string(pos.ion);
    elseif ismember('atom', pos.Properties.VariableNames)
        base = string(pos.atom);
    else
        error('visualisationProfileResolve:missingGroupField', ...
            'pos must contain ion or atom column.');
    end
else
    if ismember('atom', pos.Properties.VariableNames)
        base = string(pos.atom);
    elseif ismember('ion', pos.Properties.VariableNames)
        base = string(pos.ion);
    else
        error('visualisationProfileResolve:missingGroupField', ...
            'pos must contain ion or atom column.');
    end
end

keep = true(size(base));
if showUnranged
    base(ismissing(base)) = "unranged";
else
    keep = ~ismissing(base);
end

baseWork = base(keep);
rowIdx = find(keep);

groupKey = baseWork;
if splitIsotope && ismember('isotope', pos.Properties.VariableNames)
    iso = string(pos.isotope(rowIdx));
    iso(ismissing(iso)) = "NaN";
    groupKey = groupKey + "-" + iso;
end
if splitCharge && ismember('chargeState', pos.Properties.VariableNames)
    cs = pos.chargeState(rowIdx);
    csString = strings(size(cs));
    csString(:) = "";
    valid = ~isnan(cs);
    for k = 1:numel(cs)
        if valid(k)
            n = round(cs(k));
            if n > 0
                csString(k) = repmat("+", 1, n);
            end
        end
    end
    groupKey = groupKey + csString;
end

[groupNames, ~, groupIdx] = unique(groupKey, 'stable');
baseNames = strings(size(groupNames));
displayNames = strings(size(groupNames));
groupIndices = cell(numel(groupNames), 1);
groupCounts = zeros(numel(groupNames), 1);

for i = 1:numel(groupNames)
    idxLocal = find(groupIdx == i);
    idxGlobal = rowIdx(idxLocal);
    groupIndices{i} = idxGlobal;
    groupCounts(i) = numel(idxGlobal);
    baseNames(i) = baseWork(idxLocal(1));
    displayNames(i) = formatDisplayNameLocal(baseWork(idxLocal(1)), pos, idxGlobal(1), splitIsotope, splitCharge);
end

end

function displayName = formatDisplayNameLocal(baseName, pos, rowIdx, splitIsotope, splitCharge)
name = string(baseName);
if name == "" || name == "missing"
    name = "unranged";
end
if splitCharge
    name = stripChargeFromNameLocal(name);
end

isoPart = "";
if splitIsotope && ismember('isotope', pos.Properties.VariableNames)
    iso = pos.isotope(rowIdx);
    if ~ismissing(iso)
        isoPart = string(iso);
    end
end

chargePart = "";
if splitCharge && ismember('chargeState', pos.Properties.VariableNames)
    cs = pos.chargeState(rowIdx);
    if ~isnan(cs)
        n = round(cs);
        if n > 0
            chargePart = repmat("+", 1, n);
        end
    end
end

displayName = isoPart + name + chargePart;
end

function nameOut = stripChargeFromNameLocal(nameIn)
nameOut = string(nameIn);
nameOut = regexprep(nameOut, '(\\d*\\+)+$', '');
nameOut = regexprep(nameOut, '\\++$', '');
nameOut = strtrim(nameOut);
end

function colors = mapColorsLocal(baseNames, colorScheme)
n = numel(baseNames);
defaultColors = lines(max(n, 7));
colors = zeros(n, 3);

hasScheme = istable(colorScheme) && ...
    all(ismember({'ion', 'color'}, colorScheme.Properties.VariableNames));
if hasScheme
    ionNames = string(colorScheme.ion);
else
    ionNames = string.empty(0, 1);
end

for i = 1:n
    name = baseNames(i);
    if hasScheme
        idx = find(ionNames == name, 1, 'first');
        if isempty(idx) && name == "unranged"
            idx = find(ionNames == "background", 1, 'first');
        end
        if ~isempty(idx)
            colors(i, :) = colorScheme.color(idx, :);
            continue;
        end
    end
    colors(i, :) = defaultColors(mod(i - 1, size(defaultColors, 1)) + 1, :);
end

end
