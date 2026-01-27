function [colorScheme, info] = colorSchemeGenerate(pos, options)
% colorSchemeGenerate creates a new color scheme from a pos table.
%
% colorScheme = colorSchemeGenerate(pos)
% [colorScheme, info] = colorSchemeGenerate(pos, options)
%
% INPUT
% pos:      pos table containing ion or atom column
% options:  struct with fields
%   ionField:        "auto" | "ion" | "atom" (default: "auto")
%   candidateCount:  number of HSV candidates (default: 12000)
%   hueMax:          upper hue limit (default: 0.8)
%   satMin:          lower saturation limit (default: 0.25)
%   valMin:          lower value limit (default: 0.6)
%   valMax:          upper value limit (default: 1.0)
%   seedPalette:     true/false use color-blind safe seed palette (default: true)
%
% OUTPUT
% colorScheme: table with columns ion and color (RGB)
% info:        struct with ion names and seed palette size
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.ionField (1,1) string = "auto"
    options.candidateCount (1,1) double {mustBePositive} = 12000
    options.hueMax (1,1) double {mustBeNonnegative} = 0.8
    options.satMin (1,1) double {mustBeNonnegative} = 0.25
    options.valMin (1,1) double {mustBeNonnegative} = 0.6
    options.valMax (1,1) double {mustBeNonnegative} = 1.0
    options.seedPalette (1,1) logical = true
end

ionNames = extractIonNames(pos, options.ionField);
numIon = numel(ionNames);
if numIon == 0
    error('colorSchemeGenerate:emptyIonList', ...
        'No ions found in the provided table.');
end

candidateRGB = sampleHSVColors(options.candidateCount, options.hueMax, ...
    options.satMin, options.valMin, options.valMax);

if options.seedPalette
    seedRGB = okabeItoPalette();
else
    seedRGB = zeros(0, 3);
end

colors = selectColorsMaximin(candidateRGB, numIon, seedRGB);

colorScheme = table(categorical(ionNames), colors, 'VariableNames', {'ion','color'});

info = struct();
info.ionNames = ionNames;
info.seedCount = min(size(seedRGB, 1), numIon);
end

function ionNames = extractIonNames(pos, ionField)
    field = ionField;
    if field == "auto"
        if ismember('ion', pos.Properties.VariableNames)
            field = "ion";
        elseif ismember('atom', pos.Properties.VariableNames)
            field = "atom";
        else
            error('colorSchemeGenerate:missingIonField', ...
                'pos must contain ion or atom column.');
        end
    end

    if ~ismember(char(field), pos.Properties.VariableNames)
        error('colorSchemeGenerate:missingIonField', ...
            'pos does not contain the specified ion field.');
    end

    ionValues = pos.(field);
    ionValues = string(ionValues);
    ionValues = ionValues(~ismissing(ionValues) & ionValues ~= "");
    ionNames = unique(ionValues, 'stable');
end

function colors = selectColorsMaximin(candidateRGB, numColors, seedRGB)
    if numColors <= size(seedRGB, 1)
        colors = seedRGB(1:numColors, :);
        return;
    end

    chosen = seedRGB;
    remaining = candidateRGB;
    weightSets = [1 1 1; 1 0.25 1; 1 1 0.25];

    while size(chosen, 1) < numColors
        nextIdx = selectNextColor(remaining, chosen, weightSets);
        chosen(end+1, :) = remaining(nextIdx, :);
        remaining(nextIdx, :) = [];
    end

    colors = chosen(1:numColors, :);
end

function idx = selectNextColor(candidates, chosen, weightSets)
    if isempty(chosen)
        idx = 1;
        return;
    end

    candLab = rgb2lab(candidates);
    chosenLab = rgb2lab(chosen);

    minDist = inf(size(candidates, 1), 1);
    for w = 1:size(weightSets, 1)
        weights = weightSets(w, :);
        dist = pdist2(chosenLab .* weights, candLab .* weights, "euclidean");
        minDist = min(minDist, min(dist, [], 1)');
    end

    [~, idx] = max(minDist);
end

function colors = sampleHSVColors(n, hueMax, satMin, valMin, valMax)
    hueMax = max(0.1, min(1, hueMax));
    satMin = max(0, min(1, satMin));
    valMin = max(0, min(1, valMin));
    valMax = max(valMin, min(1, valMax));

    nH = max(4, ceil(n^(1/3)));
    nS = max(4, ceil(n^(1/3)));
    nV = max(2, ceil(n / (nH * nS)));

    h = linspace(0, hueMax, nH);
    s = linspace(satMin, 1, nS);
    v = linspace(valMin, valMax, nV);

    [H, S, V] = ndgrid(h, s, v);
    hsv = [H(:), S(:), V(:)];
    if size(hsv, 1) > n
        hsv = hsv(1:n, :);
    end
    colors = hsv2rgb(hsv);
end

function colors = okabeItoPalette()
    colors = [ ...
        0.9020 0.6240 0.0000; ... % orange
        0.3370 0.7060 0.9140; ... % sky blue
        0.0000 0.6200 0.4510; ... % bluish green
        0.9410 0.8940 0.2590; ... % yellow
        0.0000 0.4470 0.6980; ... % blue
        0.8350 0.3690 0.0000; ... % vermillion
        0.8000 0.4750 0.6550; ... % reddish purple
        0.4000 0.4000 0.4000; ... % gray
        ];
end
