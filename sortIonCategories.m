function [catsOut, varargout] = sortIonCategories(cats, sortBy, varargin)
% sortIonCategories sorts ion/isotope/atom category names.
%
% catsOut = sortIonCategories(cats, sortBy)
% [catsOut, data1, data2, ...] = sortIonCategories(cats, sortBy, data1, data2, ...)
%
% Sorts category names by atomic number or mass-to-charge ratio.
% Optionally reorders associated data arrays using the same sort order.
%
% INPUT
% cats:     cell array of category names (e.g., {'Fe', 'O', 'C', 'unranged'})
% sortBy:   'atomic' | 'weight' | 'none' | 'name'
%           - 'atomic': sort by atomic number (primary) and isotope (secondary)
%           - 'weight': sort by mass-to-charge ratio
%           - 'name': alphabetical sort
%           - 'none': no sorting, return as-is
% data1, data2, ...: optional data arrays to reorder (must be same length as cats)
%
% OUTPUT
% catsOut:  sorted category names
% data1, data2, ...: reordered data arrays (if provided)
%
% EXAMPLE
% cats = {'O', 'Fe', 'C', 'unranged'};
% counts = [100, 500, 50, 20];
%
% % Sort by atomic number
% [sortedCats, sortedCounts] = sortIonCategories(cats, 'atomic', counts);
% % Result: {'C', 'O', 'Fe', 'unranged'}, [50, 100, 500, 20]
%
% % Sort by mass
% [sortedCats, sortedCounts] = sortIonCategories(cats, 'weight', counts);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

catsOut = cats;
varargout = varargin;

sortBy = lower(string(sortBy));
if sortBy == "none"
    return;
end

nCats = numel(cats);
keys = zeros(nCats, 1);
secondary = zeros(nCats, 1);

for i = 1:nCats
    name = string(cats{i});

    % Special categories go to the end
    if name == "unranged" || name == "unknown" || name == "" || ismissing(name)
        keys(i) = inf;
        secondary(i) = inf;
        continue;
    end

    switch sortBy
        case "atomic"
            [z, iso] = atomicNumberFromName(name);
            keys(i) = z;
            secondary(i) = iso;

        case "weight"
            keys(i) = ionWeightFromName(name);
            secondary(i) = 0;

        case "name"
            % Use ASCII values for alphabetical sort
            keys(i) = i;  % Will be replaced by sortrows on name
            secondary(i) = 0;

        otherwise
            keys(i) = i;
            secondary(i) = 0;
    end
end

if sortBy == "name"
    % Alphabetical sort, keeping special categories at end
    isSpecial = keys == inf;
    [~, order] = sort(string(cats));
    % Move special categories to end
    specialIdx = order(isSpecial(order));
    normalIdx = order(~isSpecial(order));
    order = [normalIdx; specialIdx];
else
    [~, order] = sortrows([keys, secondary]);
end

catsOut = cats(order);

% Reorder all additional data arrays
for k = 1:numel(varargin)
    data = varargin{k};
    if numel(data) == nCats
        if isrow(data)
            varargout{k} = data(order);
        else
            varargout{k} = data(order);
        end
    else
        varargout{k} = data;  % Don't reorder if size doesn't match
    end
end

end

%% Helper Functions

function [z, iso] = atomicNumberFromName(name)
% Extract atomic number and isotope from ion/atom name

    z = inf;
    iso = inf;
    try
        [ionTable, ~] = ionConvertName(char(name));
        if istable(ionTable) && height(ionTable) > 0
            zVals = zeros(height(ionTable), 1);
            isoVals = ionTable.isotope;
            for k = 1:height(ionTable)
                zVals(k) = symbolConvertAtomicNumber(char(ionTable.element(k)));
            end
            z = min(zVals);
            if any(~isnan(isoVals))
                iso = min(isoVals(~isnan(isoVals)));
            else
                iso = inf;
            end
        end
    catch
        z = inf;
        iso = inf;
    end
end

function w = ionWeightFromName(name)
% Calculate mass-to-charge ratio from ion name

    w = inf;
    try
        [ionTable, chargeState] = ionConvertName(char(name));
        if isempty(ionTable) || ~istable(ionTable) || height(ionTable) == 0
            return;
        end
        total = 0;
        for k = 1:height(ionTable)
            z = symbolConvertAtomicNumber(char(ionTable.element(k)));
            if isnan(z)
                continue;
            end
            iso = ionTable.isotope(k);
            total = total + isotopeWeightForElement(z, iso);
        end
        if ~isnan(chargeState) && chargeState ~= 0
            total = total / abs(chargeState);
        end
        w = total;
    catch
        w = inf;
    end
end

function w = isotopeWeightForElement(z, iso)
% Get isotope weight for element

    w = NaN;
    try
        weights = isotopicWeight(z);
        if isempty(weights)
            return;
        end
        if isnan(iso)
            % Use abundance-weighted average
            abund = weights(:, 3);
            w = sum(weights(:, 2) .* abund) / sum(abund);
        else
            idx = weights(:, 1) == iso;
            if any(idx)
                w = weights(find(idx, 1, 'first'), 2);
            else
                abund = weights(:, 3);
                w = sum(weights(:, 2) .* abund) / sum(abund);
            end
        end
    catch
        w = NaN;
    end
end
