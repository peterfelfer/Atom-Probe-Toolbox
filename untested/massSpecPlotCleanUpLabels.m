function [labelHandles, leaderHandles, info] = massSpecPlotCleanUpLabels(spec, varargin)
% massSpecPlotCleanUpLabels rearranges mass spectrum labels to avoid overlap.
%
% [labelHandles, leaderHandles, info] = massSpecPlotCleanUpLabels(spec)
% [labelHandles, leaderHandles, info] = massSpecPlotCleanUpLabels(spec, name, value)
%
% INPUT
% spec:    mass spectrum handle (figure, axes, or plot handle)
% name-value options:
%   'labelOffsetFactor' base multiplier for log axis (default: 1.4)
%   'labelOffset'       base offset (fraction of y-range) for linear axis (default: 0.05)
%   'tierSpacingFactor' tier multiplier for log axis (default: 1.25)
%   'tierSpacing'       tier spacing (fraction of y-range) for linear axis (default: 0.06)
%   'underlineScale'    divide y by this for log underline (default: 1.05)
%   'underlineOffset'   underline offset (fraction of y-range) for linear axis (default: 0.02)
%   'stubFactor'        vertical stub multiplier for log axis (default: 1.15)
%   'stubLength'        vertical stub length (fraction of y-range) for linear axis (default: 0.04)
%   'xPadding'          horizontal padding (fraction of x-range) (default: 0.01)
%   'maxXShift'         maximum x-shift allowed before using next tier (fraction of x-range) (default: 0.2)
%   'maxTiers'          maximum tiers to try (default: 12)
%
% OUTPUT
% labelHandles:  text handles that were rearranged
% leaderHandles: line handles for leader lines
% info:          struct with fields labelsMoved and tiersUsed
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if nargin < 1 || isempty(spec)
    spec = gcf;
end

options = struct( ...
    'labelOffsetFactor', 1.4, ...
    'labelOffset', 0.05, ...
    'tierSpacingFactor', 1.25, ...
    'tierSpacing', 0.06, ...
    'underlineScale', 1.05, ...
    'underlineOffset', 0.02, ...
    'stubFactor', 1.15, ...
    'stubLength', 0.04, ...
    'xPadding', 0.01, ...
    'maxXShift', 0.2, ...
    'maxTiers', 12);

if ~isempty(varargin)
    options = parseOptions(options, varargin{:});
end

ax = resolveAxes(spec);
labels = findLabelHandles(ax);
labelHandles = labels;

if isempty(labels)
    leaderHandles = gobjects(0, 1);
    info = struct('labelsMoved', 0, 'tiersUsed', 0);
    return;
end

delete(findobj(ax, 'Tag', 'massSpecLabelLeader'));

rangeInfo = collectRangeInfo(ax);
isLog = strcmpi(ax.YScale, 'log');

xRange = diff(ax.XLim);
yRange = diff(ax.YLim);
padX = options.xPadding * xRange;

[labelNames, labelAnchors] = mapLabelsToRanges(labels, rangeInfo, options.labelOffsetFactor, isLog);

[~, order] = sort(labelAnchors(:, 1));
labels = labels(order);
labelNames = labelNames(order);
labelAnchors = labelAnchors(order, :);

tiers = cell(1, options.maxTiers);
leaderHandles = gobjects(0, 1);
labelsMoved = 0;
tiersUsed = 0;

for i = 1:numel(labels)
    h = labels(i);
    anchor = labelAnchors(i, :);
    baseY = computeBaseY(ax.YLim(2), yRange, options, isLog);
    xCenter = anchor(1);

    placed = false;
    for tier = 1:options.maxTiers
        yLabel = offsetTier(baseY, tier, yRange, options, isLog);
        set(h, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Rotation', 90, 'Position', [xCenter, yLabel, 0]);
        e = h.Extent;
        labelWidth = max(e(3), padX);
        minCenter = ax.XLim(1) + labelWidth / 2;
        maxCenter = ax.XLim(2) - labelWidth / 2;
        maxShift = options.maxXShift * xRange;
        centerMin = max(xCenter - maxShift, minCenter);
        centerMax = min(xCenter + maxShift, maxCenter);

        [candidateX, ok] = placeInTier(tiers{tier}, xCenter, centerMin, centerMax, labelWidth, padX);
        if ok
            set(h, 'Position', [candidateX, yLabel, 0]);
            e = h.Extent;
            tiers{tier} = [tiers{tier}; ...
                candidateX - labelWidth / 2 - padX, candidateX + labelWidth / 2 + padX];
            placed = true;
            tiersUsed = max(tiersUsed, tier);
            break;
        end
    end

    if ~placed
        tier = options.maxTiers;
        yLabel = offsetTier(baseY, tier, yRange, options, isLog);
        set(h, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'Rotation', 90, 'Position', [xCenter, yLabel, 0]);
        e = h.Extent;
        labelWidth = max(e(3), padX);
        tiers{tier} = [tiers{tier}; ...
            xCenter - labelWidth / 2 - padX, xCenter + labelWidth / 2 + padX];
    end

    labelsMoved = labelsMoved + 1;
    leaderHandles = [leaderHandles; drawLeader(ax, h, anchor, e, options, yRange, xRange, isLog)]; %#ok<AGROW>
end

info = struct('labelsMoved', labelsMoved, 'tiersUsed', tiersUsed);
end

function ax = resolveAxes(spec)
    ax = spec;
    if isgraphics(spec, 'figure')
        axesList = findobj(spec, 'Type', 'axes');
        if ~isempty(axesList)
            ax = axesList(1);
        end
    elseif ~isgraphics(spec, 'axes')
        ax = ancestor(spec, 'axes');
    end
    if isempty(ax) || ~isgraphics(ax, 'axes')
        error('massSpecPlotCleanUpLabels:invalidHandle', ...
            'Could not resolve an axes from the provided handle.');
    end
end

function labels = findLabelHandles(ax)
    labels = findobj(ax, 'Type', 'text');
    keep = false(size(labels));
    for i = 1:numel(labels)
        try
            keep(i) = labels(i).UserData.plotType == "text";
        catch
            keep(i) = false;
        end
    end
    labels = labels(keep);
end

function info = collectRangeInfo(ax)
    plots = ax.Children;
    info = struct('name', string.empty(0, 1), 'xPeak', [], 'yPeak', []);
    for pl = 1:numel(plots)
        try
            type = plots(pl).UserData.plotType;
        catch
            type = "unknown";
        end
        if type ~= "range"
            continue;
        end
        [xPeak, yPeak] = rangePeak(plots(pl));
        name = rangeName(plots(pl));
        info.name(end+1, 1) = string(name);
        info.xPeak(end+1, 1) = xPeak;
        info.yPeak(end+1, 1) = yPeak;
    end
end

function [names, anchors] = mapLabelsToRanges(labels, rangeInfo, offsetFactor, isLog)
    names = strings(numel(labels), 1);
    anchors = zeros(numel(labels), 2);
    for i = 1:numel(labels)
        h = labels(i);
        name = string(h.DisplayName);
        if name == "" && ~isempty(h.String)
            name = string(h.String);
        end
        names(i) = name;
        idx = find(rangeInfo.name == name, 1, 'first');
        if ~isempty(idx)
            anchors(i, :) = [rangeInfo.xPeak(idx), rangeInfo.yPeak(idx)];
        else
            pos = h.Position;
            anchors(i, :) = pos(1:2);
            if isLog
                anchors(i, 2) = anchors(i, 2) / offsetFactor;
            end
        end
    end
end

function baseY = computeBaseY(topY, yRange, options, isLog)
    if isLog
        baseY = topY / options.labelOffsetFactor;
    else
        baseY = topY - options.labelOffset * yRange;
    end
end

function yLabel = offsetTier(baseY, tier, yRange, options, isLog)
    if isLog
        yLabel = baseY / (options.tierSpacingFactor ^ (tier - 1));
    else
        yLabel = baseY - (tier - 1) * options.tierSpacing * yRange;
    end
end

function [candidate, ok] = placeInTier(intervals, center, centerMin, centerMax, labelWidth, padX)
    half = labelWidth / 2 + padX;
    if centerMin > centerMax
        candidate = center;
        ok = false;
        return;
    end
    if isempty(intervals)
        candidate = min(max(center, centerMin), centerMax);
        ok = true;
        return;
    end

    intervals = sortrows(intervals, 1);
    forbidden = [intervals(:, 1) - half, intervals(:, 2) + half];
    forbidden = mergeIntervals(forbidden);

    free = computeFreeIntervals([centerMin, centerMax], forbidden);
    if isempty(free)
        candidate = center;
        ok = false;
        return;
    end

    candidate = pickClosestCenter(free, center);
    ok = true;
end

function merged = mergeIntervals(intervals)
    if isempty(intervals)
        merged = intervals;
        return;
    end
    intervals = sortrows(intervals, 1);
    merged = intervals(1, :);
    for i = 2:size(intervals, 1)
        if intervals(i, 1) <= merged(end, 2)
            merged(end, 2) = max(merged(end, 2), intervals(i, 2));
        else
            merged = [merged; intervals(i, :)]; %#ok<AGROW>
        end
    end
end

function free = computeFreeIntervals(bounds, forbidden)
    free = [];
    current = bounds(1);
    for i = 1:size(forbidden, 1)
        if forbidden(i, 1) > current
            free = [free; current, min(forbidden(i, 1), bounds(2))]; %#ok<AGROW>
        end
        current = max(current, forbidden(i, 2));
        if current >= bounds(2)
            break;
        end
    end
    if current < bounds(2)
        free = [free; current, bounds(2)]; %#ok<AGROW>
    end
end

function center = pickClosestCenter(intervals, target)
    d = abs(intervals - target);
    [~, idx] = min(min(d, [], 2));
    center = min(max(target, intervals(idx, 1)), intervals(idx, 2));
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('massSpecPlotCleanUpLabels:invalidOptions', ...
            'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        value = varargin{k+1};
        switch name
            case "labeloffsetfactor"
                opts.labelOffsetFactor = max(0.1, value);
            case "labeloffset"
                opts.labelOffset = max(0, value);
            case "tierspacingfactor"
                opts.tierSpacingFactor = max(0.1, value);
            case "tierspacing"
                opts.tierSpacing = max(0, value);
            case "underlinescale"
                opts.underlineScale = max(0.1, value);
            case "underlineoffset"
                opts.underlineOffset = max(0, value);
            case "stubfactor"
                opts.stubFactor = max(0.1, value);
            case "stublength"
                opts.stubLength = max(0, value);
            case "xpadding"
                opts.xPadding = max(0, value);
            case "maxxshift"
                opts.maxXShift = max(0, value);
            case "maxtiers"
                opts.maxTiers = max(1, round(value));
            otherwise
                error('massSpecPlotCleanUpLabels:invalidOption', ...
                    'Unknown option "%s".', name);
        end
    end
end

function leaders = drawLeader(ax, h, anchor, extent, options, yRange, xRange, isLog)
    xLeft = extent(1);
    xRight = extent(1) + extent(3);
    yBottom = extent(2);
    yTop = extent(2) + extent(4);
    underlineX = xRight;

    if isLog
        stubLength = max(extent(4), options.stubLength * yRange);
        stubTop = anchor(2) + stubLength;
    else
        stubLength = max(extent(4), options.stubLength * yRange);
        stubTop = anchor(2) + stubLength;
    end

    if stubTop >= yBottom
        stubTop = yBottom - 0.01 * yRange;
    end
    if stubTop <= anchor(2)
        stubTop = anchor(2) + 0.01 * yRange;
    end

    leaders = gobjects(3, 1);
    leaders(1) = line(ax, [underlineX underlineX], [yBottom yTop], ...
        'Color', [0 0 0], 'LineWidth', 0.75, 'Tag', 'massSpecLabelLeader');
    leaders(2) = line(ax, [anchor(1) anchor(1)], [anchor(2) stubTop], ...
        'Color', [0 0 0], 'LineWidth', 0.75, 'Tag', 'massSpecLabelLeader');
    leaders(3) = line(ax, [anchor(1) underlineX], [stubTop yBottom], ...
        'Color', [0 0 0], 'LineWidth', 0.75, 'Tag', 'massSpecLabelLeader');

    for k = 1:3
        leaders(k).UserData.plotType = "labelLeader";
        leaders(k).UserData.label = h;
    end
end

function [xPeak, yPeak] = rangePeak(h)
    x = h.XData;
    y = h.YData;
    if isempty(x) || isempty(y)
        xPeak = h.XData(1);
        yPeak = h.YData(1);
        return;
    end
    [yPeak, idx] = max(y);
    xPeak = x(idx);
end

function name = rangeName(h)
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
