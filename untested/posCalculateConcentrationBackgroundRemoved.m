function [conc, info] = posCalculateConcentrationBackgroundRemoved(pos, detEff, excludeList, volumeName, varargin)
% posCalculateConcentrationBackgroundRemoved calculates concentrations with background correction.
%
% conc = posCalculateConcentrationBackgroundRemoved(pos, detEff, excludeList, volumeName)
% [conc, info] = posCalculateConcentrationBackgroundRemoved(pos, detEff, excludeList, volumeName, name, value)
%
% INPUT
% pos:          pos table with mc and ion/atom information
% detEff:       detector efficiency (fraction or percent)
% excludeList:  cell array of species to exclude (optional)
% volumeName:   volume name (optional)
%
% name-value options:
%   'method'        'als' | 'minBetweenPeaks' | 'linearBetweenPeaks' | 'tofConstant'
%   'rangeTable'    range table from rangesExtractFromMassSpec (optional)
%   'massSpec'      mass spectrum handle (area plot/axes/figure) (optional)
%   'bin'           histogram bin width if no massSpec handle (default: 0.01)
%   'plotBackground' true/false (default: false)
%   'alsLambda'     smoothness for ALS (default: 1e6)
%   'alsP'          asymmetry for ALS (default: 0.001)
%   'alsIter'       ALS iterations (default: 10)
%   'alsUseUnranged' use only bins outside ranges for ALS (default: true)
%   'tofBackground' constant TOF background level (counts per time bin)
%
% OUTPUT
% conc:   concentration table (same format as posCalculateConcentrationSimple)
% info:   struct with background details
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if detEff > 1
    detEff = detEff / 100;
end

if ~exist('volumeName', 'var') || isempty(volumeName)
    volumeName = inputname(1);
end

if ~exist('excludeList', 'var')
    excludeList = {};
end

options = struct( ...
    'method', 'als', ...
    'rangeTable', [], ...
    'massSpec', [], ...
    'bin', 0.01, ...
    'plotBackground', false, ...
    'alsLambda', 1e6, ...
    'alsP', 0.001, ...
    'alsIter', 10, ...
    'alsUseUnranged', true, ...
    'tofBackground', NaN);

if ~isempty(varargin)
    options = parseOptions(options, varargin{:});
end

% Determine atomic/ionic mode
if any(ismember(pos.Properties.VariableNames, 'atom'))
    type = 'atomic';
    columnType = 'atom';
    atoms = pos.atom;
else
    type = 'ionic';
    columnType = 'ion';
    atoms = pos.ion;
end
atoms(isundefined(atoms)) = 'unranged';

% Build mass spectrum if needed for background
[x, y, ax, specHandle] = resolveMassSpec(pos, options);

% Compute background
bg = computeBackground(x, y, options, specHandle);

% Plot background if requested
bgPlotHandle = gobjects(0, 1);
if options.plotBackground && ~isempty(ax) && isgraphics(ax, 'axes')
    holdState = ishold(ax);
    hold(ax, 'on');
    yPlot = bg;
    if strcmpi(ax.YScale, 'log')
        minY = ax.YLim(1);
        if minY <= 0
            minY = min(yPlot(yPlot > 0));
        end
        yPlot(yPlot < minY) = minY;
    end
    bgPlotHandle = plot(ax, x, yPlot, 'LineWidth', 1.2, 'Color', [1 0 0], 'LineStyle', '--');
    bgPlotHandle.DisplayName = 'background';
    bgPlotHandle.UserData.plotType = "backgroundEstimate";
    if ~holdState
        hold(ax, 'off');
    end
end

% Compute background counts per range if possible
bgCountsByIon = zeros(numel(categories(atoms)), 1);
rangeTable = options.rangeTable;
if isEmptyRangeTable(rangeTable)
    rangeTable = [];
end
if isempty(rangeTable) && ~isempty(specHandle)
    rangeTable = rangesExtractFromMassSpec(specHandle);
end
if isempty(rangeTable) && any(strcmpi(options.method, {'minBetweenPeaks', 'linearBetweenPeaks'}))
    error('posCalculateConcentrationBackgroundRemoved:missingRanges', ...
        'rangeTable or massSpec handle required for %s method.', options.method);
end

if istable(rangeTable) && ~isempty(rangeTable)
    bgCountsByIon = accumulateBackgroundByIon(rangeTable, x, bg, atoms);
end

% Compute corrected counts
cats = categories(atoms);
counts = countcats(atoms);
if iscolumn(counts)
    counts = counts';
end
bgCounts = bgCountsByIon';
if numel(bgCounts) < numel(counts)
    bgCounts(end+1:numel(counts)) = 0;
end
corrCounts = counts - bgCounts;
corrCounts(corrCounts < 0) = 0;

% Exclusion handling
isExcluded = ismember(cats, excludeList);
isExcluded = isExcluded';

% Concentration + variance
concCounts = corrCounts;
concFrac = concCounts ./ sum(concCounts(~isExcluded));
concFrac(isExcluded) = 0;
variance = concFrac .* (1 - concFrac) ./ max(concCounts, 1) * (1 - detEff);

countsOut = [concCounts; concFrac; variance];

% Output table
conc = array2table(countsOut, 'VariableNames', cats');
conc.Properties.VariableDescriptions = repmat({columnType}, size(cats'));
volumeCat = categorical(repmat(string(volumeName), 3, 1));
conc = [table(volumeCat, 'VariableNames', {'volume'}), ...
    table([0; 0; 0], 'VariableNames', {'distance'}), ...
    table(categorical({type; type; type}), 'VariableNames', {'type'}), ...
    table(categorical({'count'; 'concentration'; 'variance'}), 'VariableNames', {'format'}), ...
    conc];

info = struct();
info.method = options.method;
info.massSpecX = x;
info.massSpecY = y;
info.background = bg;
info.backgroundPlot = bgPlotHandle;
info.rangeTable = rangeTable;
info.backgroundCounts = bgCounts;
info.correctedCounts = corrCounts;
end

function [x, y, ax, specHandle] = resolveMassSpec(pos, options)
    ax = [];
    specHandle = [];
    if ~isempty(options.massSpec)
        specHandle = resolveSpecHandle(options.massSpec);
        if isgraphics(specHandle)
            ax = ancestor(specHandle, 'axes');
            x = specHandle.XData;
            y = specHandle.YData;
            return;
        end
    end

    if ~ismember('mc', pos.Properties.VariableNames)
        error('posCalculateConcentrationBackgroundRemoved:missingMc', ...
            'pos must contain mc column if no massSpec handle is given.');
    end
    mc = pos.mc;
    mcmax = max(mc);
    bin = options.bin;
    x = linspace(0, mcmax, max(2, round(mcmax / bin)));
    y = hist(mc, x);
end

function spec = resolveSpecHandle(spec)
    if isgraphics(spec, 'axes')
        plots = findobj(spec, 'Type', 'area');
        spec = pickMassSpecPlot(plots);
    elseif isgraphics(spec, 'figure')
        plots = findobj(spec, 'Type', 'area');
        spec = pickMassSpecPlot(plots);
    elseif isgraphics(spec)
        spec = spec;
    else
        spec = gobjects(0, 1);
    end
end

function spec = pickMassSpecPlot(plots)
    spec = gobjects(0, 1);
    for i = 1:numel(plots)
        try
            if plots(i).UserData.plotType == "massSpectrum"
                spec = plots(i);
                return;
            end
        catch
        end
    end
    if ~isempty(plots)
        spec = plots(1);
    end
end

function bg = computeBackground(x, y, options, specHandle)
    method = lower(string(options.method));
    rangeTable = options.rangeTable;
    if isEmptyRangeTable(rangeTable)
        rangeTable = [];
    end
    if isempty(rangeTable) && ~isempty(specHandle) && any(strcmp(method, ["minbetweenpeaks", "linearbetweenpeaks", "als"]))
        rangeTable = rangesExtractFromMassSpec(specHandle);
    end
    switch method
        case "als"
            if options.alsUseUnranged && istable(rangeTable) && ~isempty(rangeTable)
                mask = computeOutsideRangeMask(x, rangeTable);
            else
                mask = [];
            end
            bg = alsBaseline(y(:), options.alsLambda, options.alsP, options.alsIter, mask);
            bg = bg(:)';
        case "minbetweenpeaks"
            bg = baselineBetweenPeaks(x, y, rangeTable, false);
        case "linearbetweenpeaks"
            if isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0
                error('posCalculateConcentrationBackgroundRemoved:missingRanges', ...
                    'rangeTable required for between-peaks background. Provide rangeTable or massSpec with ranges.');
            end
            bg = baselineBetweenPeaks(x, y, rangeTable, true);
        case "tofconstant"
            bg = tofConstantBackground(x, options.tofBackground);
        otherwise
            error('posCalculateConcentrationBackgroundRemoved:invalidMethod', ...
                'Unknown background method "%s".', options.method);
    end
    bg(bg < 0) = 0;
end

function bg = alsBaseline(y, lambda, p, niter, mask)
    y = y(:);
    L = numel(y);
    D = diff(speye(L), 2);
    w = ones(L, 1);
    if nargin < 5 || isempty(mask)
        mask = ones(L, 1);
    else
        mask = mask(:);
        if numel(mask) ~= L
            mask = ones(L, 1);
        end
    end
    for i = 1:niter
        W = spdiags(w, 0, L, L);
        z = (W + lambda * (D' * D)) \ (w .* y);
        w = (p * (y > z) + (1 - p) * (y < z)) .* mask;
    end
    bg = z;
end

function mask = computeOutsideRangeMask(x, rangeTable)
    mask = true(size(x));
    for i = 1:height(rangeTable)
        mask = mask & ~(x >= rangeTable.mcbegin(i) & x <= rangeTable.mcend(i));
    end
    mask = double(mask);
end

function bg = baselineBetweenPeaks(x, y, rangeTable, useLinear)
    if isempty(rangeTable) || ~istable(rangeTable)
        error('posCalculateConcentrationBackgroundRemoved:missingRanges', ...
            'rangeTable required for between-peaks background.');
    end

    rangeTable = sortrows(rangeTable, 'mcbegin', 'ascend');
    nRanges = height(rangeTable);
    bg = NaN(size(y));

    if useLinear
        % Pre-compute linear fits for ALL gaps between adjacent ranges
        gapFits = precomputeGapFits(x, y, rangeTable);

        for i = 1:nRanges
            mcbegin = rangeTable.mcbegin(i);
            mcend = rangeTable.mcend(i);
            inRange = x >= mcbegin & x <= mcend;
            if ~any(inRange)
                continue;
            end

            % Get background at left boundary from left gap's fit
            yBegin = evaluateGapFitAtBoundary(gapFits, i, mcbegin, 'left', x, y, rangeTable);
            % Get background at right boundary from right gap's fit
            yEnd = evaluateGapFitAtBoundary(gapFits, i, mcend, 'right', x, y, rangeTable);

            bg(inRange) = linspace(yBegin, yEnd, sum(inRange));
        end
    else
        % Original min-between-peaks logic
        for i = 1:nRanges
            mcbegin = rangeTable.mcbegin(i);
            mcend = rangeTable.mcend(i);
            inRange = x >= mcbegin & x <= mcend;
            if ~any(inRange)
                continue;
            end
            leftMin = localGapMin(x, y, rangeTable, i, -1);
            rightMin = localGapMin(x, y, rangeTable, i, 1);
            value = min([leftMin, rightMin], [], 'omitnan');
            if isempty(value) || isnan(value)
                value = min(y(inRange));
            end
            bg(inRange) = value;
        end
    end

    % Fill regions outside ranges to show continuous background
    bg = fillmissing(bg, 'linear');
    bg = fillmissing(bg, 'nearest');
end

function gapFits = precomputeGapFits(x, y, rangeTable)
    % Compute linear fits for each gap between adjacent ranges
    % gapFits{i} contains the fit for the gap BEFORE range i
    % gapFits{nRanges+1} contains the fit for the gap AFTER the last range

    nRanges = height(rangeTable);
    gapFits = cell(nRanges + 1, 1);

    % Gap before first range (from start of spectrum to first range)
    gapMask = x < rangeTable.mcbegin(1);
    touchPoint = rangeTable.mcbegin(1);  % boundary point
    gapFits{1} = fitGapLinear(x, y, gapMask, touchPoint);

    % Gaps between adjacent ranges
    for i = 1:nRanges - 1
        gapMask = x > rangeTable.mcend(i) & x < rangeTable.mcbegin(i + 1);
        % Touch point is where the two ranges meet (or would meet)
        touchPoint = (rangeTable.mcend(i) + rangeTable.mcbegin(i + 1)) / 2;
        gapFits{i + 1} = fitGapLinear(x, y, gapMask, touchPoint);
    end

    % Gap after last range (from last range to end of spectrum)
    gapMask = x > rangeTable.mcend(nRanges);
    touchPoint = rangeTable.mcend(nRanges);  % boundary point
    gapFits{nRanges + 1} = fitGapLinear(x, y, gapMask, touchPoint);
end

function fit = fitGapLinear(x, y, gapMask, touchPoint)
    fit = struct('valid', false, 'coeff', [0 0], 'minVal', NaN, 'touchPoint', touchPoint, 'touchVal', NaN);

    % Try to get spectrum value at touch point (for touching ranges)
    [~, touchIdx] = min(abs(x - touchPoint));
    if ~isempty(touchIdx)
        fit.touchVal = y(touchIdx);
    end

    if ~any(gapMask)
        return;
    end

    xGap = x(gapMask);
    yGap = y(gapMask);
    fit.minVal = min(yGap);

    if numel(xGap) >= 2
        % Fit linear polynomial
        fit.coeff = polyfit(xGap(:), yGap(:), 1);
        fit.valid = true;
    end
end

function yVal = evaluateGapFitAtBoundary(gapFits, rangeIdx, xPos, side, x, y, rangeTable)
    % For 'left' side: use the gap BEFORE this range (gapFits{rangeIdx})
    % For 'right' side: use the gap AFTER this range (gapFits{rangeIdx + 1})
    %
    % Special handling for first peak's left boundary only:
    % - First peak, left boundary: use the inner gap (after first peak) and extrapolate back
    %
    % For last peak: use normal logic (left gap for left, right gap for right)
    % This avoids extrapolating the inner gap beyond its range which can cause artifacts

    nRanges = height(rangeTable);

    if strcmp(side, 'left')
        primaryFitIdx = rangeIdx;
        altFitIdx = rangeIdx + 1;
    else
        primaryFitIdx = rangeIdx + 1;
        altFitIdx = rangeIdx;
    end

    % Only for first peak's left boundary: prefer the inner gap (after first peak)
    isFirstPeakLeft = (rangeIdx == 1) && strcmp(side, 'left');

    if isFirstPeakLeft
        % Swap primary and alt - use the inner gap for first peak's left boundary
        temp = primaryFitIdx;
        primaryFitIdx = altFitIdx;
        altFitIdx = temp;
    end

    % Try primary fit
    yVal = tryEvaluateFit(gapFits{primaryFitIdx}, xPos);

    if ~isnan(yVal)
        return;
    end

    % Try alternate fit
    yVal = tryEvaluateFit(gapFits{altFitIdx}, xPos);

    if ~isnan(yVal)
        return;
    end

    % Last resort: use minimum in the range itself
    inRange = x >= rangeTable.mcbegin(rangeIdx) & x <= rangeTable.mcend(rangeIdx);
    if any(inRange)
        yVal = min(y(inRange));
    else
        yVal = 0;
    end
end

function yVal = tryEvaluateFit(fit, xPos)
    yVal = NaN;
    if fit.valid
        yVal = polyval(fit.coeff, xPos);
        % Ensure non-negative
        if yVal < 0
            if ~isnan(fit.minVal)
                yVal = fit.minVal;
            else
                yVal = NaN;
            end
        end
    elseif ~isnan(fit.minVal)
        yVal = fit.minVal;
    elseif ~isnan(fit.touchVal)
        % Ranges are touching - use spectrum value at touch point
        yVal = fit.touchVal;
    end
end

function val = localGapMin(x, y, rangeTable, idx, direction)
    val = NaN;
    if direction < 0 && idx > 1
        gap = x >= rangeTable.mcend(idx - 1) & x <= rangeTable.mcbegin(idx);
    elseif direction > 0 && idx < height(rangeTable)
        gap = x >= rangeTable.mcend(idx) & x <= rangeTable.mcbegin(idx + 1);
    else
        gap = false(size(x));
    end
    if any(gap)
        val = min(y(gap));
    end
end

function bg = tofConstantBackground(x, tofBackground)
    if isnan(tofBackground)
        error('posCalculateConcentrationBackgroundRemoved:missingTOFBackground', ...
            'tofBackground must be provided for tofConstant method.');
    end
    bg = tofBackground ./ sqrt(max(x, eps));
end

function tf = isEmptyRangeTable(rangeTable)
    tf = isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0;
end

function bgCountsByIon = accumulateBackgroundByIon(rangeTable, x, bg, atoms)
    cats = categories(atoms);
    bgCountsByIon = zeros(numel(cats), 1);
    for i = 1:height(rangeTable)
        inRange = x >= rangeTable.mcbegin(i) & x <= rangeTable.mcend(i);
        if ~any(inRange)
            continue;
        end
        bgCount = sum(bg(inRange));
        ionName = string(rangeTable.rangeName(i));
        idx = find(string(cats) == ionName, 1, 'first');
        if ~isempty(idx)
            bgCountsByIon(idx) = bgCountsByIon(idx) + bgCount;
        end
    end
end

function options = parseOptions(options, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('posCalculateConcentrationBackgroundRemoved:invalidOptions', ...
            'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        value = varargin{k+1};
        switch name
            case "method"
                options.method = char(value);
            case "rangetable"
                options.rangeTable = value;
            case "massspec"
                options.massSpec = value;
            case "bin"
                options.bin = value;
            case "plotbackground"
                options.plotBackground = logical(value);
            case "alslambda"
                options.alsLambda = value;
            case "alsp"
                options.alsP = value;
            case "alsiter"
                options.alsIter = value;
            case "alsuseunranged"
                options.alsUseUnranged = logical(value);
            case "tofbackground"
                options.tofBackground = value;
            otherwise
                error('posCalculateConcentrationBackgroundRemoved:invalidOption', ...
                    'Unknown option "%s".', name);
        end
    end
end
