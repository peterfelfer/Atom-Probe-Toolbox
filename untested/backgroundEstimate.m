function [bg, info] = backgroundEstimate(mcCenters, countsNorm, rangeTable, varargin)
% backgroundEstimate estimates background in mass spectrum data.
%
% [bg, info] = backgroundEstimate(mcCenters, countsNorm, rangeTable)
% [bg, info] = backgroundEstimate(..., 'name', value)
%
% This function estimates the background contribution in a mass spectrum
% using one of two methods. The background can then be subtracted from
% peak counts to improve concentration accuracy.
%
% INPUT
%   mcCenters   - mass-to-charge bin centers [Da] (1 x N or N x 1 array)
%   countsNorm  - normalized counts at each bin center (counts/Da/totalIons)
%   rangeTable  - table with mcbegin, mcend columns defining peak regions
%
% NAME-VALUE OPTIONS
%   'method'          - Background fitting method:
%                       'massSpecInvSqrt' (default) - fits B(m/c) = A/sqrt(m/c)
%                       'linearBetweenPeaks' - linear interpolation between gaps
%
%   'minPeakDistance' - Distance from peaks to exclude from fitting [Da]
%                       (default: 0.3). Only used by 'massSpecInvSqrt'.
%
%   'fitLimits'       - Nx2 matrix of [begin, end] m/c ranges to use for fitting
%                       (optional). Example: [5, 20; 40, 60] fits only in
%                       5-20 Da and 40-60 Da regions.
%
% OUTPUT
%   bg   - background values at mcCenters (counts/Da/totalIons)
%          Same size as mcCenters.
%
%   info - struct with fit details:
%          .method      - method used ('massSpecInvSqrt' or 'linearBetweenPeaks')
%          .coefficient - (massSpecInvSqrt) the A parameter in A/sqrt(m/c)
%          .rsquared    - (massSpecInvSqrt) R-squared of fit
%          .fitRegions  - (massSpecInvSqrt) [x, y] points used for fitting
%          .gapFits     - (linearBetweenPeaks) cell array of gap fit structs
%
% METHODS
%   'massSpecInvSqrt' - Fits B(m/c) = A / sqrt(m/c) to unranged regions.
%       Physics basis: A constant background in TOF space transforms to
%       1/sqrt(m/c) in mass spectrum space due to the m/c ~ t^2 relationship.
%       The analytical integral for background counts in a range [a,b] is:
%           bgCounts = 2 * A * (sqrt(b) - sqrt(a)) * totalIons
%
%   'linearBetweenPeaks' - Linear interpolation between gaps.
%       Fits linear models in each gap between ranges and interpolates
%       under peaks. Useful when background shape doesn't follow 1/sqrt.
%
% EXAMPLE
%   % Build histogram from pos data
%   edges = 0.5:0.01:100;
%   counts = histcounts(pos.mc, edges);
%   mcCenters = edges(1:end-1) + 0.005;
%   countsNorm = counts / (0.01 * height(pos));
%
%   % Estimate background
%   [bg, info] = backgroundEstimate(mcCenters, countsNorm, rangeTable, ...
%       'method', 'massSpecInvSqrt', 'fitLimits', [10 80]);
%
%   % Plot
%   figure; semilogy(mcCenters, countsNorm, 'k', mcCenters, bg, 'r--');
%   legend('Spectrum', 'Background');
%
% SEE ALSO
%   posCalculateConcentrationBackgroundRemoved, posCalculateConcentrationDeconvolved
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Parse options
options = struct(...
    'method', 'massSpecInvSqrt', ...
    'minPeakDistance', 0.3, ...
    'fitLimits', []);

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "method"
            options.method = char(value);
        case "minpeakdistance"
            options.minPeakDistance = double(value);
        case "fitlimits"
            options.fitLimits = value;
        otherwise
            error('backgroundEstimate:invalidOption', ...
                'Unknown option "%s".', name);
    end
end

%% Compute background based on method
method = lower(string(options.method));
info = struct();
info.method = char(method);

switch method
    case "massspecinvsqrt"
        [bg, fitInfo] = massSpecInvSqrtBackground(mcCenters, countsNorm, rangeTable, ...
            options.minPeakDistance, options.fitLimits);
        info.coefficient = fitInfo.coefficient;
        info.rsquared = fitInfo.rsquared;
        info.fitRegions = fitInfo.regions;

    case "linearbetweenpeaks"
        [bg, fitInfo] = linearBetweenPeaksBackground(mcCenters, countsNorm, rangeTable);
        info.gapFits = fitInfo.gapFits;

    otherwise
        error('backgroundEstimate:invalidMethod', ...
            'Method must be ''massSpecInvSqrt'' or ''linearBetweenPeaks''.');
end

end

%% ========== Background Methods ==========

function [bg, info] = massSpecInvSqrtBackground(x, y, rangeTable, minPeakDistance, fitLimits)
% Fit B(m/c) = A / sqrt(m/c) to unranged regions

    info = struct();
    info.coefficient = 0;
    info.rsquared = 0;
    info.regions = [];

    % Find unranged regions (valleys between peaks)
    bgMask = true(size(x));
    for i = 1:height(rangeTable)
        inPeak = x >= (rangeTable.mcbegin(i) - minPeakDistance) & ...
                 x <= (rangeTable.mcend(i) + minPeakDistance);
        bgMask = bgMask & ~inPeak;
    end

    % Apply fitLimits if provided
    if nargin >= 5 && ~isempty(fitLimits) && size(fitLimits, 2) >= 2
        inLimits = false(size(x));
        for i = 1:size(fitLimits, 1)
            inLimits = inLimits | (x >= fitLimits(i, 1) & x <= fitLimits(i, 2));
        end
        bgMask = bgMask & inLimits;
    end

    % Exclude very low m/c where 1/sqrt diverges
    bgMask = bgMask & (x > 1);

    xBg = x(bgMask);
    yBg = y(bgMask);

    if numel(xBg) < 3
        warning('backgroundEstimate:insufficientBgData', ...
            'Not enough background data points for fit. Using median.');
        bg = ones(size(x)) * median(y);
        return;
    end

    info.regions = [xBg(:), yBg(:)];

    % Fit model: y = A / sqrt(x)
    % Linear regression: y = A * x^(-0.5)
    xTrans = 1 ./ sqrt(xBg(:));

    % Simple least squares fit through origin: A = sum(x.*y) / sum(x.^2)
    A = sum(xTrans .* yBg(:)) / sum(xTrans.^2);
    if A < 0
        A = 0;
    end

    info.coefficient = A;

    % Calculate R-squared
    yPred = A ./ sqrt(xBg(:));
    ssRes = sum((yBg(:) - yPred).^2);
    ssTot = sum((yBg(:) - mean(yBg)).^2);
    if ssTot > 0
        info.rsquared = 1 - ssRes / ssTot;
    end

    % Generate background for all x values
    bg = A ./ sqrt(max(x, eps));

    if info.rsquared < 0.3
        warning('backgroundEstimate:poorFit', ...
            'Background fit has low R^2 (%.2f). Results may be unreliable.', info.rsquared);
    end
end

function [bg, info] = linearBetweenPeaksBackground(x, y, rangeTable)
% Linear interpolation between peaks

    rangeTable = sortrows(rangeTable, 'mcbegin', 'ascend');
    nRanges = height(rangeTable);
    bg = NaN(size(y));

    % Pre-compute linear fits for gaps
    gapFits = precomputeGapFits(x, y, rangeTable);

    for i = 1:nRanges
        mcbegin = rangeTable.mcbegin(i);
        mcend = rangeTable.mcend(i);
        inRange = x >= mcbegin & x <= mcend;
        if ~any(inRange)
            continue;
        end

        yBegin = evaluateGapFitAtBoundary(gapFits, i, mcbegin, 'left', x, y, rangeTable);
        yEnd = evaluateGapFitAtBoundary(gapFits, i, mcend, 'right', x, y, rangeTable);

        bg(inRange) = linspace(yBegin, yEnd, sum(inRange));
    end

    % Fill regions outside ranges
    bg = fillmissing(bg, 'linear');
    bg = fillmissing(bg, 'nearest');

    info.gapFits = gapFits;
end

%% ========== Helper Functions ==========

function gapFits = precomputeGapFits(x, y, rangeTable)
% Compute linear fits for each gap between adjacent ranges

    nRanges = height(rangeTable);
    gapFits = cell(nRanges + 1, 1);

    % Gap before first range
    gapMask = x < rangeTable.mcbegin(1);
    gapFits{1} = fitGapLinear(x, y, gapMask, rangeTable.mcbegin(1));

    % Gaps between adjacent ranges
    for i = 1:nRanges - 1
        gapMask = x > rangeTable.mcend(i) & x < rangeTable.mcbegin(i + 1);
        touchPoint = (rangeTable.mcend(i) + rangeTable.mcbegin(i + 1)) / 2;
        gapFits{i + 1} = fitGapLinear(x, y, gapMask, touchPoint);
    end

    % Gap after last range
    gapMask = x > rangeTable.mcend(nRanges);
    gapFits{nRanges + 1} = fitGapLinear(x, y, gapMask, rangeTable.mcend(nRanges));
end

function fit = fitGapLinear(x, y, gapMask, touchPoint)
% Fit linear model to gap data

    fit = struct('valid', false, 'coeff', [0 0], 'minVal', NaN, 'touchPoint', touchPoint, 'touchVal', NaN);

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
        fit.coeff = polyfit(xGap(:), yGap(:), 1);
        fit.valid = true;
    end
end

function yVal = evaluateGapFitAtBoundary(gapFits, rangeIdx, xPos, side, x, y, rangeTable)
% Evaluate gap fit at range boundary

    nRanges = height(rangeTable);

    if strcmp(side, 'left')
        primaryFitIdx = rangeIdx;
        altFitIdx = rangeIdx + 1;
    else
        primaryFitIdx = rangeIdx + 1;
        altFitIdx = rangeIdx;
    end

    if (rangeIdx == 1) && strcmp(side, 'left')
        temp = primaryFitIdx;
        primaryFitIdx = altFitIdx;
        altFitIdx = temp;
    end

    yVal = tryEvaluateFit(gapFits{primaryFitIdx}, xPos);
    if ~isnan(yVal)
        return;
    end

    yVal = tryEvaluateFit(gapFits{altFitIdx}, xPos);
    if ~isnan(yVal)
        return;
    end

    inRange = x >= rangeTable.mcbegin(rangeIdx) & x <= rangeTable.mcend(rangeIdx);
    if any(inRange)
        yVal = min(y(inRange));
    else
        yVal = 0;
    end
end

function yVal = tryEvaluateFit(fit, xPos)
% Try to evaluate fit at position

    yVal = NaN;
    if fit.valid
        yVal = polyval(fit.coeff, xPos);
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
        yVal = fit.touchVal;
    end
end
