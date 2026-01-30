function [bg, info] = backgroundEstimate(x, y, method, varargin)
% backgroundEstimate estimates background for mass spectrum data.
%
% [bg, info] = backgroundEstimate(x, y, method)
% [bg, info] = backgroundEstimate(x, y, method, Name, Value)
%
% Estimates background signal in mass spectrum data using various methods.
%
% INPUT
% x:        mass-to-charge values (1xN or Nx1)
% y:        spectrum intensity values (same size as x)
% method:   'als' | 'minBetweenPeaks' | 'linearBetweenPeaks' | 'tofConstant' | 'tofFit' | 'massSpecInvSqrt'
%
% Name-Value Options:
%   'rangeTable'    range table with mcbegin/mcend columns (required for
%                   minBetweenPeaks, linearBetweenPeaks, and tofFit methods)
%   'alsLambda'     smoothness parameter for ALS (default: 1e6)
%   'alsP'          asymmetry parameter for ALS (default: 0.001)
%   'alsIter'       number of ALS iterations (default: 10)
%   'alsUseUnranged' use only bins outside ranges for ALS (default: true)
%   'tofBackground' constant TOF background level (required for tofConstant)
%   'tofBgResult'   result struct from tofSpecBackgroundDetermination (optional)
%   'tofBlockIdx'   block index to use from tofBgResult (default: mean of all)
%   'tofConversion' conversion factor k where m/c = k * t^2 [Da/ns^2] (default: 5e-4)
%   'minPeakDistance' minimum distance from peaks for tofFit [Da] (default: 0.3)
%   'massSpecBgResult' result table from massSpecBackgroundDetermination (optional)
%   'massSpecBlockIdx' block index to use from massSpecBgResult (default: mean of all)
%   'massSpecIonIdx'  ion index to use from massSpecBgResult (nearest)
%
% OUTPUT
% bg:       background estimate (same size as y)
% info:     struct with method-specific details
%
% METHODS
% 'als' - Asymmetric Least Squares baseline estimation
%         Iteratively fits a smooth baseline that stays below the signal.
%         Good for general background estimation.
%
% 'minBetweenPeaks' - Minimum value between adjacent ranges
%         Uses the minimum spectrum value in gaps between ranges.
%         Simple and robust but may underestimate.
%
% 'linearBetweenPeaks' - Linear interpolation between gaps
%         Fits linear models in each gap and interpolates under peaks.
%         More accurate for sloping backgrounds.
%
% 'tofConstant' - Time-of-flight constant background
%         Models background as constant in TOF space (1/sqrt(m/c)).
%         Physics-based model for detector noise.
%         Requires 'tofBackground' parameter.
%
% 'tofFit' - Fit 1/sqrt(m/c) model to valleys
%         Fits B(m/c) = A / sqrt(m/c) + C to background regions.
%         Automatically detects valleys between peaks using rangeTable.
%         Physics-based model that accounts for TOF-to-mass conversion.
%
% 'massSpecInvSqrt' - Use invSqrt coefficient from massSpecBackgroundDetermination
%         Uses B(m/c) = A / sqrt(m/c) with A taken from the provided table.
%
% EXAMPLE
% % Using ALS method
% [bg, info] = backgroundEstimate(x, y, 'als', 'alsLambda', 1e5);
%
% % Using linear between peaks
% [bg, info] = backgroundEstimate(x, y, 'linearBetweenPeaks', 'rangeTable', ranges);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

% Parse options
options = struct( ...
    'rangeTable', [], ...
    'alsLambda', 1e6, ...
    'alsP', 0.001, ...
    'alsIter', 10, ...
    'alsUseUnranged', true, ...
    'tofBackground', NaN, ...
    'tofBgResult', [], ...
    'tofBlockIdx', [], ...
    'tofConversion', 5e-4, ...
    'minPeakDistance', 0.3, ...
    'massSpecBgResult', [], ...
    'massSpecBlockIdx', [], ...
    'massSpecIonIdx', []);

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "rangetable"
            options.rangeTable = value;
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
        case "tofbgresult"
            options.tofBgResult = value;
        case "tofblockidx"
            options.tofBlockIdx = value;
        case "tofconversion"
            options.tofConversion = value;
        case "minpeakdistance"
            options.minPeakDistance = value;
        case "massspecbgresult"
            options.massSpecBgResult = value;
        case "massspecblockidx"
            options.massSpecBlockIdx = value;
        case "massspecionidx"
            options.massSpecIonIdx = value;
    end
end

method = lower(string(method));
if method == "invsqrt" || method == "massspec"
    method = "massspecinvsqrt";
end
rangeTable = options.rangeTable;

if isEmptyRangeTable(rangeTable)
    rangeTable = [];
end

info = struct();
info.method = char(method);

switch method
    case "als"
        if options.alsUseUnranged && ~isempty(rangeTable)
            mask = computeOutsideRangeMask(x, rangeTable);
        else
            mask = [];
        end
        bg = alsBaseline(y(:), options.alsLambda, options.alsP, options.alsIter, mask);
        bg = bg(:)';
        info.lambda = options.alsLambda;
        info.p = options.alsP;
        info.iterations = options.alsIter;

    case "minbetweenpeaks"
        if isempty(rangeTable)
            error('backgroundEstimate:missingRanges', ...
                'rangeTable required for minBetweenPeaks method.');
        end
        bg = baselineBetweenPeaks(x, y, rangeTable, false);

    case "linearbetweenpeaks"
        if isempty(rangeTable)
            error('backgroundEstimate:missingRanges', ...
                'rangeTable required for linearBetweenPeaks method.');
        end
        [bg, gapFits] = baselineBetweenPeaks(x, y, rangeTable, true);
        info.gapFits = gapFits;

    case "tofconstant"
        k = options.tofConversion;
        % Check if tofBgResult is provided
        if ~isempty(options.tofBgResult)
            bgResult = options.tofBgResult;
            if ~isempty(options.tofBlockIdx)
                % Use specific block
                blockIdx = options.tofBlockIdx;
                if blockIdx > numel(bgResult.tofRate)
                    error('backgroundEstimate:invalidBlockIdx', ...
                        'tofBlockIdx %d exceeds number of blocks %d.', blockIdx, numel(bgResult.tofRate));
                end
                tofRate = bgResult.tofRate(blockIdx);
            else
                % Use mean rate
                tofRate = bgResult.meanRate;
            end
            % Convert TOF rate to background using correct formula:
            % B_mc = B_tof / (2 * sqrt(k * m/c))
            % tofRate is in counts/ns/normalization, output is counts/Da/normalization
            bg = tofRate ./ (2 * sqrt(k * max(x, eps)));
            info.tofRate = tofRate;
            if isfield(bgResult, 'normalization')
                info.tofNormalization = bgResult.normalization;
                if strcmpi(bgResult.normalization, 'pulse')
                    info.tofUnits = 'counts/ns/pulse';
                else
                    info.tofUnits = 'counts/ns/totalIons';
                end
            else
                info.tofNormalization = 'ion';
                info.tofUnits = 'counts/ns/totalIons';
            end
            info.tofBlockIdx = options.tofBlockIdx;
            info.tofConversion = k;
        elseif ~isnan(options.tofBackground)
            bg = tofConstantBackground(x, options.tofBackground, k);
            info.tofBackground = options.tofBackground;
            info.tofConversion = k;
        else
            error('backgroundEstimate:missingTOFBackground', ...
                'tofBackground or tofBgResult must be provided for tofConstant method.');
        end

    case "toffit"
        if isempty(rangeTable)
            error('backgroundEstimate:missingRanges', ...
                'rangeTable required for tofFit method.');
        end
        [bg, fitInfo] = tofFitBackground(x, y, rangeTable, options.minPeakDistance);
        info.fitCoefficients = fitInfo.coefficients;
        info.fitRsquared = fitInfo.rsquared;
        info.fitRegions = fitInfo.regions;
        info.fitResidual = fitInfo.residual;

    case "massspecinvsqrt"
        bgTable = options.massSpecBgResult;
        if isempty(bgTable) || ~istable(bgTable) || ~all(ismember({'ionIdx','coeff'}, bgTable.Properties.VariableNames))
            error('backgroundEstimate:missingMassSpecBg', ...
                'massSpecBgResult table with ionIdx and coeff is required for massSpecInvSqrt.');
        end
        if ~isempty(options.massSpecBlockIdx)
            idx = options.massSpecBlockIdx;
            if idx > height(bgTable)
                error('backgroundEstimate:invalidBlockIdx', ...
                    'massSpecBlockIdx %d exceeds number of blocks %d.', idx, height(bgTable));
            end
            coeff = bgTable.coeff(idx);
        elseif ~isempty(options.massSpecIonIdx)
            [~, idx] = min(abs(bgTable.ionIdx - options.massSpecIonIdx));
            coeff = bgTable.coeff(idx);
        else
            coeff = mean(bgTable.coeff, 'omitnan');
        end
        bg = coeff ./ sqrt(max(x, eps));
        info.coeff = coeff;
        info.massSpecBgResult = true;

    otherwise
        error('backgroundEstimate:invalidMethod', ...
            'Unknown background method "%s". Use als, minBetweenPeaks, linearBetweenPeaks, tofConstant, tofFit, or massSpecInvSqrt.', method);
end

bg(bg < 0) = 0;
info.rangeTable = rangeTable;

end

%% Helper Functions

function tf = isEmptyRangeTable(rangeTable)
    tf = isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0;
end

function mask = computeOutsideRangeMask(x, rangeTable)
    mask = true(size(x));
    for i = 1:height(rangeTable)
        mask = mask & ~(x >= rangeTable.mcbegin(i) & x <= rangeTable.mcend(i));
    end
    mask = double(mask);
end

function bg = alsBaseline(y, lambda, p, niter, mask)
% Asymmetric Least Squares baseline estimation
% Based on Eilers & Boelens (2005)

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

function [bg, gapFits] = baselineBetweenPeaks(x, y, rangeTable, useLinear)
% Estimate background using values between peaks

    rangeTable = sortrows(rangeTable, 'mcbegin', 'ascend');
    nRanges = height(rangeTable);
    bg = NaN(size(y));
    gapFits = [];

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
        % Min-between-peaks logic
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

    nRanges = height(rangeTable);
    gapFits = cell(nRanges + 1, 1);

    % Gap before first range
    gapMask = x < rangeTable.mcbegin(1);
    touchPoint = rangeTable.mcbegin(1);
    gapFits{1} = fitGapLinear(x, y, gapMask, touchPoint);

    % Gaps between adjacent ranges
    for i = 1:nRanges - 1
        gapMask = x > rangeTable.mcend(i) & x < rangeTable.mcbegin(i + 1);
        touchPoint = (rangeTable.mcend(i) + rangeTable.mcbegin(i + 1)) / 2;
        gapFits{i + 1} = fitGapLinear(x, y, gapMask, touchPoint);
    end

    % Gap after last range
    gapMask = x > rangeTable.mcend(nRanges);
    touchPoint = rangeTable.mcend(nRanges);
    gapFits{nRanges + 1} = fitGapLinear(x, y, gapMask, touchPoint);
end

function fit = fitGapLinear(x, y, gapMask, touchPoint)
% Fit linear model to gap data

    fit = struct('valid', false, 'coeff', [0 0], 'minVal', NaN, 'touchPoint', touchPoint, 'touchVal', NaN);

    % Get spectrum value at touch point (for touching ranges)
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

    % For first peak's left boundary: prefer the inner gap
    isFirstPeakLeft = (rangeIdx == 1) && strcmp(side, 'left');
    if isFirstPeakLeft
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

function val = localGapMin(x, y, rangeTable, idx, direction)
% Find minimum value in gap adjacent to range

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

function bg = tofConstantBackground(x, tofBackground, k)
% TOF-constant background model
% Converts constant TOF background to m/c space using:
%   B_mc = B_tof / (2 * sqrt(k * m/c))
%
% Input:
%   x - mass-to-charge values [Da]
%   tofBackground - background rate in TOF space [counts/ns/ion]
%   k - TOF conversion factor where m/c = k * t^2 [Da/ns^2]

    if nargin < 3 || isempty(k)
        k = 5e-4;  % Default estimate
    end
    bg = tofBackground ./ (2 * sqrt(k * max(x, eps)));
end


function [bg, fitInfo] = tofFitBackground(x, y, rangeTable, minPeakDistance)
% Fit 1/sqrt(m/c) model to background regions between peaks
%
% Model: B(m/c) = A / sqrt(m/c) + C
%
% The physics: constant background in TOF space transforms to 1/sqrt(m/c)
% in mass spectrum space due to the Jacobian of the m/c = k*t^2 transformation.

    fitInfo = struct();
    fitInfo.coefficients = [0, 0];
    fitInfo.rsquared = 0;
    fitInfo.regions = [];
    fitInfo.residual = 0;

    % Sort ranges by mcbegin
    rangeTable = sortrows(rangeTable, 'mcbegin', 'ascend');
    nRanges = height(rangeTable);

    % Find background regions (valleys between peaks)
    bgMask = true(size(x));

    % Exclude regions within minPeakDistance of any range
    for i = 1:nRanges
        inPeak = x >= (rangeTable.mcbegin(i) - minPeakDistance) & ...
                 x <= (rangeTable.mcend(i) + minPeakDistance);
        bgMask = bgMask & ~inPeak;
    end

    % Also exclude very low m/c values where 1/sqrt diverges
    bgMask = bgMask & (x > 1);

    % Get background data points
    xBg = x(bgMask);
    yBg = y(bgMask);

    if numel(xBg) < 3
        warning('backgroundEstimate:insufficientBgData', ...
            'Not enough background data points for tofFit. Using constant background.');
        bg = ones(size(x)) * median(y);
        return;
    end

    % Store regions used for fitting
    fitInfo.regions = [xBg(:), yBg(:)];

    % Fit model: y = A / sqrt(x) + C
    % Linearize: y = A * x^(-0.5) + C
    % Design matrix: [x^(-0.5), 1]

    xTrans = 1 ./ sqrt(xBg(:));
    X = [xTrans, ones(numel(xBg), 1)];
    yFit = yBg(:);

    % Use non-negative least squares to ensure A >= 0
    % (background should decrease with increasing m/c)
    try
        % Try lsqnonneg for non-negative A
        coeffs = lsqnonneg(X, yFit);
    catch
        % Fallback to regular least squares
        coeffs = X \ yFit;
        if coeffs(1) < 0
            coeffs(1) = 0;
        end
    end

    A = coeffs(1);
    C = coeffs(2);
    if C < 0
        C = 0;
    end

    fitInfo.coefficients = [A, C];

    % Calculate R-squared
    yPred = A ./ sqrt(xBg(:)) + C;
    ssRes = sum((yBg(:) - yPred).^2);
    ssTot = sum((yBg(:) - mean(yBg)).^2);
    if ssTot > 0
        fitInfo.rsquared = 1 - ssRes / ssTot;
    else
        fitInfo.rsquared = 0;
    end
    fitInfo.residual = sqrt(ssRes / numel(yBg));

    % Generate background for all x values
    bg = A ./ sqrt(max(x, eps)) + C;
    bg = reshape(bg, size(x));

    % Warn if fit is poor
    if fitInfo.rsquared < 0.5
        warning('backgroundEstimate:poorTofFit', ...
            'TOF background fit has low R^2 (%.2f). Consider using a different method.', ...
            fitInfo.rsquared);
    end
end
