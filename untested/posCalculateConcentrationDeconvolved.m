function [conc, info] = posCalculateConcentrationDeconvolved(pos, detEff, excludeList, volumeName, varargin)
% posCalculateConcentrationDeconvolved estimates concentrations using overlap deconvolution.
%
% conc = posCalculateConcentrationDeconvolved(pos, detEff, excludeList, volumeName)
% [conc, info] = posCalculateConcentrationDeconvolved(pos, detEff, excludeList, volumeName, name, value)
%
%--------------------------------------------------------------------------
% DESCRIPTION
%--------------------------------------------------------------------------
% This function addresses a fundamental challenge in atom probe tomography:
% overlapping mass-to-charge peaks from different ionic species. When multiple
% ions have isotopes that fall within the same mass range, simple peak
% integration cannot distinguish their contributions. This function uses
% non-negative least squares (NNLS) fitting to deconvolve these overlapping
% peaks based on known natural isotope abundance patterns.
%
%--------------------------------------------------------------------------
% MATHEMATICAL BACKGROUND
%--------------------------------------------------------------------------
% The deconvolution problem is formulated as a linear system:
%
%   c = F * a + residual
%
% where:
%   c (m x 1) = observed counts in each of m mass ranges
%   F (m x n) = fraction matrix, where F(r,i) is the expected fraction of
%               ion i's signal that falls into range r, computed from natural
%               isotope abundances
%   a (n x 1) = unknown true counts for each of n ionic species
%
% The solution is found by minimizing ||F*a - c||^2 subject to a >= 0
% (non-negativity constraint, since counts cannot be negative).
%
% UNCERTAINTY QUANTIFICATION:
% The variance of the NNLS solution is estimated from the residuals:
%   sigma^2 = sum(residuals.^2) / (m - k)
% where k is the number of non-zero components in the solution.
% The covariance matrix is approximated as:
%   Cov(a) = sigma^2 * pinv(F' * F)
% This propagates through mode aggregation (ionic->isotopic->atomic) by
% summing variances of independent variables.
%
%--------------------------------------------------------------------------
% ALGORITHM STEPS
%--------------------------------------------------------------------------
%   1. BUILD FRACTION MATRIX: For each candidate ion, compute the isotope
%      distribution using natural abundances. Map each isotope's mass/charge
%      to the appropriate range, accumulating abundance fractions in F.
%
%   2. IDENTIFY CONNECTED COMPONENTS: Ranges and ions form a bipartite graph.
%      Connected components are groups of ranges/ions that share overlap.
%      Solving each component independently improves numerical stability.
%
%   3. SOLVE NNLS: For each component, solve min ||F_sub * a - c_sub||^2
%      subject to a >= 0 using MATLAB's lsqnonneg. Optionally apply Tikhonov
%      regularization for ill-conditioned problems:
%        min ||F*a - c||^2 + lambda*||a||^2
%
%   4. UNCERTAINTY QUANTIFICATION: Estimate variance of the solution from
%      the residuals and the pseudo-inverse of F'*F (covariance matrix).
%
%   5. GOODNESS-OF-FIT: Compute chi-square statistics and R-squared to
%      assess how well the model fits the observed data.
%
%--------------------------------------------------------------------------
% WHEN TO USE THIS FUNCTION
%--------------------------------------------------------------------------
%   - Fe-Ni systems at 28 Da (56Fe++ overlaps with 58Ni++)
%   - Ti-Cr systems at 24-26 Da
%   - Al-Fe systems (27Al+ overlaps with 54Fe++)
%   - Any system with known isotopic overlaps
%   - When accurate quantification is critical and overlaps are significant
%
%--------------------------------------------------------------------------
% WHEN NOT TO USE THIS FUNCTION
%--------------------------------------------------------------------------
%   - Simple systems without overlaps (use posCalculateConcentrationSimple)
%   - Unknown or non-standard isotope ratios (e.g., enriched samples)
%   - Peaks with significant molecular ion contributions not in the ion table
%   - Very low count situations where Poisson noise dominates
%
%--------------------------------------------------------------------------
% INPUT
%--------------------------------------------------------------------------
% pos:          pos table with mc column (mass-to-charge values)
%               Must contain 'mc' column. Optionally contains 'ion', 'atom',
%               'isotope' columns for raw count comparison.
%
% detEff:       detector efficiency as fraction (0-1) or percent (1-100)
%               Used for variance calculation. Typical values: 0.37-0.80
%               depending on the instrument and operating conditions.
%
% excludeList:  cell array of species names to exclude from concentration
%               calculation (optional, default: {}). Excluded species are
%               still counted but set to 0 in concentration output.
%               Example: {'Ga', 'O'} to exclude gallium and oxygen.
%
% volumeName:   string identifier for this volume/ROI (optional)
%               Used in the output table's 'volume' column.
%               Default: variable name of pos input.
%
%--------------------------------------------------------------------------
% NAME-VALUE OPTIONS
%--------------------------------------------------------------------------
%   'massSpec'       Graphics handle to mass spectrum (area plot, axes, or
%                    figure). Used to extract range table and ion table if
%                    not provided separately. The mass spectrum should have
%                    ranges and ions defined using the standard toolbox
%                    annotation methods.
%
%   'rangeTable'     Table defining mass ranges with columns:
%                    - mcbegin: lower bound of range (Da)
%                    - mcend: upper bound of range (Da)
%                    - rangeName: identifier for the range
%                    - chargeState: charge state of ranged ions
%                    Overrides ranges extracted from massSpec if provided.
%
%   'ionTable'       Table of candidate ions with columns:
%                    - ionName: chemical formula (e.g., 'Fe', 'FeO')
%                    - chargeState: ionic charge state
%                    - ion: (optional) ion categorical for labeling
%                    Extracted from massSpec if not provided.
%
%   'isotopeTable'   Table of isotope data or path to .mat file containing
%                    isotope table. Must have columns for element, mass
%                    number, mass, and natural abundance.
%                    Default: 'isotopeTable_naturalAbundances.mat'
%
%   'mode'           Output aggregation mode:
%                    - 'ionic': report by ionic species (e.g., Fe++, FeO+)
%                    - 'isotopic': report by isotope (e.g., 56Fe, 58Fe)
%                    - 'atomic': report by element (e.g., Fe, Ni, O)
%                    Default: 'atomic' if pos has 'atom' column, else 'ionic'
%
%   'plotFits'       Logical flag to generate diagnostic plot (default: false)
%                    Creates stacked bar chart showing each ion's contribution
%                    to each range, overlaid with observed counts. If background
%                    correction is enabled, background is shown as the bottom layer.
%
%   'plotBackground' Logical flag to plot background on mass spectrum (default: false)
%                    Requires 'massSpec' handle to be provided.
%
%   'colorScheme'    Color scheme for plotting (optional). Can be:
%                    - A colorScheme struct from colorSchemeCreate()
%                    - If not provided, colors are extracted from massSpec
%                    - If neither available, default colors are used
%
%   'regularization' Tikhonov regularization parameter lambda (default: 0)
%                    When lambda > 0, solves:
%                      min ||F*a - c||^2 + lambda*||a||^2
%                    Use small values (0.001-0.1) to stabilize ill-conditioned
%                    problems. Larger values bias toward smaller counts.
%                    A warning is issued when condition number > 100.
%
%   'directCountIons' Cell array of ion names to count directly without
%                    deconvolution (default: {}). Use this for:
%                    - Ions with non-natural isotope ratios (e.g., deuterium
%                      tracers, isotopically enriched samples)
%                    - Ions with no overlaps where deconvolution adds noise
%                    Example: {'H', 'H2'} to directly count hydrogen species
%                    The ion names should match the ionName in ionTable.
%
%   'tracerPeaks'    Cell array of isotopic ion labels that are tracer peaks
%                    (default: {}). For each tracer peak:
%                    1. The peak is excluded from NNLS deconvolution
%                    2. The natural contribution is calculated from other
%                       isotopes of the same element using abundance ratios
%                    3. The tracer excess (measured - natural) is reported
%                       in a separate column with suffix " (tracer)"
%                    Example: {'2H+'} for deuterium tracer experiments
%                    This properly separates natural D from experimental D.
%
%   'backgroundMethod' Background correction method before deconvolution:
%                    - 'none': no background correction (default)
%                    - 'linearBetweenPeaks': linear interpolation between gaps
%                    - 'massSpecInvSqrt': fit B(m/c) = A/sqrt(m/c) to unranged regions
%                    Background is subtracted from range counts before NNLS.
%
%   'minPeakDistance' Minimum distance from peaks for background fitting [Da] (default: 0.3)
%
%   'fitLimits'      Nx2 matrix of [begin, end] m/c ranges for background fitting (optional)
%                    Example: [5, 20; 40, 60] fits only in 5-20 Da and 40-60 Da
%
%   'bin'            Histogram bin width for background estimation (default: 0.01)
%
%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
% conc:  Table with concentration results. Columns are species names,
%        rows are different metrics:
%
%        Row 'countsRaw':     Raw counts from unique detector hits (via pos.ionIdx)
%                             These are direct range counts without deconvolution
%        Row 'count':         Deconvolved counts (after NNLS fitting)
%        Row 'concentration': Atomic/ionic fraction (excluding excludeList)
%        Row 'variance':      Sampling variance of concentration
%                             Formula: p*(1-p)/n * (1-detEff)
%                             where p=concentration, n=counts
%        Row 'countVariance': Variance of deconvolved counts from NNLS
%                             covariance matrix (propagated through mode
%                             aggregation). Represents deconvolution uncertainty.
%        Row 'standardError': Square root of countVariance. Use for error bars.
%        Row 'relativeError': standardError / count (coefficient of variation)
%                             Values > 0.1 (10%) suggest high uncertainty.
%
%        Additional columns:
%        - volume: volume/ROI identifier
%        - distance: placeholder for proxigram distance (0 here)
%        - type: 'ionic', 'isotopic', or 'atomic'
%        - format: row identifier
%
% info:  Struct with detailed diagnostic information:
%
%        .ionCounts       (n x 1) Deconvolved counts for each ion
%        .ionVariance     (n x 1) Variance of ionCounts from NNLS covariance
%        .rangeCounts     (m x 1) Observed counts in each range
%        .fitCounts       (m x 1) Fitted counts (F * ionCounts)
%        .residualCounts  (m x 1) Residuals (rangeCounts - fitCounts)
%        .fractionMatrix  (m x n) The F matrix used for deconvolution
%        .ionLabels       (n x 1) String labels for each ion
%        .components      Struct array of connected components:
%                         .ranges - indices of ranges in component
%                         .ions   - indices of ions in component
%        .conditionNumbers (k x 1) Condition number of F for each component
%                         Values > 100 indicate potential numerical issues
%        .chiSquare       (k x 1) Pearson chi-square statistic per component
%                         Measures goodness-of-fit
%        .pValue          (k x 1) P-value for chi-square test
%                         Low values (< 0.05) suggest poor fit
%        .rSquared        Scalar, overall R-squared (1 - SS_res/SS_tot)
%                         Values close to 1 indicate good fit
%        .isDirectCount   (n x 1) Logical array indicating which ions were
%                         counted directly (bypassed NNLS deconvolution)
%        .backgroundMethod String indicating background method used
%        .backgroundCounts (m x 1) Background counts subtracted per range
%        .detectionLimit  (n x 1) 3-sigma detection limit per ion
%        .isAboveDetectionLimit (n x 1) Logical, true if count > detection limit
%        .rSquaredWeighted Weighted R-squared (by expected counts)
%        .summary         Struct with quick diagnostic statistics:
%                         - totalRangedCounts, totalDeconvolvedCounts
%                         - totalBackgroundCounts, countBalance
%                         - numIllConditioned, meanConditionNumber
%                         - numBelowDetectionLimit
%                         - overallChiSquare, overallDOF
%
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
%   % Basic usage with mass spectrum handle
%   [conc, info] = posCalculateConcentrationDeconvolved(pos, 0.57, {}, ...
%       'myROI', 'massSpec', massSpecHandle);
%
%   % With explicit tables and regularization for ill-conditioned system
%   [conc, info] = posCalculateConcentrationDeconvolved(pos, 0.57, {'Ga'}, ...
%       'ROI1', 'rangeTable', rangeTable, 'ionTable', ionTable, ...
%       'regularization', 0.01, 'plotFits', true);
%
%   % Extract concentration with error bars
%   concRow = conc(conc.format == 'concentration', :);
%   errRow = conc(conc.format == 'standardError', :);
%   elements = conc.Properties.VariableNames(5:end);
%   bar(1:numel(elements), table2array(concRow(1, 5:end)));
%   hold on;
%   errorbar(1:numel(elements), table2array(concRow(1, 5:end)), ...
%       table2array(errRow(1, 5:end)), 'k.');
%   set(gca, 'XTickLabel', elements);
%
%   % Examine fit quality
%   fprintf('R-squared: %.4f\n', info.rSquared);
%   if any(info.conditionNumbers > 100)
%       warning('Some components are ill-conditioned - consider regularization');
%   end
%
%   % Check chi-square p-values for each component
%   for c = 1:numel(info.components)
%       fprintf('Component %d: cond=%.1f, chi2=%.2f, p=%.3f\n', ...
%           c, info.conditionNumbers(c), info.chiSquare(c), info.pValue(c));
%   end
%
%   % Plot residuals to check for systematic errors
%   figure; bar(info.residualCounts);
%   xlabel('Range index'); ylabel('Residual counts');
%   title('Positive residuals may indicate missing ions');
%
%   % Deuterium tracer example using tracerPeaks:
%   % This separates natural D (from H isotope ratios) from experimental D
%   [conc, info] = posCalculateConcentrationDeconvolved(pos, 0.57, {}, 'ROI', ...
%       'rangeTable', rangeTable, 'ionTable', ionTable, ...
%       'tracerPeaks', {'2H+'});
%
%   % The output table will have:
%   % - 'H' column: total H including natural D contribution
%   % - 'H (tracer)' column: excess D above natural abundance (experimental tracer)
%   %
%   % Check tracer info:
%   tracerIdx = find(info.isTracer);
%   for i = 1:numel(tracerIdx)
%       k = tracerIdx(i);
%       fprintf('%s: natural=%.0f, excess=%.0f\n', ...
%           info.ionLabels(k), info.tracerNaturalCounts(k), info.tracerExcessCounts(k));
%   end
%
%   % Alternative: use directCountIons if you want to bypass deconvolution entirely
%   % (e.g., when the reference isotope is also unreliable)
%   [conc, info] = posCalculateConcentrationDeconvolved(pos, 0.57, {}, 'ROI', ...
%       'rangeTable', rangeTable, 'ionTable', ionTable, ...
%       'directCountIons', {'1H+', '2H+'});
%
%--------------------------------------------------------------------------
% SEE ALSO
%--------------------------------------------------------------------------
%   posCalculateConcentrationSimple, ionsCreateIsotopeList, lsqnonneg,
%   rangesExtractFromMassSpec, ionsExtractFromMassSpec
%
%--------------------------------------------------------------------------
% REFERENCES
%--------------------------------------------------------------------------
%   - Lawson, C.L. and Hanson, R.J. "Solving Least Squares Problems."
%     SIAM, 1995. (NNLS algorithm)
%   - Hudson, D., et al. "Optimisation of mass ranging for atom probe
%     microanalysis." Ultramicroscopy (2011).
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
    'massSpec', [], ...
    'rangeTable', [], ...
    'ionTable', [], ...
    'isotopeTable', [], ...
    'mode', '', ...
    'plotFits', false, ...
    'plotBackground', false, ...
    'colorScheme', [], ...
    'regularization', 0, ...
    'directCountIons', {{}}, ...
    'tracerPeaks', {{}}, ...
    'backgroundMethod', 'none', ...
    'minPeakDistance', 0.3, ...
    'fitLimits', [], ...
    'bin', 0.01);

if ~isempty(varargin)
    options = parseOptions(options, varargin{:});
end

% Determine output mode
hasAtomColumn = any(ismember(pos.Properties.VariableNames, 'atom'));
mode = lower(string(options.mode));
if mode == ""
    if hasAtomColumn
        mode = "atomic";
    else
        mode = "ionic";
    end
end

switch mode
    case "ionic"
        columnType = 'ion';
        type = 'ionic';
    case "isotopic"
        columnType = 'isotope';
        type = 'isotopic';
    case "atomic"
        columnType = 'atom';
        type = 'atomic';
    otherwise
        error('posCalculateConcentrationDeconvolved:invalidMode', ...
            'Mode must be ''ionic'', ''isotopic'', or ''atomic''.');
end

if ~ismember('mc', pos.Properties.VariableNames)
    error('posCalculateConcentrationDeconvolved:missingMc', ...
        'pos must contain mc column.');
end

% Resolve spec/ranges/ions
[specHandle, axHandle] = resolveSpecHandle(options.massSpec);
rangeTable = options.rangeTable;
if isEmptyRangeTable(rangeTable) && ~isempty(specHandle)
    rangeTable = rangesExtractFromMassSpec(specHandle);
end

if isEmptyRangeTable(rangeTable)
    error('posCalculateConcentrationDeconvolved:missingRanges', ...
        'rangeTable or massSpec handle required.');
end

ionTable = options.ionTable;
if (isempty(ionTable) || ~istable(ionTable) || height(ionTable) == 0) && ~isempty(specHandle)
    ionTable = ionsExtractFromMassSpec(specHandle);
end

% Counts per range (from pos.mc)
rangeCounts = countsPerRange(pos.mc, rangeTable);
rangeCenters = (rangeTable.mcbegin + rangeTable.mcend) / 2;

% Background correction if requested
backgroundCounts = zeros(size(rangeCounts));
rangeCountsRaw = rangeCounts;  % Keep raw counts for output
bgMethod = lower(string(options.backgroundMethod));
bgInfo = struct();

if bgMethod ~= "none" && bgMethod ~= ""
    % Build internal histogram from pos data
    mcMax = max(pos.mc);
    mcMin = max(0.5, min(pos.mc));
    binWidth = options.bin;
    edges = mcMin:binWidth:mcMax;
    histCounts = histcounts(pos.mc, edges);
    mcCenters = edges(1:end-1) + binWidth/2;

    % Normalize to counts/Da/totalIons for fitting
    totalIons = height(pos);
    countsNorm = histCounts / (binWidth * totalIons);

    % Use shared background estimation function
    [bg, bgFitInfo] = backgroundEstimate(mcCenters, countsNorm, rangeTable, ...
        'method', options.backgroundMethod, ...
        'minPeakDistance', options.minPeakDistance, ...
        'fitLimits', options.fitLimits);

    bgInfo.method = bgFitInfo.method;
    if isfield(bgFitInfo, 'coefficient')
        bgInfo.fitCoefficient = bgFitInfo.coefficient;
    end
    if isfield(bgFitInfo, 'rsquared')
        bgInfo.fitRsquared = bgFitInfo.rsquared;
    end
    if isfield(bgFitInfo, 'fitRegions')
        bgInfo.fitRegions = bgFitInfo.fitRegions;
    end
    if isfield(bgFitInfo, 'gapFits')
        bgInfo.gapFits = bgFitInfo.gapFits;
    end
    bgInfo.fitLimits = options.fitLimits;
    bgInfo.mcCenters = mcCenters;
    bgInfo.background = bg;  % counts/Da/totalIons
    bgInfo.binWidth = binWidth;
    bgInfo.totalIons = totalIons;

    % Calculate background counts per range
    if strcmpi(bgInfo.method, 'massSpecInvSqrt') && isfield(bgInfo, 'fitCoefficient')
        % Use analytical integral: ∫ a/√mc dmc = 2a√mc
        coeff = bgInfo.fitCoefficient;
        for r = 1:height(rangeTable)
            mcBegin = rangeTable.mcbegin(r);
            mcEnd = rangeTable.mcend(r);
            backgroundCounts(r) = 2 * coeff * (sqrt(mcEnd) - sqrt(mcBegin)) * totalIons;
        end
    else
        % Use histogram summing for linearBetweenPeaks
        for r = 1:height(rangeTable)
            inRange = mcCenters >= rangeTable.mcbegin(r) & mcCenters <= rangeTable.mcend(r);
            if any(inRange)
                backgroundCounts(r) = sum(bg(inRange)) * binWidth * totalIons;
            end
        end
    end

    % Subtract background from range counts (ensure non-negative)
    rangeCounts = max(0, rangeCounts - backgroundCounts);
end

% Unranged counts
inAnyRange = atomsInAnyRange(pos.mc, rangeTable);
unrangedCount = sum(~inAnyRange);

% Build candidate ions
isotopeTable = resolveIsotopeTable(options.isotopeTable);
[ionLabels, ionBaseNames, chargeStates] = buildIonLabels(ionTable, rangeTable);

numIons = numel(ionLabels);
numRanges = height(rangeTable);

% No ion candidates -> use raw range counts
if numIons == 0
    ionLabels = string(rangeTable.rangeName);
    ionBaseNames = ionLabels;
    chargeStates = rangeTable.chargeState;
    numIons = numel(ionLabels);
    ionCounts = rangeCounts;
    ionVariance = rangeCounts;  % Poisson variance for raw counts
    F = eye(numRanges, numIons);
    fitCounts = rangeCounts;
    residualCounts = zeros(numRanges, 1);
    components = struct('ranges', {}, 'ions', {});
    condNums = [];
    chiSquare = [];
    pValue = [];
    isDirectCount = false(numIons, 1);  % No direct count ions when no ion table
    isTracer = false(numIons, 1);
    tracerExcessCounts = zeros(numIons, 1);
    tracerNaturalCounts = zeros(numIons, 1);
else
    % Build fraction matrix F (ranges x ions)
    F = zeros(numRanges, numIons);
    ionIsotope = cell(numIons, 1);
    ionAbundance = cell(numIons, 1);
    ionWeight = cell(numIons, 1);

    for k = 1:numIons
        ionName = string(ionBaseNames(k));
        chargeState = chargeStates(k);
        if ~isfinite(chargeState) || chargeState == 0
            continue;
        end
        [isoCombos, abundance, weight] = ionsCreateIsotopeList(ionName, isotopeTable);
        ionIsotope{k} = isoCombos;
        ionAbundance{k} = abundance(:);
        ionWeight{k} = weight(:) ./ abs(chargeState);

        for r = 1:numRanges
            inRange = ionWeight{k} >= rangeTable.mcbegin(r) & ionWeight{k} <= rangeTable.mcend(r);
            if any(inRange)
                F(r, k) = sum(ionAbundance{k}(inRange));
            end
        end
    end

    % Handle direct count ions (bypass deconvolution)
    % These ions are counted directly from their ranges without NNLS
    % Match against full ion labels (e.g., "2H+", "D+", "56Fe++")
    directCountIons = options.directCountIons;
    isDirectCount = false(numIons, 1);
    ionCounts = zeros(numIons, 1);
    ionVariance = zeros(numIons, 1);
    rangeCountsAdj = rangeCounts;  % Adjusted range counts after direct counting

    if ~isempty(directCountIons)
        % Convert to string array for matching
        directCountList = string(directCountIons);

        for k = 1:numIons
            % Match against full ion label (includes isotope and charge, e.g., "2H+")
            % Also try base name for convenience
            ionLabel = string(ionLabels(k));
            ionBase = string(ionBaseNames(k));

            % Check if this ion should be directly counted
            % Match: exact label, label without spaces, or base name
            labelNoSpace = regexprep(ionLabel, '\s+', '');
            isMatch = any(strcmpi(ionLabel, directCountList)) || ...
                      any(strcmpi(labelNoSpace, directCountList)) || ...
                      any(strcmpi(ionBase, directCountList));

            if isMatch
                isDirectCount(k) = true;

                % Find ranges where this ion contributes
                ionRanges = find(F(:, k) > 0);
                if ~isempty(ionRanges)
                    % Direct count: sum range counts weighted by ion's fraction
                    % For single-ion ranges, this gives the full count
                    % For shared ranges, this gives the ion's expected share
                    for r = ionRanges'
                        totalFrac = sum(F(r, :));
                        if totalFrac > 0
                            ionFrac = F(r, k) / totalFrac;
                            ionCounts(k) = ionCounts(k) + rangeCounts(r) * ionFrac;
                            % Subtract from available counts for deconvolution
                            rangeCountsAdj(r) = rangeCountsAdj(r) - rangeCounts(r) * ionFrac;
                        end
                    end
                    % Poisson variance for direct counts
                    ionVariance(k) = ionCounts(k);
                end

                % Remove this ion from F so it's not included in NNLS
                F(:, k) = 0;
            end
        end
    end

    % Handle tracer peaks (separate natural contribution from tracer excess)
    % Tracer peaks are excluded from NNLS, counted directly, then natural
    % contribution is calculated from reference isotopes
    tracerPeaks = options.tracerPeaks;
    isTracer = false(numIons, 1);
    tracerDirectCounts = zeros(numIons, 1);  % Direct counts from range
    tracerNaturalCounts = zeros(numIons, 1); % Calculated natural contribution
    tracerExcessCounts = zeros(numIons, 1);  % Tracer excess (direct - natural)
    tracerRefIdx = zeros(numIons, 1);        % Index of reference isotope ion
    tracerAbundanceRatio = zeros(numIons, 1); % Ratio: tracer_abundance / ref_abundance

    if ~isempty(tracerPeaks)
        tracerList = string(tracerPeaks);

        for k = 1:numIons
            ionLabel = string(ionLabels(k));
            labelNoSpace = regexprep(ionLabel, '\s+', '');

            % Check if this ion is a tracer peak
            isMatch = any(strcmpi(ionLabel, tracerList)) || ...
                      any(strcmpi(labelNoSpace, tracerList));

            if isMatch
                isTracer(k) = true;

                % Count directly from ranges
                ionRanges = find(F(:, k) > 0);
                if ~isempty(ionRanges)
                    for r = ionRanges'
                        totalFrac = sum(F(r, :));
                        if totalFrac > 0
                            ionFrac = F(r, k) / totalFrac;
                            tracerDirectCounts(k) = tracerDirectCounts(k) + rangeCounts(r) * ionFrac;
                            rangeCountsAdj(r) = rangeCountsAdj(r) - rangeCounts(r) * ionFrac;
                        end
                    end
                end

                % Find reference isotope of the same element/molecule
                % Parse the tracer ion to get element and charge
                [tracerElement, tracerIsotope, tracerCharge] = parseIonLabel(ionLabel, ionBaseNames(k), chargeStates(k));
                tracerBaseName = string(ionBaseNames(k));

                % For molecular ions, we need to match the molecule type (e.g., H2 with H2)
                % Extract molecule pattern from baseName (e.g., "H2" from "2H2" or "H2")
                tracerMolecule = regexprep(char(tracerBaseName), '^\d+', '');  % Remove leading isotope number

                % Find other ions of the same element/molecule with same charge (reference isotopes)
                bestRefIdx = 0;
                bestRefAbund = 0;

                for j = 1:numIons
                    if j == k || isTracer(j) || isDirectCount(j)
                        continue;
                    end
                    [refElement, refIsotope, refCharge] = parseIonLabel(ionLabels(j), ionBaseNames(j), chargeStates(j));
                    refBaseName = string(ionBaseNames(j));
                    refMolecule = regexprep(char(refBaseName), '^\d+', '');

                    % Match criteria:
                    % 1. Same element
                    % 2. Same charge state
                    % 3. Same molecule type (H2 matches H2, not H)
                    isSameElement = strcmpi(tracerElement, refElement);
                    isSameCharge = tracerCharge == refCharge;
                    isSameMolecule = strcmpi(tracerMolecule, refMolecule);

                    if isSameElement && isSameCharge && isSameMolecule
                        % Get abundance of this reference isotope
                        refAbund = getIsotopeAbundance(isotopeTable, refElement, refIsotope);

                        % Prefer the most abundant reference isotope
                        if refAbund > bestRefAbund
                            bestRefAbund = refAbund;
                            bestRefIdx = j;
                        end
                    end
                end

                if bestRefIdx > 0
                    tracerRefIdx(k) = bestRefIdx;
                    [refElement, refIsotope, ~] = parseIonLabel(ionLabels(bestRefIdx), ionBaseNames(bestRefIdx), chargeStates(bestRefIdx));

                    % Get natural abundance ratio from isotope table
                    tracerAbund = getIsotopeAbundance(isotopeTable, tracerElement, tracerIsotope);
                    refAbund = getIsotopeAbundance(isotopeTable, refElement, refIsotope);

                    if refAbund > 0
                        tracerAbundanceRatio(k) = tracerAbund / refAbund;
                    end
                end

                % Remove tracer from F so it's not included in NNLS
                F(:, k) = 0;
            end
        end
    end

    % Build connected components between ranges and ions (excluding direct count and tracer ions)
    [components, rangeHasIon] = buildComponents(F);

    fitCounts = zeros(numRanges, 1);
    residualCounts = zeros(numRanges, 1);

    % Diagnostic arrays
    numComponents = numel(components);
    condNums = zeros(numComponents, 1);
    chiSquare = zeros(numComponents, 1);
    pValue = zeros(numComponents, 1);

    lambda = options.regularization;

    for c = 1:numComponents
        rIdx = components(c).ranges;
        iIdx = components(c).ions;
        if isempty(iIdx)
            residualCounts(rIdx) = rangeCountsAdj(rIdx);
            continue;
        end
        Fsub = F(rIdx, iIdx);
        Csub = rangeCountsAdj(rIdx);  % Use adjusted counts (after direct counting)

        % Calculate condition number and warn if ill-conditioned
        condNum = cond(Fsub);
        condNums(c) = condNum;
        if condNum > 100
            warning('posCalculateConcentrationDeconvolved:illConditioned', ...
                'Component %d has condition number %.1f - results may be unreliable.', c, condNum);
        end

        % Handle edge case: single ion per range (no deconvolution needed)
        if numel(iIdx) == 1 && numel(rIdx) == 1
            % Direct assignment - no NNLS needed
            if Fsub > 0
                a = Csub / Fsub;
                % Variance from Poisson statistics
                ionVariance(iIdx) = ionVariance(iIdx) + Csub / (Fsub^2);
            else
                a = 0;
            end
        else
            % Solve NNLS with optional regularization
            if lambda > 0
                Freg = [Fsub; sqrt(lambda) * eye(numel(iIdx))];
                Creg = [Csub; zeros(numel(iIdx), 1)];
                a = lsqnonneg(Freg, Creg);
            else
                a = lsqnonneg(Fsub, Csub);
            end

            % Uncertainty quantification via residual-based covariance
            fit = Fsub * a;
            residuals = Csub - fit;
            dof = max(numel(Csub) - sum(a > 0), 1);  % degrees of freedom
            sigma2 = sum(residuals.^2) / dof;

            % Covariance approximation for active set (non-zero components)
            activeIdx = a > 0;
            if any(activeIdx)
                Factive = Fsub(:, activeIdx);
                FtF = Factive' * Factive;
                if rcond(FtF) > eps
                    covA = sigma2 * pinv(FtF);
                    ionVariance(iIdx(activeIdx)) = ionVariance(iIdx(activeIdx)) + diag(covA);
                end
            end
        end

        ionCounts(iIdx) = ionCounts(iIdx) + a;
        fit = Fsub * a;
        fitCounts(rIdx) = fitCounts(rIdx) + fit;
        residualCounts(rIdx) = residualCounts(rIdx) + (Csub - fit);

        % Chi-square goodness-of-fit test
        expected = fit;
        observed = Csub;
        validBins = expected > 0;
        if any(validBins)
            chiSq = sum((observed(validBins) - expected(validBins)).^2 ./ max(expected(validBins), 1));
            chiSquare(c) = chiSq;
            dofChi = sum(validBins) - sum(a > 0);
            if dofChi > 0
                % Use chi2cdf if Statistics Toolbox available, otherwise approximate
                try
                    pValue(c) = 1 - chi2cdf(chiSq, dofChi);
                catch
                    % Approximate using normal distribution for large dof
                    z = (chiSq - dofChi) / sqrt(2 * dofChi);
                    pValue(c) = 1 - 0.5 * (1 + erf(z / sqrt(2)));
                end
            else
                pValue(c) = NaN;
            end
        end
    end

    % Ranges with no ion candidates keep full residual
    noIonRanges = ~rangeHasIon;
    residualCounts(noIonRanges) = rangeCounts(noIonRanges);

    % Calculate tracer contributions now that we have deconvolved reference counts
    for k = 1:numIons
        if isTracer(k) && tracerRefIdx(k) > 0
            refIdx = tracerRefIdx(k);
            refCount = ionCounts(refIdx);
            refVar = ionVariance(refIdx);

            % Natural contribution = reference count * abundance ratio
            tracerNaturalCounts(k) = refCount * tracerAbundanceRatio(k);

            % Tracer excess = direct count - natural contribution
            tracerExcessCounts(k) = max(0, tracerDirectCounts(k) - tracerNaturalCounts(k));

            % The ion count for the tracer is the natural contribution
            % (tracer excess is reported separately)
            ionCounts(k) = tracerNaturalCounts(k);

            % Propagate variance: Var(natural) = Var(ref) * ratio^2
            % This properly accounts for the uncertainty in the reference count
            ionVariance(k) = refVar * tracerAbundanceRatio(k)^2;
        elseif isTracer(k)
            % No reference found - treat as direct count
            ionCounts(k) = tracerDirectCounts(k);
            ionVariance(k) = tracerDirectCounts(k);  % Poisson variance for direct count
            tracerExcessCounts(k) = 0;  % Can't separate without reference
            warning('posCalculateConcentrationDeconvolved:noTracerReference', ...
                'No reference isotope found for tracer %s. Using direct count.', ionLabels(k));
        end
    end
end

% Build info output struct
info = struct();
info.ionCounts = ionCounts;
info.ionVariance = ionVariance;
info.rangeCounts = rangeCounts;
info.fitCounts = fitCounts;
info.residualCounts = residualCounts;
info.fractionMatrix = F;
info.ionLabels = ionLabels;
info.components = components;
info.conditionNumbers = condNums;
info.chiSquare = chiSquare;
info.pValue = pValue;
info.isDirectCount = isDirectCount;  % Which ions were counted directly (not deconvolved)
info.isTracer = isTracer;            % Which ions are tracer peaks
info.tracerExcessCounts = tracerExcessCounts;  % Tracer excess (above natural)
info.tracerNaturalCounts = tracerNaturalCounts; % Natural contribution to tracer peaks

% Background correction info
info.backgroundMethod = char(bgMethod);
info.backgroundCounts = backgroundCounts;

% Detection limit estimation (3-sigma criterion)
info.detectionLimit = 3 * sqrt(max(ionVariance, 1));
info.isAboveDetectionLimit = ionCounts > info.detectionLimit;

% Calculate overall R-squared
totalSS = sum((rangeCounts - mean(rangeCounts)).^2);
residualSS = sum(residualCounts.^2);
if totalSS > 0
    info.rSquared = 1 - residualSS / totalSS;
else
    info.rSquared = NaN;
end

% Weighted R-squared (by expected counts - more meaningful for count data)
weights = max(fitCounts, 1);
totalSS_w = sum(weights .* (rangeCounts - sum(rangeCounts .* weights) / sum(weights)).^2);
residualSS_w = sum(weights .* residualCounts.^2);
if totalSS_w > 0
    info.rSquaredWeighted = 1 - residualSS_w / totalSS_w;
else
    info.rSquaredWeighted = NaN;
end

% Summary statistics for quick diagnostics
info.summary = struct(...
    'totalRangedCounts', sum(rangeCounts), ...
    'totalDeconvolvedCounts', sum(ionCounts), ...
    'totalBackgroundCounts', sum(backgroundCounts), ...
    'countBalance', sum(ionCounts) / max(sum(rangeCounts), 1), ...  % Should be ~1
    'numIllConditioned', sum(condNums > 100), ...
    'meanConditionNumber', mean(condNums(condNums > 0)), ...
    'numBelowDetectionLimit', sum(~info.isAboveDetectionLimit), ...
    'overallChiSquare', sum(chiSquare), ...
    'overallDOF', max(0, sum(rangeCounts > 0) - sum(ionCounts > 0)));

if options.plotFits
    plotStackedFit(rangeCenters, rangeCountsRaw, F, ionCounts, ionVariance, ionLabels, backgroundCounts, specHandle, axHandle, options.colorScheme);
end

if options.plotBackground && ~isempty(specHandle) && isfield(bgInfo, 'background')
    plotBackgroundOnMassSpec(specHandle, bgInfo.mcCenters, bgInfo.background, bgInfo.binWidth, bgInfo.totalIons);
end

% Aggregate counts for output mode
ionCat = categorical(ionLabels);
countByIon = ionCounts(:)';
varByIon = ionVariance(:)';

% Add unranged (Poisson variance = count for unranged)
ionCat = addcats(ionCat, 'unranged');
ionCat(end+1) = 'unranged';
countByIon(end+1) = unrangedCount;
varByIon(end+1) = unrangedCount;  % Poisson variance for unranged counts

switch mode
    case "ionic"
        outCat = ionCat;
        outCounts = countByIon;
        outVar = varByIon;
    otherwise
        % Convert ionic categories to output mode
        convertedCat = ionConvertMode(ionCat, char(mode));
        [outCat, outCounts] = aggregateCounts(convertedCat, countByIon);
        % Variance aggregation: sum variances for independent variables
        % Must aggregate varByIon using the same converted categories
        [~, outVar] = aggregateCounts(convertedCat, varByIon);
end

% Raw counts from unique detector hits
% Use pos.ionIdx to get unique hits, then apply ranges directly
[rawCat, rawCounts] = rawCountsFromUniqueHits(pos, rangeTable, mode);

% Union categories: deconvolved + raw
% This ensures we see both raw and deconvolved data even if labels differ
cats = categories(outCat);
rawCats = categories(rawCat);
for i = 1:numel(rawCats)
    if ~any(strcmp(cats, rawCats{i}))
        cats{end+1} = rawCats{i}; %#ok<AGROW>
    end
end

countsRaw = zeros(1, numel(cats));
for i = 1:numel(rawCats)
    idx = find(strcmp(cats, rawCats{i}), 1, 'first');
    if ~isempty(idx)
        countsRaw(idx) = rawCounts(i);
    end
end

counts = zeros(1, numel(cats));
countVar = zeros(1, numel(cats));
outNames = categories(outCat);
for i = 1:numel(outNames)
    idx = find(strcmp(cats, outNames{i}), 1, 'first');
    if ~isempty(idx)
        counts(idx) = sum(outCounts(outCat == outNames{i}));
        countVar(idx) = sum(outVar(outCat == outNames{i}));
    end
end

% Add tracer excess columns
% Aggregate by converted category name (combines multiple charge states/molecules in atomic mode)
tracerCatMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
tracerVarMap = containers.Map('KeyType', 'char', 'ValueType', 'double');

if any(isTracer)
    tracerIdx = find(isTracer);
    for t = 1:numel(tracerIdx)
        k = tracerIdx(t);
        % Create tracer label based on output mode
        tracerLabel = string(ionLabels(k));
        switch mode
            case "ionic"
                tracerCatName = char(tracerLabel + " (tracer)");
            case "isotopic"
                % Convert to isotopic label - handles molecular ions
                converted = ionConvertMode(categorical(tracerLabel), 'isotopic');
                % For molecules, this might return multiple isotopes - take first
                convStr = string(converted);
                if strlength(convStr) > 0
                    tracerCatName = char(convStr + " (tracer)");
                else
                    tracerCatName = char(tracerLabel + " (tracer)");
                end
            case "atomic"
                % Convert to atomic label - handles molecular ions
                converted = ionConvertMode(categorical(tracerLabel), 'atomic');
                convStr = string(converted);
                if strlength(convStr) > 0
                    tracerCatName = char(convStr + " (tracer)");
                else
                    tracerCatName = char(tracerLabel + " (tracer)");
                end
        end

        % Aggregate: sum counts for same category (e.g., 2H+ and 2H++ both -> "H (tracer)")
        if isKey(tracerCatMap, tracerCatName)
            tracerCatMap(tracerCatName) = tracerCatMap(tracerCatName) + tracerExcessCounts(k);
            tracerVarMap(tracerCatName) = tracerVarMap(tracerCatName) + tracerExcessCounts(k);
        else
            tracerCatMap(tracerCatName) = tracerExcessCounts(k);
            tracerVarMap(tracerCatName) = tracerExcessCounts(k);  % Poisson variance
        end
    end
end

% Convert maps to arrays
tracerCats = keys(tracerCatMap);
tracerCounts = zeros(1, numel(tracerCats));
tracerVar = zeros(1, numel(tracerCats));
for i = 1:numel(tracerCats)
    tracerCounts(i) = tracerCatMap(tracerCats{i});
    tracerVar(i) = tracerVarMap(tracerCats{i});
end

% Append tracer categories to main category list
cats = [cats, tracerCats];
counts = [counts, tracerCounts];
countVar = [countVar, tracerVar];
countsRaw = [countsRaw, zeros(1, numel(tracerCats))];  % No raw counts for tracer excess

% Sort categories by mass (unranged and tracer categories go at the end)
[cats, counts, countVar, countsRaw] = sortCategoriesByMass(cats, counts, countVar, countsRaw, isotopeTable, mode);

% Exclusion handling
isExcluded = ismember(cats, excludeList);
isExcluded = isExcluded';

concCounts = counts;
concFrac = concCounts ./ sum(concCounts(~isExcluded));
concFrac(isExcluded) = 0;
variance = concFrac .* (1 - concFrac) ./ max(concCounts, 1) * (1 - detEff);

% Deconvolution uncertainty metrics
countVariance = countVar;  % Variance from NNLS covariance, propagated through mode
standardError = sqrt(countVariance);  % Standard error of counts
relativeError = standardError ./ max(concCounts, 1);  % Coefficient of variation
relativeError(concCounts == 0) = 0;  % Avoid division artifacts for zero counts

% Calculate background counts per category
% Distribute total background proportionally to deconvolved counts
totalBackground = sum(backgroundCounts);
bgCounts = zeros(size(counts));
if totalBackground > 0 && sum(concCounts) > 0
    % Distribute proportionally to deconvolved counts
    bgCounts = totalBackground * (concCounts / sum(concCounts));
end
bgCounts(strcmpi(cats, 'unranged')) = 0;  % No background for unranged

% Background fraction (background / raw counts)
bgFraction = bgCounts ./ max(countsRaw, 1);
bgFraction(countsRaw == 0) = 0;

countsOut = [countsRaw; bgCounts; bgFraction; concCounts; concFrac; variance; countVariance; standardError; relativeError];

% Output table
numRows = 9;
conc = array2table(countsOut, 'VariableNames', cats');
conc.Properties.VariableDescriptions = repmat({columnType}, size(cats'));
volumeCat = categorical(repmat(string(volumeName), numRows, 1));

% Type column: countsRaw is always 'ionic' (detector hit counts), others use output mode
typeCol = categorical({...
    'ionic'; ...      % countsRaw - always ionic (raw detector hits per range)
    type; ...         % backgroundCounts
    type; ...         % backgroundFraction
    type; ...         % count
    type; ...         % concentration
    type; ...         % variance
    type; ...         % countVariance
    type; ...         % standardError
    type});           % relativeError

conc = [table(volumeCat, 'VariableNames', {'volume'}), ...
    table(zeros(numRows, 1), 'VariableNames', {'distance'}), ...
    table(typeCol, 'VariableNames', {'type'}), ...
    table(categorical({'countsRaw'; 'backgroundCounts'; 'backgroundFraction'; 'count'; 'concentration'; 'variance'; ...
        'countVariance'; 'standardError'; 'relativeError'}), 'VariableNames', {'format'}), ...
    conc];

end

function options = parseOptions(options, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('posCalculateConcentrationDeconvolved:invalidOptions', ...
            'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        value = varargin{k+1};
        switch name
            case "massspec"
                options.massSpec = value;
            case "rangetable"
                options.rangeTable = value;
            case "iontable"
                options.ionTable = value;
            case "isotopetable"
                options.isotopeTable = value;
            case "mode"
                options.mode = char(value);
            case "plotfits"
                options.plotFits = logical(value);
            case "plotbackground"
                options.plotBackground = logical(value);
            case "colorscheme"
                options.colorScheme = value;
            case "regularization"
                options.regularization = double(value);
            case "directcountions"
                if ischar(value) || isstring(value)
                    options.directCountIons = {char(value)};
                else
                    options.directCountIons = value;
                end
            case "tracerpeaks"
                if ischar(value) || isstring(value)
                    options.tracerPeaks = {char(value)};
                else
                    options.tracerPeaks = value;
                end
            case "backgroundmethod"
                options.backgroundMethod = char(value);
            case "minpeakdistance"
                options.minPeakDistance = double(value);
            case "fitlimits"
                options.fitLimits = value;
            case "bin"
                options.bin = double(value);
            otherwise
                error('posCalculateConcentrationDeconvolved:invalidOption', ...
                    'Unknown option "%s".', name);
        end
    end
end

function [spec, ax] = resolveSpecHandle(spec)
    ax = [];
    if isempty(spec)
        spec = gobjects(0, 1);
        return;
    end
    if isgraphics(spec, 'axes')
        ax = spec;
        plots = findobj(spec, 'Type', 'area');
        spec = pickMassSpecPlot(plots);
        return;
    elseif isgraphics(spec, 'figure')
        plots = findobj(spec, 'Type', 'area');
        spec = pickMassSpecPlot(plots);
    elseif isgraphics(spec)
        % spec is already a graphics handle, keep as-is
    else
        spec = gobjects(0, 1);
    end
    if isgraphics(spec)
        ax = ancestor(spec, 'axes');
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

function tf = isEmptyRangeTable(rangeTable)
    tf = isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0;
end

function rangeCounts = countsPerRange(mc, rangeTable)
    numRanges = height(rangeTable);
    rangeCounts = zeros(numRanges, 1);
    for r = 1:numRanges
        inRange = mc >= rangeTable.mcbegin(r) & mc <= rangeTable.mcend(r);
        rangeCounts(r) = sum(inRange);
    end
end

function mask = atomsInAnyRange(mc, rangeTable)
    mask = false(size(mc));
    for r = 1:height(rangeTable)
        mask = mask | (mc >= rangeTable.mcbegin(r) & mc <= rangeTable.mcend(r));
    end
end

function isotopeTable = resolveIsotopeTable(isotopeTableInput)
    if istable(isotopeTableInput)
        isotopeTable = isotopeTableInput;
        return;
    end

    if isstring(isotopeTableInput) || ischar(isotopeTableInput)
        fileName = char(isotopeTableInput);
    else
        fileName = 'isotopeTable_naturalAbundances.mat';
    end

    if ~exist(fileName, 'file')
        localFile = fullfile(fileparts(mfilename('fullpath')), '..', fileName);
        if exist(localFile, 'file')
            fileName = localFile;
        else
            error('posCalculateConcentrationDeconvolved:missingIsotopeTable', ...
                'Could not find isotope table file: %s', fileName);
        end
    end

    data = load(fileName);
    fields = fieldnames(data);
    isotopeTable = [];
    for k = 1:numel(fields)
        val = data.(fields{k});
        if istable(val)
            isotopeTable = val;
            break;
        end
    end
    if isempty(isotopeTable)
        error('posCalculateConcentrationDeconvolved:invalidIsotopeTable', ...
            'No table found in isotope table file.');
    end
end

function [labels, baseNames, chargeStates] = buildIonLabels(ionTable, rangeTable)
    labels = string.empty(0, 1);
    baseNames = string.empty(0, 1);
    chargeStates = [];

    if isempty(ionTable) || ~istable(ionTable) || height(ionTable) == 0
        return;
    end

    numIons = height(ionTable);
    labels = strings(numIons, 1);
    baseNames = strings(numIons, 1);
    chargeStates = zeros(numIons, 1);

    for k = 1:numIons
        chargeStates(k) = ionTable.chargeState(k);
        baseNames(k) = string(ionTable.ionName(k));
        try
            if ismember('ion', ionTable.Properties.VariableNames)
                ionEntry = ionTable.ion(k);
                if iscell(ionTable.ion)
                    ionEntry = ionTable.ion{k};
                end
                labels(k) = string(ionConvertName(ionEntry, ionTable.chargeState(k)));
            else
                labels(k) = appendCharge(baseNames(k), ionTable.chargeState(k));
            end
        catch
            labels(k) = appendCharge(baseNames(k), ionTable.chargeState(k));
        end
    end

    % remove duplicates while preserving order
    [labels, ia] = unique(labels, 'stable');
    baseNames = baseNames(ia);
    chargeStates = chargeStates(ia);

    if isempty(labels)
        labels = string(rangeTable.rangeName);
        baseNames = labels;
        chargeStates = rangeTable.chargeState;
    end
end

function out = appendCharge(name, chargeState)
    out = string(name);
    if ~isfinite(chargeState) || chargeState == 0
        return;
    end
    if chargeState < 0
        sym = '-';
    else
        sym = '+';
    end
    out = strtrim(out + " " + repmat(sym, 1, abs(chargeState)));
end

function [components, rangeHasIon] = buildComponents(F)
    numRanges = size(F, 1);
    numIons = size(F, 2);
    rangeVisited = false(numRanges, 1);
    ionVisited = false(numIons, 1);
    % Pre-allocate with upper bound (each range as separate component)
    components = struct('ranges', cell(1, numRanges), 'ions', cell(1, numRanges));
    compIdx = 0;

    rangeHasIon = any(F > 0, 2);

    for r = 1:numRanges
        if rangeVisited(r) || ~rangeHasIon(r)
            continue;
        end
        compIdx = compIdx + 1;
        % Pre-allocate arrays with upper bounds
        ranges = zeros(numRanges, 1);
        ions = zeros(numIons, 1);
        rangeCount = 0;
        ionCount = 0;

        queueRanges = zeros(numRanges, 1);
        queueHead = 1;
        queueTail = 1;
        queueRanges(queueTail) = r;
        queueTail = queueTail + 1;
        rangeVisited(r) = true;

        while queueHead < queueTail
            curRange = queueRanges(queueHead);
            queueHead = queueHead + 1;
            rangeCount = rangeCount + 1;
            ranges(rangeCount) = curRange;

            ionIdx = find(F(curRange, :) > 0);
            for i = ionIdx
                if ~ionVisited(i)
                    ionVisited(i) = true;
                    ionCount = ionCount + 1;
                    ions(ionCount) = i;
                    linkedRanges = find(F(:, i) > 0);
                    for lr = linkedRanges'
                        if ~rangeVisited(lr)
                            rangeVisited(lr) = true;
                            queueRanges(queueTail) = lr;
                            queueTail = queueTail + 1;
                        end
                    end
                end
            end
        end

        components(compIdx).ranges = ranges(1:rangeCount);
        components(compIdx).ions = ions(1:ionCount);
    end

    % Trim unused components
    if compIdx < numRanges
        components = components(1:compIdx);
    end
end

function plotStackedFit(rangeCenters, rangeCounts, F, ionCounts, ionVariance, ionLabels, backgroundCounts, specHandle, axHandle, colorScheme)
    % Calculate contributions per ion per range
    contrib = F .* ionCounts';
    numRanges = numel(rangeCenters);
    numIons = numel(ionLabels);

    % Calculate fit and uncertainty per range (ion contributions only, not background)
    fitCounts = sum(contrib, 2);

    % Propagate variance to range level: Var(sum) = sum of variances (weighted by F^2)
    % For each range: variance = sum_i (F(r,i)^2 * ionVariance(i))
    rangeVariance = zeros(numRanges, 1);
    for r = 1:numRanges
        rangeVariance(r) = sum((F(r, :)'.^2) .* ionVariance);
    end
    rangeStdErr = sqrt(rangeVariance);

    % Get ion colors: try massSpec first, then colorScheme, then defaults
    ionColors = getIonColorsFromSpec(specHandle, ionLabels);
    if isempty(ionColors) && ~isempty(colorScheme)
        ionColors = getIonColorsFromScheme(colorScheme, ionLabels);
    end

    % Create figure
    fig = figure('Name', 'Deconvolution fit');
    ax = axes(fig);

    % Build stacked data: background at bottom, then ion contributions
    hasBackground = any(backgroundCounts > 0);
    if hasBackground
        stackedData = [backgroundCounts(:), contrib];
        stackLabels = ["Background"; ionLabels(:)];
    else
        stackedData = contrib;
        stackLabels = ionLabels(:);
    end

    % Plot stacked bars
    bh = bar(ax, rangeCenters, stackedData, 'stacked', 'EdgeColor', 'none');

    % Apply colors: gray for background, then ion colors
    if hasBackground
        bh(1).FaceColor = [0.6 0.6 0.6];  % Gray for background
        if ~isempty(ionColors)
            for i = 1:min(numel(bh)-1, size(ionColors, 1))
                if all(~isnan(ionColors(i, :)))
                    bh(i+1).FaceColor = ionColors(i, :);
                end
            end
        end
    else
        if ~isempty(ionColors)
            for i = 1:min(numel(bh), size(ionColors, 1))
                if all(~isnan(ionColors(i, :)))
                    bh(i).FaceColor = ionColors(i, :);
                end
            end
        end
    end

    hold(ax, 'on');

    % Plot observed counts (raw, before background subtraction)
    plot(ax, rangeCenters, rangeCounts, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'Observed');

    % Plot fit with uncertainty band (total = background + ion contributions)
    % Sort by range center for proper line plotting
    [sortedCenters, sortIdx] = sort(rangeCenters);
    totalFit = fitCounts + backgroundCounts(:);
    sortedFit = totalFit(sortIdx);
    sortedStdErr = rangeStdErr(sortIdx);

    % Uncertainty band (±1 sigma)
    fill(ax, [sortedCenters; flipud(sortedCenters)], ...
         [sortedFit + sortedStdErr; flipud(sortedFit - sortedStdErr)], ...
         [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'none', ...
         'DisplayName', 'Fit ±1σ');

    % Fit line
    plot(ax, sortedCenters, sortedFit, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Fit');

    % Residuals
    residuals = rangeCounts - totalFit;

    hold(ax, 'off');

    xlabel(ax, 'mass-to-charge [Da]');
    ylabel(ax, 'counts [#]');
    if hasBackground
        title(ax, 'Deconvolution: observed vs. fit (background + ions)');
    else
        title(ax, 'Deconvolution: observed vs. fit');
    end

    % Legend: show background (if present), ions (limited), plus observed/fit
    numLabels = numel(stackLabels);
    if numLabels <= 20
        legendEntries = [stackLabels; "Observed"; "Fit ±1σ"; "Fit"];
        legend(ax, legendEntries, 'Interpreter', 'none', 'Location', 'best');
    else
        if hasBackground
            legend(ax, {'Background', 'Observed', 'Fit ±1σ', 'Fit'}, 'Location', 'best');
        else
            legend(ax, {'Observed', 'Fit ±1σ', 'Fit'}, 'Location', 'best');
        end
    end

    if ~isempty(axHandle) && isgraphics(axHandle, 'axes')
        ax.XLim = axHandle.XLim;
    end

    % Add second subplot for residuals
    ax2 = axes(fig, 'Position', [0.13 0.11 0.775 0.15]);
    stem(ax2, rangeCenters, residuals, 'k', 'Marker', 'none', 'LineWidth', 1);
    hold(ax2, 'on');
    yline(ax2, 0, 'k--');
    % Show ±1 sigma bounds
    plot(ax2, sortedCenters, sortedStdErr, 'r--', 'LineWidth', 0.5);
    plot(ax2, sortedCenters, -sortedStdErr, 'r--', 'LineWidth', 0.5);
    hold(ax2, 'off');
    xlabel(ax2, 'mass-to-charge [Da]');
    ylabel(ax2, 'residual [#]');
    if ~isempty(axHandle) && isgraphics(axHandle, 'axes')
        ax2.XLim = axHandle.XLim;
    end

    % Adjust main axes position to make room for residuals
    ax.Position = [0.13 0.35 0.775 0.55];
end

function ionColors = getIonColorsFromSpec(specHandle, ionLabels)
% Extract ion colors from mass spectrum handle
% Returns Nx3 matrix of RGB colors, or empty if not found

    ionColors = [];
    numIons = numel(ionLabels);

    if isempty(specHandle) || ~isgraphics(specHandle)
        return;
    end

    % Initialize with NaN (will use default colors)
    ionColors = nan(numIons, 3);

    try
        % Get the axes containing the spec
        if isgraphics(specHandle, 'axes')
            ax = specHandle;
        else
            ax = ancestor(specHandle, 'axes');
        end

        if isempty(ax)
            return;
        end

        % Find ion annotations (text objects or other annotations)
        % Look for objects with UserData containing ion info
        allChildren = findall(ax);

        for i = 1:numel(allChildren)
            obj = allChildren(i);
            try
                ud = obj.UserData;
                if isstruct(ud)
                    % Check for ion name in UserData
                    if isfield(ud, 'ionName') || isfield(ud, 'ion')
                        if isfield(ud, 'ionName')
                            objIonName = string(ud.ionName);
                        else
                            objIonName = string(ud.ion);
                        end

                        % Find matching ion in our list
                        for j = 1:numIons
                            ionLabel = string(ionLabels(j));
                            labelNoSpace = regexprep(ionLabel, '\s+', '');

                            if strcmpi(objIonName, ionLabel) || strcmpi(objIonName, labelNoSpace)
                                % Get color from object
                                if isprop(obj, 'Color')
                                    ionColors(j, :) = obj.Color;
                                elseif isprop(obj, 'FaceColor') && ~ischar(obj.FaceColor)
                                    ionColors(j, :) = obj.FaceColor;
                                elseif isfield(ud, 'color')
                                    ionColors(j, :) = ud.color;
                                end
                                break;
                            end
                        end
                    end

                    % Also check for range-based color info
                    if isfield(ud, 'rangeName') && isfield(ud, 'color')
                        rangeName = string(ud.rangeName);
                        for j = 1:numIons
                            if strcmpi(rangeName, ionLabels(j))
                                ionColors(j, :) = ud.color;
                                break;
                            end
                        end
                    end
                end
            catch
                % Skip objects that don't have the expected properties
            end
        end

        % Also try to get colors from range patches/areas
        patches = findobj(ax, 'Type', 'patch');
        areas = findobj(ax, 'Type', 'area');
        allPatches = [patches; areas];

        for i = 1:numel(allPatches)
            obj = allPatches(i);
            try
                ud = obj.UserData;
                if isstruct(ud) && isfield(ud, 'rangeName')
                    rangeName = string(ud.rangeName);
                    for j = 1:numIons
                        ionLabel = string(ionLabels(j));
                        if strcmpi(rangeName, ionLabel) || contains(ionLabel, rangeName, 'IgnoreCase', true)
                            if isprop(obj, 'FaceColor') && isnumeric(obj.FaceColor)
                                ionColors(j, :) = obj.FaceColor;
                            end
                            break;
                        end
                    end
                end
            catch
            end
        end
    catch
        % Return empty if anything fails
        ionColors = [];
    end
end

function ionColors = getIonColorsFromScheme(colorScheme, ionLabels)
% Get ion colors from a colorScheme struct
% Returns Nx3 matrix of RGB colors, or empty if not found

    ionColors = [];
    numIons = numel(ionLabels);

    if isempty(colorScheme) || ~isstruct(colorScheme)
        return;
    end

    ionColors = nan(numIons, 3);

    % Check if colorScheme has ions field (from colorSchemeCreate)
    if isfield(colorScheme, 'ions') && istable(colorScheme.ions)
        ionTable = colorScheme.ions;
        if ismember('ionName', ionTable.Properties.VariableNames) && ...
           ismember('color', ionTable.Properties.VariableNames)
            schemeIonNames = string(ionTable.ionName);
            schemeColors = ionTable.color;

            for j = 1:numIons
                ionLabel = string(ionLabels(j));
                labelNoSpace = regexprep(ionLabel, '\s+', '');

                for k = 1:numel(schemeIonNames)
                    if strcmpi(schemeIonNames(k), ionLabel) || strcmpi(schemeIonNames(k), labelNoSpace)
                        if iscell(schemeColors)
                            ionColors(j, :) = schemeColors{k};
                        else
                            ionColors(j, :) = schemeColors(k, :);
                        end
                        break;
                    end
                end
            end
        end
    end

    % Check if all NaN (no colors found)
    if all(isnan(ionColors(:)))
        ionColors = [];
    end
end

function [catsOut, countsOut] = aggregateCounts(cats, counts)
    cats = removecats(cats);
    catNames = categories(cats);
    countsOut = zeros(1, numel(catNames));
    for i = 1:numel(catNames)
        countsOut(i) = sum(counts(cats == catNames{i}));
    end
    catsOut = categorical(catNames);
end

function atoms = buildIsotopeLabels(pos)
    atomStr = string(pos.atom);
    isoStr = string(pos.isotope);
    isAtomMissing = ismissing(atomStr) | atomStr == "" | atomStr == "<undefined>";
    isIsoMissing = ismissing(isoStr) | isoStr == "" | isoStr == "NaN";
    labels = atomStr;
    mask = ~isAtomMissing & ~isIsoMissing;
    labels(mask) = isoStr(mask) + atomStr(mask);
    labels(isAtomMissing) = "unranged";
    atoms = categorical(labels);
end

function [rawCat, rawCounts] = rawCountsFromUniqueHits(pos, rangeTable, mode)
% Extract raw counts from unique detector hits
% Uses pos.ionIdx to identify unique hits, applies ranges, and counts per ion
%
% This gives the true "raw" ion counts before any deconvolution or
% decomposition of molecular ions.

    % Check if ionIdx column exists
    hasIonIdx = ismember('ionIdx', pos.Properties.VariableNames);
    hasIonColumn = ismember('ion', pos.Properties.VariableNames);

    if ~hasIonIdx
        % Fallback: use existing ion column if no ionIdx
        if hasIonColumn
            warning('posCalculateConcentrationDeconvolved:noIonIdx', ...
                'pos.ionIdx not found. Using pos.ion for raw counts (may include decomposed molecular ions).');
            rawCat = pos.ion;
            rawCat(isundefined(rawCat)) = 'unranged';
            rawCat = removecats(rawCat);
            rawCounts = countcats(rawCat);
        else
            % No ionIdx and no ion column - return unranged counts
            warning('posCalculateConcentrationDeconvolved:noIonData', ...
                'Neither pos.ionIdx nor pos.ion found. Cannot determine raw counts per species.');
            rawCat = categorical({'unranged'});
            rawCounts = height(pos);
        end
        return;
    end

    % Get unique detector hits
    [uniqueIdx, firstOccurrence] = unique(pos.ionIdx);
    mcUnique = pos.mc(firstOccurrence);
    numHits = numel(uniqueIdx);

    % Allocate each unique hit to a range
    numRanges = height(rangeTable);
    hitRangeIdx = zeros(numHits, 1);  % 0 = unranged

    for r = 1:numRanges
        inRange = mcUnique >= rangeTable.mcbegin(r) & mcUnique <= rangeTable.mcend(r);
        hitRangeIdx(inRange) = r;
    end

    % Build ion labels from range table
    ionLabels = strings(numHits, 1);
    ionLabels(:) = "unranged";

    for r = 1:numRanges
        mask = hitRangeIdx == r;
        if any(mask)
            ionLabels(mask) = string(rangeTable.rangeName(r));
        end
    end

    % Convert to categorical
    rawIonCat = categorical(ionLabels);

    % Convert to requested mode
    mode = lower(string(mode));
    switch mode
        case "ionic"
            rawCat = rawIonCat;
        otherwise
            % Convert ionic labels to isotopic or atomic
            rawCat = ionConvertMode(rawIonCat, char(mode));
    end

    % Count per category
    rawCat = removecats(rawCat);
    rawCounts = countcats(rawCat);
end

function [element, isotope, chargeState] = parseIonLabel(ionLabel, baseName, defaultCharge)
% Parse an ion label like "2H+", "56Fe++", "H2+", "2H2+" to extract element, isotope, charge
%
% Returns:
%   element: element symbol (e.g., "H", "Fe") - primary element for molecules
%   isotope: mass number (e.g., 2, 56), or 0 if not specified
%   chargeState: charge (e.g., 1, 2)
%
% For molecular ions like "2H2+", returns the element (H) and isotope (2)
% The molecule count (2 in H2) is ignored for reference matching

    ionLabel = string(ionLabel);
    baseName = string(baseName);

    % Default values
    element = "";
    isotope = 0;
    chargeState = defaultCharge;

    labelStr = char(ionLabel);

    % Remove charge symbols and count them
    plusCount = sum(labelStr == '+');
    minusCount = sum(labelStr == '-');
    if plusCount > 0
        chargeState = plusCount;
    elseif minusCount > 0
        chargeState = -minusCount;
    end

    % Remove charge symbols and spaces
    labelStr = regexprep(labelStr, '[+\-\s]', '');

    % Try multiple patterns:
    % 1. Isotope + Element + optional molecule count: "2H", "2H2", "56Fe"
    % 2. Element + molecule count: "H2", "Fe"

    % Pattern for isotope-prefixed: "2H", "2H2", "56Fe"
    tokens = regexp(labelStr, '^(\d+)([A-Z][a-z]?)(\d*)$', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        tok = tokens{1};
        isotope = str2double(tok{1});
        element = string(tok{2});
        return;
    end

    % Pattern for non-isotope-prefixed: "H", "H2", "Fe", "FeO"
    tokens = regexp(labelStr, '^([A-Z][a-z]?)(\d*)(.*)$', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        tok = tokens{1};
        element = string(tok{1});
        % isotope remains 0 (not specified) - will need to find most abundant
        return;
    end

    % Fallback: try to extract element from baseName
    baseStr = char(baseName);
    % Remove numbers and get first element
    tokens = regexp(baseStr, '([A-Z][a-z]?)', 'tokens');
    if ~isempty(tokens)
        element = string(tokens{1}{1});
    end
end

function abundance = getIsotopeAbundance(isotopeTable, element, massNumber)
% Get natural abundance of a specific isotope from the isotope table
%
% Returns abundance as fraction (0-1), or 0 if not found

    abundance = 0;

    if isempty(element) || massNumber == 0
        return;
    end

    element = string(element);

    % Find matching isotope in table
    % Expected columns: element (or symbol), massNumber (or A), abundance

    % Try different column name conventions
    elementCol = '';
    massCol = '';
    abundCol = '';

    colNames = isotopeTable.Properties.VariableNames;

    % Find element column
    for name = ["element", "symbol", "Element", "Symbol"]
        if any(strcmpi(colNames, name))
            elementCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    % Find mass number column
    for name = ["massNumber", "A", "mass_number", "MassNumber", "isotope"]
        if any(strcmpi(colNames, name))
            massCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    % Find abundance column
    for name = ["abundance", "naturalAbundance", "Abundance", "relativeAbundance"]
        if any(strcmpi(colNames, name))
            abundCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    if isempty(elementCol) || isempty(massCol) || isempty(abundCol)
        return;
    end

    % Search for matching isotope
    elements = string(isotopeTable.(elementCol));
    masses = isotopeTable.(massCol);
    abundances = isotopeTable.(abundCol);

    mask = strcmpi(elements, element) & masses == massNumber;
    if any(mask)
        abundance = abundances(find(mask, 1, 'first'));
        % Convert from percentage to fraction if needed
        if abundance > 1
            abundance = abundance / 100;
        end
    end
end

function [catsSorted, countsSorted, countVarSorted, countsRawSorted] = sortCategoriesByMass(cats, counts, countVar, countsRaw, isotopeTable, mode)
% Sort categories by atomic/ionic mass
% Special categories (unranged, tracer) are placed at the end
%
% For atomic mode: sort by atomic mass of the element
% For isotopic mode: sort by isotope mass
% For ionic mode: sort by mass-to-charge ratio

    numCats = numel(cats);
    masses = inf(numCats, 1);  % Default to inf so unknowns go to end
    isSpecial = false(numCats, 1);  % Track unranged and tracer categories

    for i = 1:numCats
        catName = string(cats{i});

        % Check for special categories (put at end)
        if strcmpi(catName, 'unranged') || contains(catName, '(tracer)', 'IgnoreCase', true)
            isSpecial(i) = true;
            % For tracer, extract base name for sorting
            if contains(catName, '(tracer)', 'IgnoreCase', true)
                catName = strtrim(regexprep(catName, '\(tracer\)', '', 'ignorecase'));
            else
                continue;  % unranged stays at inf
            end
        end

        % Parse the category name to get element and mass info
        mass = getCategoryMass(catName, isotopeTable, mode);
        if mass > 0
            masses(i) = mass;
        end
    end

    % Sort: first by isSpecial (false=0 before true=1), then by mass
    [~, sortIdx] = sortrows([isSpecial, masses], [1, 2]);

    catsSorted = cats(sortIdx);
    countsSorted = counts(sortIdx);
    countVarSorted = countVar(sortIdx);
    countsRawSorted = countsRaw(sortIdx);
end

function mass = getCategoryMass(catName, isotopeTable, mode)
% Get the mass associated with a category name
%
% For atomic mode: returns average atomic mass of the element
% For isotopic mode: returns isotope mass
% For ionic mode: returns mass-to-charge ratio

    mass = 0;
    catStr = char(catName);

    % Remove charge symbols and count them
    chargeState = sum(catStr == '+') - sum(catStr == '-');
    if chargeState == 0
        chargeState = 1;  % Assume +1 for atomic/isotopic mode
    end
    catStr = regexprep(catStr, '[+\-\s]', '');

    % Try to parse isotope + element: "56Fe", "2H", etc.
    tokens = regexp(catStr, '^(\d+)?([A-Z][a-z]?)(\d*)$', 'tokens');

    element = '';
    isotope = 0;

    if ~isempty(tokens) && ~isempty(tokens{1})
        tok = tokens{1};
        if ~isempty(tok{1})
            isotope = str2double(tok{1});
        end
        element = tok{2};
    else
        % Try just element: "Fe", "H"
        tokens = regexp(catStr, '^([A-Z][a-z]?)$', 'tokens');
        if ~isempty(tokens) && ~isempty(tokens{1})
            element = tokens{1}{1};
        end
    end

    if isempty(element)
        return;
    end

    % Get mass from isotope table
    mass = getElementMass(isotopeTable, element, isotope);

    % For ionic mode, divide by charge state
    if strcmpi(mode, 'ionic') && chargeState > 0
        mass = mass / chargeState;
    end
end

function mass = getElementMass(isotopeTable, element, isotope)
% Get mass of an element/isotope from the isotope table
% If isotope=0, returns the average atomic mass (abundance-weighted)

    mass = 0;
    element = string(element);

    if isempty(isotopeTable)
        return;
    end

    % Find column names
    colNames = isotopeTable.Properties.VariableNames;

    elementCol = '';
    massNumCol = '';
    massCol = '';
    abundCol = '';

    for name = ["element", "symbol", "Element", "Symbol"]
        if any(strcmpi(colNames, name))
            elementCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    for name = ["massNumber", "A", "mass_number", "MassNumber", "isotope"]
        if any(strcmpi(colNames, name))
            massNumCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    for name = ["mass", "atomicMass", "Mass", "weight"]
        if any(strcmpi(colNames, name))
            massCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    for name = ["abundance", "naturalAbundance", "Abundance", "relativeAbundance"]
        if any(strcmpi(colNames, name))
            abundCol = colNames{strcmpi(colNames, name)};
            break;
        end
    end

    if isempty(elementCol) || isempty(massCol)
        return;
    end

    elements = string(isotopeTable.(elementCol));
    masses = isotopeTable.(massCol);

    % Find rows for this element
    elementMask = strcmpi(elements, element);

    if ~any(elementMask)
        return;
    end

    if isotope > 0 && ~isempty(massNumCol)
        % Specific isotope requested
        massNums = isotopeTable.(massNumCol);
        mask = elementMask & massNums == isotope;
        if any(mask)
            mass = masses(find(mask, 1, 'first'));
        end
    else
        % Return average atomic mass (abundance-weighted) or most abundant isotope mass
        if ~isempty(abundCol)
            abundances = isotopeTable.(abundCol);
            % Normalize abundances if in percentage
            if any(abundances > 1)
                abundances = abundances / 100;
            end
            elemMasses = masses(elementMask);
            elemAbund = abundances(elementMask);
            % Abundance-weighted average
            if sum(elemAbund) > 0
                mass = sum(elemMasses .* elemAbund) / sum(elemAbund);
            else
                mass = mean(elemMasses);
            end
        else
            % Just use average of isotope masses
            mass = mean(masses(elementMask));
        end
    end
end

%% ========== Background Plotting ==========

function h = plotBackgroundOnMassSpec(massSpec, mcCenters, bg, binWidth, totalIons)
% Plot background on an existing mass spectrum
% massSpec can be an area plot handle, axes, or figure

    h = gobjects(0, 1);

    % Resolve the area plot and axes
    [ax, areaPlot] = resolveBgPlotHandle(massSpec);

    if isempty(ax) || ~isgraphics(ax, 'axes')
        return;
    end

    % Get the spectrum's bin width from the area plot if available
    specBinWidth = binWidth;
    if ~isempty(areaPlot) && isgraphics(areaPlot)
        specX = areaPlot.XData;
        specBinWidth = median(diff(specX));
        if ~isfinite(specBinWidth) || specBinWidth <= 0
            specBinWidth = binWidth;
        end
    end

    % Detect spectrum units from y-axis label
    specUnits = detectBgPlotUnits(ax);

    % Scale background to match spectrum units
    % bg is in counts/Da/totalIons
    if specUnits == "normalised"
        % Normalized spectrum shows counts/Da/totalIons - bg matches directly
        bgPlot = bg;
    else
        % Spectrum is counts per bin - convert bg to counts/bin
        bgPlot = bg * specBinWidth * totalIons;
    end

    % Handle log scale
    if strcmpi(ax.YScale, 'log')
        minY = ax.YLim(1);
        if minY <= 0
            minY = min(bgPlot(bgPlot > 0));
        end
        bgPlot(bgPlot < minY) = minY;
    end

    % Plot
    holdState = ishold(ax);
    hold(ax, 'on');
    h = plot(ax, mcCenters, bgPlot, 'LineWidth', 1.5, 'Color', [1 0 0], 'LineStyle', '--');
    h.DisplayName = 'background';
    h.UserData.plotType = "backgroundEstimate";
    if ~holdState
        hold(ax, 'off');
    end
end

function [ax, areaPlot] = resolveBgPlotHandle(massSpec)
% Resolve mass spectrum handle to axes and area plot for background plotting

    ax = [];
    areaPlot = [];

    if isgraphics(massSpec, 'axes')
        ax = massSpec;
        plots = findobj(ax, 'Type', 'area');
        areaPlot = pickBgMassSpecPlot(plots);
    elseif isgraphics(massSpec, 'figure')
        plots = findobj(massSpec, 'Type', 'area');
        areaPlot = pickBgMassSpecPlot(plots);
        if ~isempty(areaPlot)
            ax = ancestor(areaPlot, 'axes');
        end
    elseif isgraphics(massSpec, 'area')
        % Direct area plot handle
        areaPlot = massSpec;
        ax = ancestor(massSpec, 'axes');
    elseif isgraphics(massSpec)
        % Other graphics object - try to get axes
        ax = ancestor(massSpec, 'axes');
        if ~isempty(ax)
            plots = findobj(ax, 'Type', 'area');
            areaPlot = pickBgMassSpecPlot(plots);
        end
    end
end

function spec = pickBgMassSpecPlot(plots)
% Pick the mass spectrum area plot from a list of area plots

    spec = [];
    for i = 1:numel(plots)
        try
            if isfield(plots(i).UserData, 'plotType') && plots(i).UserData.plotType == "massSpectrum"
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

function units = detectBgPlotUnits(ax)
    units = "count";
    try
        ylab = lower(string(ax.YLabel.String));
        if contains(ylab, "cts / da") || contains(ylab, "cts/da") || contains(ylab, "totcts") || contains(ylab, "/da")
            units = "normalised";
        end
    catch
    end
end
