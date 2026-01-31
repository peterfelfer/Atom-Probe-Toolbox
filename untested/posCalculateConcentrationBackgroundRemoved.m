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
%   'method'        'linearBetweenPeaks' | 'massSpecInvSqrt' (default: 'linearBetweenPeaks')
%   'mode'          'ionic' | 'isotopic' | 'atomic'
%                   - 'ionic': each unique ion species counted separately
%                   - 'isotopic': groups by isotopic composition (ignores charge)
%                   - 'atomic': groups by element (ignores isotope and charge)
%                   Default: 'atomic' if 'atom' column exists, 'ionic' otherwise
%   'rangeTable'    range table (optional if massSpec is provided)
%   'massSpec'      mass spectrum area plot handle (used to extract ranges and for plotting)
%   'bin'           histogram bin width for internal calculations [Da] (default: 0.01)
%   'plotBackground' true/false (default: false)
%   'sortBy'        'auto' | 'none' | 'atomic' | 'weight' (default: 'auto')
%   'minPeakDistance' minimum distance from peaks for background fitting [Da] (default: 0.3)
%   'fitLimits'     Nx2 matrix of [begin, end] m/c ranges to use for fitting (optional)
%                   Example: [5, 20; 40, 60] fits only in 5-20 Da and 40-60 Da
%
% OUTPUT
% conc:   concentration table (same format as posCalculateConcentrationSimple)
% info:   struct with background details
%
% METHODS
% 'linearBetweenPeaks' - Linear interpolation between gaps
%         Fits linear models in each gap and interpolates under peaks.
%
% 'massSpecInvSqrt' - Fits B(m/c) = A / sqrt(m/c) to unranged regions
%         Physics-based model: constant background in TOF space transforms
%         to 1/sqrt(m/c) in mass spectrum space.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Input validation
if detEff > 1
    detEff = detEff / 100;
end

if ~exist('volumeName', 'var') || isempty(volumeName)
    volumeName = inputname(1);
end

if ~exist('excludeList', 'var')
    excludeList = {};
end

%% Parse options
options = struct( ...
    'method', 'linearBetweenPeaks', ...
    'mode', '', ...
    'rangeTable', [], ...
    'massSpec', [], ...
    'bin', 0.01, ...
    'plotBackground', false, ...
    'sortBy', 'auto', ...
    'minPeakDistance', 0.3, ...
    'fitLimits', []);

for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "method"
            options.method = char(value);
        case "mode"
            options.mode = strtrim(char(value));
        case "rangetable"
            options.rangeTable = value;
        case "massspec"
            options.massSpec = value;
        case "bin"
            options.bin = value;
        case "plotbackground"
            options.plotBackground = logical(value);
        case "sortby"
            options.sortBy = char(value);
        case "minpeakdistance"
            options.minPeakDistance = value;
        case "fitlimits"
            options.fitLimits = value;
        otherwise
            error('posCalculateConcentrationBackgroundRemoved:invalidOption', ...
                'Unknown option "%s".', name);
    end
end

%% Get range table (from options or extract from mass spectrum)
rangeTable = options.rangeTable;
if isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0
    % Try to extract from mass spectrum
    if ~isempty(options.massSpec)
        rangeTable = rangesExtractFromMassSpec(options.massSpec);
    end
end

if isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0
    error('posCalculateConcentrationBackgroundRemoved:missingRangeTable', ...
        'No ranges found. Provide rangeTable or a massSpec with ranges.');
end

%% Determine concentration mode
hasAtomColumn = any(ismember(pos.Properties.VariableNames, 'atom'));
hasIonColumn = any(ismember(pos.Properties.VariableNames, 'ion'));
requestedMode = lower(string(options.mode));
mode = requestedMode;

if isempty(mode) || mode == ""
    if hasAtomColumn
        type = 'atomic';
        columnType = 'atom';
        atoms = pos.atom;
    else
        type = 'ionic';
        columnType = 'ion';
        atoms = pos.ion;
    end
else
    switch mode
        case "ionic"
            if ~hasIonColumn
                error('posCalculateConcentrationBackgroundRemoved:missingColumn', ...
                    'Ionic mode requires ''ion'' column in pos table.');
            end
            type = 'ionic';
            columnType = 'ion';
            atoms = pos.ion;
        case "isotopic"
            type = 'isotopic';
            if hasAtomColumn && ismember('isotope', pos.Properties.VariableNames)
                columnType = 'isotope';
                atoms = buildIsotopeLabels(pos);
            elseif hasIonColumn
                columnType = 'ion';
                atoms = ionConvertMode(pos.ion, 'isotopic');
            else
                error('posCalculateConcentrationBackgroundRemoved:missingColumn', ...
                    'Isotopic mode requires ''atom''+''isotope'' or ''ion'' column.');
            end
        case "atomic"
            if hasAtomColumn
                type = 'atomic';
                columnType = 'atom';
                atoms = pos.atom;
            elseif hasIonColumn
                type = 'atomic';
                columnType = 'ion';
                atoms = ionConvertMode(pos.ion, 'atomic');
            else
                error('posCalculateConcentrationBackgroundRemoved:missingColumn', ...
                    'Atomic mode requires ''atom'' or ''ion'' column.');
            end
        otherwise
            error('posCalculateConcentrationBackgroundRemoved:invalidMode', ...
                'Mode must be ''ionic'', ''isotopic'', or ''atomic''.');
    end
end
atoms(isundefined(atoms)) = 'unranged';

%% Build internal histogram from pos data
mc = pos.mc;
mcMax = max(mc);
mcMin = max(0.5, min(mc));
binWidth = options.bin;
edges = mcMin:binWidth:mcMax;
counts = histcounts(mc, edges);
mcCenters = edges(1:end-1) + binWidth/2;

% Normalize to counts/Da/totalIons for fitting
totalIons = height(pos);
countsNorm = counts / (binWidth * totalIons);

%% Compute background using shared helper function
[bg, bgInfo] = backgroundEstimate(mcCenters, countsNorm, rangeTable, ...
    'method', options.method, ...
    'minPeakDistance', options.minPeakDistance, ...
    'fitLimits', options.fitLimits);

info = struct();
info.method = bgInfo.method;

% Copy method-specific fields
if isfield(bgInfo, 'coefficient')
    info.fitCoefficient = bgInfo.coefficient;
end
if isfield(bgInfo, 'rsquared')
    info.fitRsquared = bgInfo.rsquared;
end
if isfield(bgInfo, 'fitRegions')
    info.fitRegions = bgInfo.fitRegions;
end
if isfield(bgInfo, 'gapFits')
    info.gapFits = bgInfo.gapFits;
end
info.fitLimits = options.fitLimits;

% Convert background to counts per bin
bgCounts_perBin = bg * binWidth * totalIons;

info.mcCenters = mcCenters;
info.countsNorm = countsNorm;
info.background = bgCounts_perBin;  % counts per bin
info.backgroundNorm = bg;            % counts/Da/totalIons (for plotting on normalized spectra)
info.binWidth = binWidth;
info.totalIons = totalIons;

%% Plot background if requested
bgPlotHandle = gobjects(0, 1);
if options.plotBackground && ~isempty(options.massSpec)
    % Pass normalized bg for plotting (will be scaled to match spectrum units)
    bgPlotHandle = plotBackgroundOnMassSpec(options.massSpec, mcCenters, bg, binWidth, totalIons);
    info.backgroundPlot = bgPlotHandle;
end

%% Compute background counts per range using analytical integral for invSqrt
bgCountsByIon = zeros(numel(categories(atoms)), 1);
cats = categories(atoms);

if method == "massspecinvsqrt" && isfield(info, 'fitCoefficient')
    % Use analytical integral: ∫ a/√mc dmc = 2a√mc
    coeff = info.fitCoefficient;
    for i = 1:height(rangeTable)
        mcBegin = rangeTable.mcbegin(i);
        mcEnd = rangeTable.mcend(i);
        % Integral of coeff/sqrt(mc) from mcBegin to mcEnd, scaled to counts
        bgCount = 2 * coeff * (sqrt(mcEnd) - sqrt(mcBegin)) * totalIons;

        % Get ion/atom label for this range
        label = getLabelForRange(rangeTable, i, type, columnType);
        idx = find(string(cats) == label, 1, 'first');
        if ~isempty(idx) && label ~= "unranged"
            bgCountsByIon(idx) = bgCountsByIon(idx) + bgCount;
        end
    end
else
    % Use histogram summing for linearBetweenPeaks
    for i = 1:height(rangeTable)
        mcBegin = rangeTable.mcbegin(i);
        mcEnd = rangeTable.mcend(i);
        inRange = mcCenters >= mcBegin & mcCenters <= mcEnd;
        if ~any(inRange)
            continue;
        end
        % bg is in counts/Da/totalIons, sum and convert to counts
        bgCount = sum(bg(inRange)) * binWidth * totalIons;

        label = getLabelForRange(rangeTable, i, type, columnType);
        idx = find(string(cats) == label, 1, 'first');
        if ~isempty(idx) && label ~= "unranged"
            bgCountsByIon(idx) = bgCountsByIon(idx) + bgCount;
        end
    end
end

%% Compute corrected counts and concentrations
rawCounts = countcats(atoms);
if iscolumn(rawCounts)
    rawCounts = rawCounts';
end
bgCounts = bgCountsByIon';
if numel(bgCounts) < numel(rawCounts)
    bgCounts(end+1:numel(rawCounts)) = 0;
end

% Sort columns
sortMode = options.sortBy;
if strcmpi(sortMode, 'auto')
    if strcmpi(type, 'atomic') || strcmpi(type, 'isotopic')
        sortMode = 'atomic';
    else
        sortMode = 'weight';
    end
end
[cats, rawCounts, bgCounts] = sortIonCategories(cats, sortMode, rawCounts, bgCounts);
bgCounts(strcmpi(cats, 'unranged')) = 0;

corrCounts = rawCounts - bgCounts;
corrCounts(corrCounts < 0) = 0;

% Exclusion handling
isExcluded = ismember(cats, excludeList)';

% Uncorrected concentration (without background removal)
concFracUncorr = rawCounts ./ sum(rawCounts(~isExcluded));
concFracUncorr(isExcluded) = 0;

% Corrected concentration (with background removal)
concFracCorr = corrCounts ./ sum(corrCounts(~isExcluded));
concFracCorr(isExcluded) = 0;

% Background fraction (background / raw counts)
bgFraction = bgCounts ./ max(rawCounts, 1);
bgFraction(rawCounts == 0) = 0;

% Variance for corrected concentration
variance = concFracCorr .* (1 - concFracCorr) ./ max(corrCounts, 1) * (1 - detEff);

%% Build output table
countsOut = [rawCounts; bgCounts; bgFraction; corrCounts; concFracUncorr; concFracCorr; variance];
conc = array2table(countsOut, 'VariableNames', cats');
conc.Properties.VariableDescriptions = repmat({columnType}, size(cats'));
nRows = 7;
volumeCat = categorical(repmat(string(volumeName), nRows, 1));
conc = [table(volumeCat, 'VariableNames', {'volume'}), ...
    table(zeros(nRows, 1), 'VariableNames', {'distance'}), ...
    table(categorical(repmat({type}, nRows, 1)), 'VariableNames', {'type'}), ...
    table(categorical({'count'; 'backgroundCounts'; 'backgroundFraction'; 'countCorrected'; 'concentration'; 'concentrationCorrected'; 'variance'}), 'VariableNames', {'format'}), ...
    conc];

info.rangeTable = rangeTable;
info.rawCounts = rawCounts;
info.backgroundCounts = bgCounts;
info.correctedCounts = corrCounts;
info.categories = cats;

end

%% ========== Plotting ==========

function h = plotBackgroundOnMassSpec(massSpec, mcCenters, bg, binWidth, totalIons)
% Plot background on an existing mass spectrum
% massSpec can be an area plot handle, axes, or figure

    h = gobjects(0, 1);

    % Resolve the area plot and axes
    [ax, areaPlot] = resolveSpectrumHandle(massSpec);

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
    specUnits = detectSpectrumUnits(ax);

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

function [ax, areaPlot] = resolveSpectrumHandle(massSpec)
% Resolve mass spectrum handle to axes and area plot

    ax = [];
    areaPlot = [];

    if isgraphics(massSpec, 'axes')
        ax = massSpec;
        plots = findobj(ax, 'Type', 'area');
        areaPlot = pickMassSpecPlot(plots);
    elseif isgraphics(massSpec, 'figure')
        plots = findobj(massSpec, 'Type', 'area');
        areaPlot = pickMassSpecPlot(plots);
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
            areaPlot = pickMassSpecPlot(plots);
        end
    end
end

function spec = pickMassSpecPlot(plots)
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

function units = detectSpectrumUnits(ax)
    units = "count";
    try
        ylab = lower(string(ax.YLabel.String));
        if contains(ylab, "cts / da") || contains(ylab, "cts/da") || contains(ylab, "totcts") || contains(ylab, "/da")
            units = "normalised";
        end
    catch
    end
end

%% ========== Helper Functions ==========

function label = getLabelForRange(rangeTable, idx, mode, columnType)
% Get the appropriate label for a range based on mode

    label = "";

    % Try to get from ion table in range
    if ismember('ion', rangeTable.Properties.VariableNames)
        ionEntry = rangeTable.ion{idx};
        if istable(ionEntry) && ~isempty(ionEntry)
            if mode == "atomic" && ismember('element', ionEntry.Properties.VariableNames)
                elements = string(ionEntry.element);
                elements = elements(elements ~= "" & elements ~= "<undefined>");
                if ~isempty(elements)
                    label = elements(1);  % Take first element
                    return;
                end
            elseif mode == "isotopic"
                if all(ismember({'element', 'isotope'}, ionEntry.Properties.VariableNames))
                    el = string(ionEntry.element(1));
                    iso = ionEntry.isotope(1);
                    if el ~= "" && el ~= "<undefined>" && ~isnan(iso)
                        label = string(iso) + el;
                        return;
                    end
                end
            end
            % Fall back to ion name
            if ismember('chargeState', rangeTable.Properties.VariableNames)
                label = string(ionConvertName(ionEntry, rangeTable.chargeState(idx)));
            else
                label = string(ionConvertName(ionEntry));
            end
            if mode ~= "ionic"
                label = string(ionConvertMode(categorical(label), mode));
            end
            return;
        end
    end

    % Fall back to rangeName
    if ismember('rangeName', rangeTable.Properties.VariableNames)
        label = string(rangeTable.rangeName(idx));
    end
end

function atoms = buildIsotopeLabels(pos)
% Build isotope labels from atom and isotope columns

    atomCol = pos.atom;
    isoCol = pos.isotope;

    if iscategorical(atomCol)
        atomStr = string(atomCol);
    else
        atomStr = string(atomCol);
    end

    if isnumeric(isoCol)
        isoStr = strings(size(isoCol));
        validIso = ~isnan(isoCol) & isoCol > 0;
        isoStr(validIso) = string(isoCol(validIso));
    elseif iscategorical(isoCol)
        isoStr = string(isoCol);
    else
        isoStr = string(isoCol);
    end

    isAtomMissing = ismissing(atomStr) | atomStr == "" | atomStr == "<undefined>";
    isIsoMissing = ismissing(isoStr) | isoStr == "" | isoStr == "NaN" | isoStr == "0";

    labels = atomStr;
    mask = ~isAtomMissing & ~isIsoMissing;
    labels(mask) = isoStr(mask) + atomStr(mask);
    labels(isAtomMissing) = "unranged";

    atoms = categorical(labels);
end
