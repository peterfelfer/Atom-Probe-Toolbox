function results = peakAnalysis(massSpec, options)
% PEAKANALYSIS Background-corrected peak analysis for mass spectra
%
% results = peakAnalysis(massSpec)
% results = peakAnalysis(massSpec, 'peaks', peakRanges)
% results = peakAnalysis(massSpec, 'rangeTable', ranges)
%
% Performs background-corrected peak counting on mass spectra. Works with
% both raw count and normalized spectra. Supports multiple peaks and
% calculates fractions relative to all ranged atoms when a range table
% is provided.
%
% INPUT:
%   massSpec - Mass spectrum data, one of:
%              - Figure handle to a mass spectrum plot (from massSpecPlot)
%              - Struct with fields: .mc (mass values), .counts (bin counts)
%              - Struct with fields: .mc, .normalized (for normalized spectra)
%              - Table with 'mc' and 'counts' (or 'normalized') columns
%
% OPTIONS:
%   'peaks'        - Nx4 array of peak regions: [bgStart1 peakStart1 peakEnd1 bgEnd1; ...]
%                    Each row defines: [background_start, peak_start, peak_end, background_end]
%                    If not provided, interactive selection is enabled
%   'rangeTable'   - Range table with ranged ions. When provided, calculates
%                    fraction of peak counts relative to all ranged atoms
%   'bgMethod'     - Background fitting method (default: 'linear')
%                    'linear'     - Linear least squares fit
%                    'exponential'- Exponential decay fit
%                    'constant'   - Constant background (average)
%   'binWidth'     - Bin width for internal histogram (default: 0.01 Da)
%                    Only used when massSpec is a figure handle
%   'showPlot'     - Display results plot (default: true)
%   'showSelection'- Show interactive peak selection guide (default: true)
%   'normalized'   - Treat input as normalized spectrum (default: auto-detect)
%
% OUTPUT:
%   results - Struct array (one element per peak) with fields:
%       .peakRange      - [start, end] of peak region in Da
%       .bgRange        - [start1, end1, start2, end2] of background regions
%       .rawCounts      - Total counts in peak region (before BG correction)
%       .bgCounts       - Estimated background counts in peak region
%       .correctedCounts- Background-corrected counts (rawCounts - bgCounts)
%       .peakLocation   - Mass-to-charge of peak maximum
%       .fractionTotal  - Fraction of all atoms (corrected counts / total counts)
%       .fractionRanged - Fraction of ranged atoms (if rangeTable provided)
%       .ppm            - Concentration in ppm (fractionTotal * 1e6)
%       .uncertainty    - Counting statistics uncertainty (Poisson)
%       .bgFit          - Background fit parameters
%       .signalToBg     - Signal-to-background ratio in peak region
%
% INTERACTIVE SELECTION:
%   When 'peaks' is not provided, the function enables interactive selection:
%   1. Click to define background region BEFORE peak (2 clicks: start, end)
%   2. Click to define the PEAK region (2 clicks: start, end)
%   3. Click to define background region AFTER peak (2 clicks: start, end)
%   4. Press ENTER to finish or continue adding more peaks
%   5. Press ESC to cancel
%
% EXAMPLES:
%   % Interactive selection on current mass spectrum figure
%   results = peakAnalysis(gcf);
%
%   % Analyze specific peak with predefined regions
%   results = peakAnalysis(massSpec, 'peaks', [26.5 27 28 28.5]);
%
%   % Multiple peaks
%   peaks = [26.5 27 28 28.5; 53.5 54 55 55.5];
%   results = peakAnalysis(massSpec, 'peaks', peaks);
%
%   % With range table for fraction calculation
%   results = peakAnalysis(massSpec, 'rangeTable', ranges, 'peaks', [1.5 2 2.5 3]);
%
%   % Using normalized spectrum
%   results = peakAnalysis(massSpec, 'normalized', true);
%
% SEE ALSO:
%   massSpecPlot, rangeAdd, peakBGCorrectedCount
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    massSpec
    options.peaks (:,4) double = []
    options.rangeTable table = table()
    options.bgMethod (1,:) char {mustBeMember(options.bgMethod, {'linear', 'exponential', 'constant'})} = 'linear'
    options.binWidth (1,1) double {mustBePositive} = 0.01
    options.showPlot (1,1) logical = true
    options.showSelection (1,1) logical = true
    options.normalized (1,1) logical = false
    options.autoDetectNormalized (1,1) logical = true
end

%% Extract mass spectrum data
[mc, counts, isNormalized, totalAtoms, binWidth] = extractMassSpecData(massSpec, options);

%% Get peak regions (interactive or from input)
if isempty(options.peaks)
    peakRegions = interactiveSelection(mc, counts, options.showSelection);
    if isempty(peakRegions)
        results = [];
        return;
    end
else
    peakRegions = options.peaks;
end

nPeaks = size(peakRegions, 1);

%% Calculate ranged atom count if range table provided
rangedAtomCount = 0;
if ~isempty(options.rangeTable) && height(options.rangeTable) > 0
    rangedAtomCount = calculateRangedAtoms(mc, counts, options.rangeTable, binWidth);
end

%% Analyze each peak
results = struct([]);

for p = 1:nPeaks
    bgStart = peakRegions(p, 1);
    peakStart = peakRegions(p, 2);
    peakEnd = peakRegions(p, 3);
    bgEnd = peakRegions(p, 4);

    % Extract data for this peak region
    inRegion = mc >= bgStart & mc <= bgEnd;
    mcRegion = mc(inRegion);
    countsRegion = counts(inRegion);

    % Define background and peak masks
    inBgBefore = mcRegion >= bgStart & mcRegion < peakStart;
    inBgAfter = mcRegion > peakEnd & mcRegion <= bgEnd;
    inPeak = mcRegion >= peakStart & mcRegion <= peakEnd;
    inBg = inBgBefore | inBgAfter;

    % Fit background
    [bgFit, bgCounts, fitCurve] = fitBackground(mcRegion, countsRegion, inBg, inPeak, options.bgMethod);

    % Calculate peak statistics
    rawCounts = sum(countsRegion(inPeak));
    correctedCounts = rawCounts - bgCounts;

    % For normalized spectra, scale back to approximate counts
    if isNormalized && totalAtoms > 0
        correctedCounts = correctedCounts * totalAtoms * binWidth;
        rawCounts = rawCounts * totalAtoms * binWidth;
        bgCounts = bgCounts * totalAtoms * binWidth;
    end

    % Find peak location
    peakMc = mcRegion(inPeak);
    peakCounts = countsRegion(inPeak);
    [~, maxIdx] = max(peakCounts);
    peakLocation = peakMc(maxIdx);

    % Calculate fractions
    fractionTotal = correctedCounts / totalAtoms;
    if rangedAtomCount > 0
        fractionRanged = correctedCounts / rangedAtomCount;
    else
        fractionRanged = NaN;
    end

    % Uncertainty (Poisson statistics)
    uncertainty = sqrt(rawCounts + bgCounts) / totalAtoms;

    % Signal to background ratio
    if bgCounts > 0
        signalToBg = correctedCounts / bgCounts;
    else
        signalToBg = Inf;
    end

    % Store results
    res = struct();
    res.peakRange = [peakStart, peakEnd];
    res.bgRange = [bgStart, peakStart, peakEnd, bgEnd];
    res.rawCounts = rawCounts;
    res.bgCounts = bgCounts;
    res.correctedCounts = correctedCounts;
    res.peakLocation = peakLocation;
    res.fractionTotal = fractionTotal;
    res.fractionRanged = fractionRanged;
    res.ppm = fractionTotal * 1e6;
    res.uncertainty = uncertainty;
    res.bgFit = bgFit;
    res.signalToBg = signalToBg;

    % Store fit data for plotting
    res.mcRegion = mcRegion;
    res.countsRegion = countsRegion;
    res.fitCurve = fitCurve;
    res.inPeak = inPeak;
    res.inBg = inBg;

    if isempty(results)
        results = res;
    else
        results(end+1) = res;
    end
end

%% Plot results
if options.showPlot
    plotResults(results, isNormalized, rangedAtomCount > 0);
end

end

%% Helper Functions

function [mc, counts, isNormalized, totalAtoms, binWidth] = extractMassSpecData(massSpec, options)
    % Extract mass spectrum data from various input formats

    isNormalized = options.normalized;
    binWidth = options.binWidth;
    totalAtoms = 0;

    if ishandle(massSpec) && isgraphics(massSpec)
        % Figure or axes handle
        if strcmp(get(massSpec, 'Type'), 'figure')
            ax = get(massSpec, 'CurrentAxes');
        else
            ax = massSpec;
        end

        % Find the mass spectrum data in the axes
        children = get(ax, 'Children');
        for i = 1:length(children)
            if isprop(children(i), 'XData') && isprop(children(i), 'YData')
                mc = get(children(i), 'XData');
                counts = get(children(i), 'YData');

                % Check if this looks like a mass spectrum
                if length(mc) > 100 && max(mc) > 1
                    break;
                end
            end
        end

        % Detect if normalized based on Y-axis label or values
        yLabel = get(get(ax, 'YLabel'), 'String');
        if contains(lower(yLabel), 'normalised') || contains(lower(yLabel), 'normalized')
            isNormalized = true;
        end
        if options.autoDetectNormalized && max(counts) < 1
            isNormalized = true;
        end

        % Estimate total atoms
        if isNormalized
            totalAtoms = 1 / (max(counts) * binWidth);  % Rough estimate
        else
            totalAtoms = sum(counts);
        end

        % Estimate bin width from data
        if length(mc) > 1
            binWidth = mean(diff(mc));
        end

    elseif istable(massSpec)
        % Table input
        mc = massSpec.mc;
        if ismember('counts', massSpec.Properties.VariableNames)
            counts = massSpec.counts;
        elseif ismember('normalized', massSpec.Properties.VariableNames)
            counts = massSpec.normalized;
            isNormalized = true;
        else
            error('peakAnalysis:invalidTable', 'Table must have ''counts'' or ''normalized'' column.');
        end

        totalAtoms = sum(counts);
        if length(mc) > 1
            binWidth = mean(diff(mc));
        end

    elseif isstruct(massSpec)
        % Struct input
        mc = massSpec.mc;
        if isfield(massSpec, 'counts')
            counts = massSpec.counts;
        elseif isfield(massSpec, 'normalized')
            counts = massSpec.normalized;
            isNormalized = true;
        elseif isfield(massSpec, 'cts')
            counts = massSpec.cts;
        else
            error('peakAnalysis:invalidStruct', 'Struct must have ''counts'', ''cts'', or ''normalized'' field.');
        end

        if isfield(massSpec, 'totalAtoms')
            totalAtoms = massSpec.totalAtoms;
        else
            totalAtoms = sum(counts);
        end

        if length(mc) > 1
            binWidth = mean(diff(mc));
        end

    else
        error('peakAnalysis:invalidInput', 'massSpec must be a figure handle, table, or struct.');
    end

    % Ensure column vectors
    mc = mc(:)';
    counts = counts(:)';
end

function peakRegions = interactiveSelection(mc, counts, showGuide)
    % Interactive peak region selection with step-by-step visual guidance

    peakRegions = [];

    % Create or use existing figure
    fig = gcf;
    ax = gca;

    % Plot the spectrum if not already plotted
    if isempty(get(ax, 'Children'))
        area(ax, mc, counts, 'FaceColor', [0.9 0.9 0.9]);
        set(ax, 'YScale', 'log');
        xlabel('mass-to-charge-state ratio [Da]');
        ylabel('counts');
    end

    hold(ax, 'on');

    % Define step labels and colors
    stepLabels = {
        'Click: Background region START (left of peak)', ...
        'Click: Background region END (where peak begins)', ...
        'Click: Peak region START', ...
        'Click: Peak region END', ...
        'Click: Background region START (right of peak)', ...
        'Click: Background region END'
    };
    stepColors = {
        [0.3 0.3 0.8], ...  % Blue - bg before start
        [0.3 0.3 0.8], ...  % Blue - bg before end
        [0.8 0.2 0.2], ...  % Red - peak start
        [0.8 0.2 0.2], ...  % Red - peak end
        [0.2 0.6 0.2], ...  % Green - bg after start
        [0.2 0.6 0.2]       % Green - bg after end
    };
    stepShortLabels = {'BG start', 'BG end', 'Peak start', 'Peak end', 'BG start', 'BG end'};

    % Create instruction text box
    if showGuide
        instructionBox = annotation(fig, 'textbox', [0.25, 0.92, 0.5, 0.06], ...
            'String', stepLabels{1}, ...
            'BackgroundColor', [1 1 0.85], 'EdgeColor', 'k', ...
            'HorizontalAlignment', 'center', 'FontSize', 10, ...
            'FontWeight', 'bold', 'FitBoxToText', 'off');

        % Add help text at bottom
        helpBox = annotation(fig, 'textbox', [0.15, 0.01, 0.7, 0.04], ...
            'String', 'Press ENTER when done adding peaks  |  Press ESC to cancel', ...
            'BackgroundColor', [0.95 0.95 0.95], 'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center', 'FontSize', 8);

        % Add progress indicator
        progressBox = annotation(fig, 'textbox', [0.01, 0.92, 0.15, 0.06], ...
            'String', 'Step 1/6', ...
            'BackgroundColor', 'w', 'EdgeColor', 'k', ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end

    selecting = true;
    peakCount = 0;

    while selecting
        points = zeros(1, 6);
        tempLines = [];
        tempTexts = [];

        try
            for i = 1:6
                % Update instruction text
                if showGuide
                    set(instructionBox, 'String', stepLabels{i});
                    set(instructionBox, 'BackgroundColor', [stepColors{i} 0.3] + [0.7 0.7 0.7 0]);
                    set(progressBox, 'String', sprintf('Step %d/6', i));
                end

                [x, ~, button] = ginput(1);

                if isempty(button)
                    % Enter pressed - finish selection
                    selecting = false;
                    break;
                elseif button == 27
                    % Escape pressed - cancel
                    selecting = false;
                    peakRegions = [];
                    break;
                end

                points(i) = x;

                % Draw vertical line marker at clicked position
                yLim = get(ax, 'YLim');
                tempLines(end+1) = plot(ax, [x x], yLim, '-', ...
                    'Color', stepColors{i}, 'LineWidth', 2);

                % Add label at top of line
                yPos = 10^(log10(yLim(1)) + 0.95*(log10(yLim(2)) - log10(yLim(1))));
                tempTexts(end+1) = text(ax, x, yPos, stepShortLabels{i}, ...
                    'Color', stepColors{i}, 'FontSize', 8, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                    'BackgroundColor', 'w', 'EdgeColor', stepColors{i}, 'Margin', 1);

                % Shade completed regions
                if i == 2 && points(1) > 0
                    % Shade background before region
                    patch(ax, [points(1) points(2) points(2) points(1)], ...
                        [yLim(1) yLim(1) yLim(2) yLim(2)], stepColors{1}, ...
                        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'Tag', 'tempPatch');
                elseif i == 4 && points(3) > 0
                    % Shade peak region
                    patch(ax, [points(3) points(4) points(4) points(3)], ...
                        [yLim(1) yLim(1) yLim(2) yLim(2)], stepColors{3}, ...
                        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'Tag', 'tempPatch');
                elseif i == 6 && points(5) > 0
                    % Shade background after region
                    patch(ax, [points(5) points(6) points(6) points(5)], ...
                        [yLim(1) yLim(1) yLim(2) yLim(2)], stepColors{5}, ...
                        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'Tag', 'tempPatch');
                end
            end

            if all(points > 0)
                % Sort points to get proper order
                sortedPoints = sort(points);

                % Map to: [bgStart, peakStart, peakEnd, bgEnd]
                % After sorting: 1=leftmost bg, 2=bg/peak boundary, 3&4=peak, 5=peak/bg boundary, 6=rightmost bg
                % The user clicked in order, so after sorting:
                % sortedPoints(1) = bg before start
                % sortedPoints(2) = bg before end / peak boundary
                % sortedPoints(3) = peak start (could overlap with 2)
                % sortedPoints(4) = peak end
                % sortedPoints(5) = bg after start (could overlap with 4)
                % sortedPoints(6) = bg after end

                peakRegions(end+1, :) = [sortedPoints(1), sortedPoints(3), sortedPoints(4), sortedPoints(6)];

                peakCount = peakCount + 1;

                % Update instruction for next peak
                if showGuide
                    set(instructionBox, 'String', sprintf('Peak %d added! Click to add another peak or press ENTER to finish', peakCount));
                    set(instructionBox, 'BackgroundColor', [0.8 1 0.8]);
                    set(progressBox, 'String', sprintf('%d peak(s)', peakCount));
                end

                % Clean up temporary elements but keep final shading
                delete(tempLines(ishandle(tempLines)));
                delete(tempTexts(ishandle(tempTexts)));
                tempLines = [];
                tempTexts = [];

                % Keep a subtle permanent highlight for this peak
                yLim = get(ax, 'YLim');
                patch(ax, [sortedPoints(3) sortedPoints(4) sortedPoints(4) sortedPoints(3)], ...
                    [yLim(1) yLim(1) yLim(2) yLim(2)], 'r', ...
                    'FaceAlpha', 0.1, 'EdgeColor', 'r', 'LineWidth', 1, 'Tag', 'peakRegion');

                % Delete temporary patches
                delete(findobj(ax, 'Tag', 'tempPatch'));
            else
                % Incomplete selection - clean up
                delete(tempLines(ishandle(tempLines)));
                delete(tempTexts(ishandle(tempTexts)));
                delete(findobj(ax, 'Tag', 'tempPatch'));
            end

        catch ME
            if strcmp(ME.identifier, 'MATLAB:ginput:Interrupted')
                selecting = false;
            else
                % Clean up on error
                delete(tempLines(ishandle(tempLines)));
                delete(tempTexts(ishandle(tempTexts)));
                rethrow(ME);
            end
        end

        if ~selecting && peakCount == 0
            break;
        end
    end

    % Clean up UI elements
    if showGuide
        if exist('instructionBox', 'var') && ishandle(instructionBox)
            delete(instructionBox);
        end
        if exist('helpBox', 'var') && ishandle(helpBox)
            delete(helpBox);
        end
        if exist('progressBox', 'var') && ishandle(progressBox)
            delete(progressBox);
        end
    end

    % Clean up any remaining temporary elements
    delete(findobj(ax, 'Tag', 'tempPatch'));

    hold(ax, 'off');
end

function rangedCount = calculateRangedAtoms(mc, counts, rangeTable, binWidth)
    % Calculate total counts within all defined ranges

    rangedCount = 0;

    if ~ismember('mcbegin', rangeTable.Properties.VariableNames) || ...
       ~ismember('mcend', rangeTable.Properties.VariableNames)
        return;
    end

    for i = 1:height(rangeTable)
        mcBegin = rangeTable.mcbegin(i);
        mcEnd = rangeTable.mcend(i);

        inRange = mc >= mcBegin & mc <= mcEnd;
        rangedCount = rangedCount + sum(counts(inRange));
    end
end

function [bgFit, bgCounts, fitCurve] = fitBackground(mc, counts, inBg, inPeak, method)
    % Fit background and estimate counts under peak

    mcBg = mc(inBg);
    countsBg = counts(inBg);
    mcPeak = mc(inPeak);

    switch method
        case 'linear'
            % Linear least squares fit
            if length(mcBg) >= 2
                p = polyfit(mcBg, countsBg, 1);
                bgFit = struct('type', 'linear', 'slope', p(1), 'intercept', p(2));
                fitCurve = polyval(p, mc);
                bgCounts = sum(polyval(p, mcPeak));
            else
                bgFit = struct('type', 'constant', 'value', mean(countsBg));
                fitCurve = ones(size(mc)) * mean(countsBg);
                bgCounts = mean(countsBg) * sum(inPeak);
            end

        case 'exponential'
            % Exponential decay fit: y = a * exp(b * x)
            if length(mcBg) >= 2 && all(countsBg > 0)
                try
                    logCounts = log(countsBg);
                    p = polyfit(mcBg, logCounts, 1);
                    a = exp(p(2));
                    b = p(1);
                    bgFit = struct('type', 'exponential', 'a', a, 'b', b);
                    fitCurve = a * exp(b * mc);
                    bgCounts = sum(a * exp(b * mcPeak));
                catch
                    % Fall back to linear
                    p = polyfit(mcBg, countsBg, 1);
                    bgFit = struct('type', 'linear', 'slope', p(1), 'intercept', p(2));
                    fitCurve = polyval(p, mc);
                    bgCounts = sum(polyval(p, mcPeak));
                end
            else
                bgFit = struct('type', 'constant', 'value', mean(countsBg));
                fitCurve = ones(size(mc)) * mean(countsBg);
                bgCounts = mean(countsBg) * sum(inPeak);
            end

        case 'constant'
            % Constant background (average of background regions)
            avgBg = mean(countsBg);
            bgFit = struct('type', 'constant', 'value', avgBg);
            fitCurve = ones(size(mc)) * avgBg;
            bgCounts = avgBg * sum(inPeak);

        otherwise
            error('peakAnalysis:invalidMethod', 'Unknown background method: %s', method);
    end

    % Ensure non-negative background
    bgCounts = max(0, bgCounts);
    fitCurve = max(0, fitCurve);
end

function plotResults(results, isNormalized, hasRangedAtoms)
    % Plot peak analysis results

    nPeaks = length(results);

    % Create figure with subplots
    fig = figure('Name', 'Peak Analysis Results', 'Color', 'w');

    if nPeaks == 1
        nRows = 1;
        nCols = 1;
    elseif nPeaks <= 2
        nRows = 1;
        nCols = 2;
    elseif nPeaks <= 4
        nRows = 2;
        nCols = 2;
    else
        nRows = ceil(nPeaks / 3);
        nCols = 3;
    end

    for p = 1:nPeaks
        res = results(p);

        subplot(nRows, nCols, p);

        % Plot spectrum
        area(res.mcRegion, res.countsRegion, 'FaceColor', [0.9 0.9 0.9], ...
            'DisplayName', 'Spectrum');
        hold on;

        % Highlight background regions
        bgMask = res.inBg;
        scatter(res.mcRegion(bgMask), res.countsRegion(bgMask), 10, 'b', 'filled', ...
            'DisplayName', 'Background');

        % Highlight peak region
        peakMask = res.inPeak;
        area(res.mcRegion(peakMask), res.countsRegion(peakMask), 'FaceColor', [1 0.8 0.8], ...
            'DisplayName', 'Peak');

        % Plot background fit
        plot(res.mcRegion, res.fitCurve, 'r-', 'LineWidth', 2, 'DisplayName', 'BG Fit');

        % Mark peak location
        yLim = get(gca, 'YLim');
        plot([res.peakLocation res.peakLocation], yLim, 'k--', 'LineWidth', 1, ...
            'DisplayName', 'Peak Max');

        hold off;

        set(gca, 'YScale', 'log');
        xlabel('m/z [Da]');
        if isNormalized
            ylabel('Normalized frequency');
        else
            ylabel('Counts');
        end

        % Add text annotation with results
        title(sprintf('Peak @ %.2f Da', res.peakLocation));

        % Build annotation text
        if res.ppm >= 100
            concStr = sprintf('%.2f %%', res.fractionTotal * 100);
        else
            concStr = sprintf('%.1f ppm', res.ppm);
        end

        txt = {sprintf('Counts: %.0f', res.correctedCounts), ...
               sprintf('Fraction: %s', concStr), ...
               sprintf('S/B: %.1f', res.signalToBg)};

        if hasRangedAtoms && ~isnan(res.fractionRanged)
            if res.fractionRanged * 100 >= 0.1
                txt{end+1} = sprintf('Of ranged: %.2f %%', res.fractionRanged * 100);
            else
                txt{end+1} = sprintf('Of ranged: %.1f ppm', res.fractionRanged * 1e6);
            end
        end

        xLim = get(gca, 'XLim');
        yLim = get(gca, 'YLim');
        text(xLim(1) + 0.02*(xLim(2)-xLim(1)), ...
             10^(log10(yLim(1)) + 0.85*(log10(yLim(2))-log10(yLim(1)))), ...
             txt, 'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k');

        legend('Location', 'northeast', 'FontSize', 7);
    end

    % Add summary table if multiple peaks
    if nPeaks > 1
        % Create summary in figure title
        sgtitle(sprintf('Peak Analysis: %d peaks analyzed', nPeaks));
    end
end
