function results = massSpecAnalysis(mass, options)
% MASSSPECANALYSIS Comprehensive mass spectrum analysis for APT data
%
% results = massSpecAnalysis(mass)
% results = massSpecAnalysis(mass, 'binWidth', 0.01)
%
% Performs automated analysis of mass-to-charge spectrum including peak
% detection, background fitting, mass resolution estimation, and optional
% peak deconvolution.
%
% INPUT:
%   mass - Vector of mass-to-charge values (Da)
%
% OPTIONS:
%   'binWidth'       - Histogram bin width in Da (default: 0.01)
%   'massRange'      - [min, max] mass range to analyze (default: auto)
%   'peakThreshold'  - Minimum peak height as fraction of max (default: 0.001)
%   'bgMethod'       - Background method: 'moving', 'polynomial', 'exponential'
%                      (default: 'moving')
%   'bgWindow'       - Window size for moving background (default: 100 bins)
%   'bgPolyOrder'    - Polynomial order for background (default: 3)
%   'fitPeaks'       - Fit Gaussian to peaks (default: true)
%   'deconvolve'     - Attempt peak deconvolution (default: false)
%   'showPlot'       - Display results plot (default: true)
%
% OUTPUT:
%   results - Structure containing:
%       .spectrum - Spectrum data
%           .edges      - Bin edges
%           .centers    - Bin centers (mass-to-charge)
%           .counts     - Raw counts per bin
%           .countRate  - Counts per Da
%       .background - Background estimation
%           .values     - Background at each bin
%           .method     - Method used
%           .corrected  - Background-corrected counts
%       .peaks - Detected peaks table with columns:
%           .position   - Peak center (Da)
%           .height     - Peak height (counts)
%           .area       - Integrated peak area
%           .fwhm       - Full width at half maximum
%           .resolution - Mass resolution (M/dM)
%           .snr        - Signal-to-noise ratio
%       .statistics - Overall spectrum statistics
%           .totalCounts    - Total counts
%           .meanResolution - Mean mass resolution
%           .noiseLevel     - Estimated noise level
%
% EXAMPLES:
%   % Basic analysis
%   results = massSpecAnalysis(posTable.mc);
%
%   % With specific settings
%   results = massSpecAnalysis(mc, 'binWidth', 0.005, 'bgMethod', 'polynomial');
%
%   % Analyze specific range
%   results = massSpecAnalysis(mc, 'massRange', [20, 80], 'showPlot', true);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    mass (:,1) double
    options.binWidth (1,1) double {mustBePositive} = 0.01
    options.massRange (1,2) double = []
    options.peakThreshold (1,1) double {mustBePositive} = 0.001
    options.bgMethod (1,:) char {mustBeMember(options.bgMethod, {'moving', 'polynomial', 'exponential'})} = 'moving'
    options.bgWindow (1,1) double {mustBePositive} = 100
    options.bgPolyOrder (1,1) double {mustBePositive, mustBeInteger} = 3
    options.fitPeaks (1,1) logical = true
    options.deconvolve (1,1) logical = false
    options.showPlot (1,1) logical = true
end

% Determine mass range
if isempty(options.massRange)
    options.massRange = [0, ceil(max(mass) * 1.05)];
end

mass = mass(mass >= options.massRange(1) & mass <= options.massRange(2));

%% Create Histogram
fprintf('Creating mass spectrum histogram...\n');

edges = options.massRange(1):options.binWidth:options.massRange(2);
counts = histcounts(mass, edges);
centers = edges(1:end-1) + options.binWidth/2;

results = struct();
results.spectrum.edges = edges;
results.spectrum.centers = centers;
results.spectrum.counts = counts;
results.spectrum.countRate = counts / options.binWidth;
results.spectrum.binWidth = options.binWidth;

%% Background Estimation
fprintf('Estimating background...\n');

[background, bgMethod] = estimateBackground(counts, options);

results.background.values = background;
results.background.method = bgMethod;
results.background.corrected = max(counts - background, 0);

%% Peak Detection
fprintf('Detecting peaks...\n');

peakTable = detectPeaks(centers, counts, background, options);

results.peaks = peakTable;

%% Peak Fitting
if options.fitPeaks && ~isempty(peakTable)
    fprintf('Fitting peaks...\n');
    results.peaks = fitPeaksGaussian(centers, counts, background, peakTable, options);
end

%% Statistics
fprintf('Computing statistics...\n');

results.statistics.totalCounts = sum(counts);
results.statistics.massRange = options.massRange;
results.statistics.nPeaks = height(results.peaks);

if ~isempty(results.peaks)
    results.statistics.meanResolution = mean(results.peaks.resolution(~isnan(results.peaks.resolution)));
    results.statistics.medianSNR = median(results.peaks.snr);
else
    results.statistics.meanResolution = NaN;
    results.statistics.medianSNR = NaN;
end

% Noise level estimation
noiseRegion = counts(counts < median(counts));
results.statistics.noiseLevel = std(noiseRegion);
results.statistics.backgroundFraction = sum(background) / sum(counts);

%% Optional Deconvolution
if options.deconvolve && ~isempty(results.peaks)
    fprintf('Attempting peak deconvolution...\n');
    results.deconvolution = deconvolvePeaks(centers, counts, background, results.peaks);
end

%% Plot Results
if options.showPlot
    plotMassSpectrum(results, options);
end

fprintf('Mass spectrum analysis complete: %d peaks detected\n', results.statistics.nPeaks);

end

%% Background Estimation
function [background, method] = estimateBackground(counts, options)
    nBins = length(counts);

    switch lower(options.bgMethod)
        case 'moving'
            % Moving minimum (robust against peaks)
            windowSize = min(options.bgWindow, nBins);
            background = movmin(counts, windowSize);

            % Smooth the result
            background = movmean(background, windowSize/2);
            method = sprintf('moving_min (window=%d)', windowSize);

        case 'polynomial'
            % Polynomial fit to local minima
            x = 1:nBins;

            % Find local minima for fitting
            localMin = islocalmin(counts, 'MinProminence', std(counts)*0.5);

            % Add edges
            localMin(1:min(10, nBins)) = true;
            localMin(max(1, nBins-10):nBins) = true;

            if sum(localMin) > options.bgPolyOrder + 1
                p = polyfit(x(localMin), counts(localMin), options.bgPolyOrder);
                background = polyval(p, x);
                background = max(background, 0);
            else
                % Fall back to moving minimum
                background = movmin(counts, options.bgWindow);
            end
            method = sprintf('polynomial (order=%d)', options.bgPolyOrder);

        case 'exponential'
            % Exponential decay fit (common for APT thermal tails)
            x = (1:nBins)';

            % Fit to local minima
            localMin = islocalmin(counts, 'MinProminence', std(counts)*0.5);
            localMin(1:min(10, nBins)) = true;

            if sum(localMin) > 3
                try
                    % Fit: bg = a * exp(-b*x) + c
                    minCounts = counts(localMin);
                    minX = x(localMin);

                    fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                                   'Lower', [0, 0, 0], ...
                                   'Upper', [max(counts), 1, max(counts)], ...
                                   'StartPoint', [max(counts)/10, 0.01, min(counts)]);
                    ft = fittype('a*exp(-b*x) + c', 'options', fo);
                    fitResult = fit(minX, minCounts, ft);

                    background = fitResult.a * exp(-fitResult.b * x) + fitResult.c;
                    background = background';
                catch
                    % Fall back to moving minimum
                    background = movmin(counts, options.bgWindow);
                end
            else
                background = movmin(counts, options.bgWindow);
            end
            method = 'exponential';
    end

    % Ensure background doesn't exceed counts
    background = min(background, counts);
end

%% Peak Detection
function peakTable = detectPeaks(centers, counts, background, options)
    % Detect peaks above threshold

    corrected = counts - background;
    threshold = options.peakThreshold * max(corrected);

    % Find local maxima
    [pks, locs, widths, proms] = findpeaks(corrected, ...
        'MinPeakHeight', threshold, ...
        'MinPeakDistance', 3, ...  % Minimum 3 bins apart
        'MinPeakProminence', threshold * 0.5);

    nPeaks = length(pks);

    if nPeaks == 0
        peakTable = table();
        return;
    end

    % Calculate properties for each peak
    position = centers(locs)';
    height = pks';
    prominence = proms';

    % FWHM and area estimates
    fwhm = zeros(nPeaks, 1);
    area = zeros(nPeaks, 1);
    resolution = zeros(nPeaks, 1);
    snr = zeros(nPeaks, 1);

    binWidth = centers(2) - centers(1);

    for i = 1:nPeaks
        % Find FWHM
        halfMax = pks(i) / 2;
        peakIdx = locs(i);

        % Search left
        leftIdx = peakIdx;
        while leftIdx > 1 && corrected(leftIdx) > halfMax
            leftIdx = leftIdx - 1;
        end

        % Search right
        rightIdx = peakIdx;
        while rightIdx < length(corrected) && corrected(rightIdx) > halfMax
            rightIdx = rightIdx + 1;
        end

        % Interpolate for better precision
        if leftIdx > 1 && rightIdx < length(corrected)
            leftPos = interp1(corrected(leftIdx:leftIdx+1), centers(leftIdx:leftIdx+1), halfMax);
            rightPos = interp1(corrected(rightIdx-1:rightIdx), centers(rightIdx-1:rightIdx), halfMax);
            fwhm(i) = rightPos - leftPos;
        else
            fwhm(i) = widths(i) * binWidth;
        end

        % Resolution
        if fwhm(i) > 0
            resolution(i) = position(i) / fwhm(i);
        else
            resolution(i) = NaN;
        end

        % Area (simple sum)
        peakRange = max(1, leftIdx):min(length(corrected), rightIdx);
        area(i) = sum(corrected(peakRange)) * binWidth;

        % SNR
        localBg = mean(background(peakRange));
        localNoise = sqrt(localBg);  % Poisson noise
        if localNoise > 0
            snr(i) = pks(i) / localNoise;
        else
            snr(i) = pks(i);
        end
    end

    peakTable = table(position, height, area, fwhm, resolution, snr, prominence);
    peakTable = sortrows(peakTable, 'position');
end

%% Gaussian Peak Fitting
function peakTable = fitPeaksGaussian(centers, counts, background, peakTable, options)
    % Fit Gaussian to each peak for better parameters

    if isempty(peakTable)
        return;
    end

    nPeaks = height(peakTable);
    corrected = counts - background;
    binWidth = centers(2) - centers(1);

    % Add columns for fit parameters
    peakTable.fitCenter = peakTable.position;
    peakTable.fitHeight = peakTable.height;
    peakTable.fitSigma = peakTable.fwhm / 2.355;  % FWHM = 2.355 * sigma
    peakTable.fitR2 = zeros(nPeaks, 1);

    for i = 1:nPeaks
        pos = peakTable.position(i);
        sigma_est = peakTable.fwhm(i) / 2.355;

        % Define fitting region (3 sigma each side)
        fitRange = abs(centers - pos) < 5 * sigma_est;

        if sum(fitRange) < 5
            continue;
        end

        xFit = centers(fitRange);
        yFit = corrected(fitRange);

        try
            % Gaussian: a * exp(-((x-b)^2)/(2*c^2))
            fo = fitoptions('Method', 'NonlinearLeastSquares', ...
                           'Lower', [0, pos-sigma_est, 0], ...
                           'Upper', [max(yFit)*2, pos+sigma_est, sigma_est*5], ...
                           'StartPoint', [peakTable.height(i), pos, sigma_est]);
            ft = fittype('a*exp(-((x-b)^2)/(2*c^2))', 'options', fo);

            [fitResult, gof] = fit(xFit', yFit', ft);

            peakTable.fitCenter(i) = fitResult.b;
            peakTable.fitHeight(i) = fitResult.a;
            peakTable.fitSigma(i) = fitResult.c;
            peakTable.fitR2(i) = gof.rsquare;

            % Update FWHM and resolution from fit
            peakTable.fwhm(i) = 2.355 * fitResult.c;
            if peakTable.fwhm(i) > 0
                peakTable.resolution(i) = fitResult.b / peakTable.fwhm(i);
            end
        catch
            % Keep original values if fit fails
            peakTable.fitR2(i) = NaN;
        end
    end
end

%% Peak Deconvolution
function deconv = deconvolvePeaks(centers, counts, background, peakTable)
    % Attempt to deconvolve overlapping peaks

    deconv = struct();
    deconv.attempted = true;

    % Find potentially overlapping peaks
    nPeaks = height(peakTable);
    overlapGroups = {};

    for i = 1:nPeaks
        pos_i = peakTable.position(i);
        fwhm_i = peakTable.fwhm(i);

        for j = i+1:nPeaks
            pos_j = peakTable.position(j);

            % Check if peaks overlap (within 2*FWHM)
            if abs(pos_i - pos_j) < (fwhm_i + peakTable.fwhm(j))
                overlapGroups{end+1} = [i, j];
            end
        end
    end

    deconv.overlapGroups = overlapGroups;
    deconv.nOverlapping = length(overlapGroups);

    % For each overlap group, fit multiple Gaussians
    deconv.fits = {};

    for g = 1:length(overlapGroups)
        group = overlapGroups{g};
        positions = peakTable.position(group);

        % Define fitting region
        minPos = min(positions) - 3*max(peakTable.fwhm(group));
        maxPos = max(positions) + 3*max(peakTable.fwhm(group));
        fitRange = centers >= minPos & centers <= maxPos;

        if sum(fitRange) < 10
            continue;
        end

        xFit = centers(fitRange);
        yFit = (counts(fitRange) - background(fitRange))';

        % Multi-Gaussian fit would go here
        % This is a placeholder for more sophisticated deconvolution
        deconv.fits{g} = struct('group', group, 'success', false, 'message', 'Not implemented');
    end
end

%% Plotting
function plotMassSpectrum(results, options)
    figure('Position', [100 100 1200 600], 'Name', 'Mass Spectrum Analysis');

    % Main spectrum plot
    subplot(2,2,[1,3]);
    semilogy(results.spectrum.centers, results.spectrum.counts + 1, 'b-', 'LineWidth', 0.5);
    hold on;
    semilogy(results.spectrum.centers, results.background.values + 1, 'r-', 'LineWidth', 1);

    % Mark detected peaks
    if ~isempty(results.peaks)
        for i = 1:height(results.peaks)
            xline(results.peaks.position(i), 'g--', 'Alpha', 0.5);
        end
    end

    xlabel('Mass-to-charge (Da)');
    ylabel('Counts (log scale)');
    title(sprintf('Mass Spectrum (%d peaks detected)', results.statistics.nPeaks));
    legend('Spectrum', 'Background', 'Location', 'northeast');
    grid on;

    % Peak list
    subplot(2,2,2);
    if ~isempty(results.peaks) && height(results.peaks) > 0
        % Show top peaks by height
        topPeaks = sortrows(results.peaks, 'height', 'descend');
        topPeaks = topPeaks(1:min(15, height(topPeaks)), :);

        barh(1:height(topPeaks), topPeaks.height);
        yticks(1:height(topPeaks));
        yticklabels(arrayfun(@(x) sprintf('%.2f', x), topPeaks.position, 'UniformOutput', false));
        xlabel('Peak Height (counts)');
        ylabel('Position (Da)');
        title('Top 15 Peaks');
        grid on;
    else
        text(0.5, 0.5, 'No peaks detected', 'HorizontalAlignment', 'center');
        axis off;
    end

    % Statistics
    subplot(2,2,4);
    axis off;

    stats_text = {
        sprintf('Total counts: %d', results.statistics.totalCounts)
        sprintf('Number of peaks: %d', results.statistics.nPeaks)
        sprintf('Mean resolution: %.0f', results.statistics.meanResolution)
        sprintf('Median SNR: %.1f', results.statistics.medianSNR)
        sprintf('Background fraction: %.1f%%', results.statistics.backgroundFraction * 100)
        sprintf('Mass range: %.1f - %.1f Da', results.statistics.massRange(1), results.statistics.massRange(2))
        sprintf('Bin width: %.3f Da', results.spectrum.binWidth)
        sprintf('Background method: %s', results.background.method)
    };

    text(0.1, 0.9, 'Statistics', 'FontSize', 12, 'FontWeight', 'bold');
    for i = 1:length(stats_text)
        text(0.1, 0.85 - i*0.1, stats_text{i}, 'FontSize', 10);
    end

    sgtitle('APT Mass Spectrum Analysis');
end
