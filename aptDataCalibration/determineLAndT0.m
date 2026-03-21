function [L, t0, fitInfo] = determineLAndT0(pos, selectedPeaks, options)
% DETERMINELANDT0 Determine flight path length L and propagation delay t0.
%
% [L, t0, fitInfo] = determineLAndT0(pos, selectedPeaks)
% [L, t0, fitInfo] = determineLAndT0(pos, selectedPeaks, 'showPlot', true)
%
% Performs linear regression on the relationship:
%   tof = t0 + L * sqrt(idealMc * amu / (2 * e * VDC)) * 1e6
%
% INPUT
%   pos           - table with tof (ns), VDC (V), detx (mm), dety (mm)
%   selectedPeaks - table with columns: mcLow, mcHigh, idealMc (Da)
%
% OPTIONS
%   'detectorWindow'     - restrict to centre +/- this (mm, default 5)
%   'maxSamplePerPeak'   - max ions per peak (default 5000)
%   'flightPathLength'   - if given (mm), fix L and only fit t0 (default NaN)
%   'showPlot'           - show regression plot (default true)
%
% OUTPUT
%   L       - fitted flight path length (mm)
%   t0      - fitted propagation delay (ns)
%   fitInfo - struct with: seb_factor, seb_t, rSquared, residuals
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    selectedPeaks table
    options.detectorWindow (1,1) double = 5
    options.maxSamplePerPeak (1,1) double = 5000
    options.flightPathLength (1,1) double = NaN
    options.showPlot (1,1) logical = true
end

e   = 1.602176634e-19;
amu = 1.66053906660e-27;

sebFactorAll = [];
sebTAll      = [];

nPeaks = height(selectedPeaks);
for k = 1:nPeaks
    lo      = selectedPeaks.mcLow(k);
    hi      = selectedPeaks.mcHigh(k);
    idealMc = selectedPeaks.idealMc(k);

    % Mask: ions in peak range AND near detector centre
    mask = (pos.mc >= lo) & (pos.mc <= hi) & ...
           (abs(pos.detx) < options.detectorWindow) & ...
           (abs(pos.dety) < options.detectorWindow);

    idxK = find(mask);
    if isempty(idxK)
        fprintf('Peak %.2f-%.2f Da: no ions in centre window — skipping.\n', lo, hi);
        continue;
    end

    % Subsample if too many
    if numel(idxK) > options.maxSamplePerPeak
        idxK = idxK(randperm(numel(idxK), options.maxSamplePerPeak));
    end

    % Regression feature: sqrt(m * amu / (2 * e * V)) scaled to ns/mm
    sf = sqrt(idealMc * amu ./ (2 * e .* pos.VDC(idxK))) * 1e6;
    tK = pos.tof(idxK);

    sebFactorAll = [sebFactorAll; sf(:)];  %#ok<AGROW>
    sebTAll      = [sebTAll;      tK(:)];  %#ok<AGROW>

    fprintf('Peak %.2f–%.2f Da (idealMc=%.4f): %d ions used\n', ...
            lo, hi, idealMc, numel(idxK));
end

if numel(sebFactorAll) < 10
    error('determineLAndT0:tooFewPoints', ...
        'Only %d data points — need at least 10 for regression.', numel(sebFactorAll));
end

% Regression
if isnan(options.flightPathLength)
    % Fit both L and t0: tof = t0 + L * factor
    pFit = polyfit(sebFactorAll, sebTAll, 1);
    L  = pFit(1);   % slope
    t0 = pFit(2);   % intercept
else
    % Fix L, only fit t0
    L  = options.flightPathLength;
    t0 = mean(sebTAll) - L * mean(sebFactorAll);
    pFit = [L, t0];
end

% Goodness of fit
predicted = polyval(pFit, sebFactorAll);
residuals = sebTAll - predicted;
ssRes = sum(residuals.^2);
ssTot = sum((sebTAll - mean(sebTAll)).^2);
rSquared = 1 - ssRes / ssTot;

fitInfo = struct('seb_factor', sebFactorAll, 'seb_t', sebTAll, ...
                 'rSquared', rSquared, 'residuals', residuals, 'poly', pFit);

fprintf('\n========================================\n');
fprintf('  Fitted flight path L = %.2f mm\n', L);
fprintf('  Fitted delay      t0 = %.2f ns\n', t0);
fprintf('  R-squared            = %.6f\n', rSquared);
fprintf('========================================\n');

% Plot
if options.showPlot
    figure('Name', 'L and t0 Regression', 'NumberTitle', 'off');
    scatter(sebFactorAll, sebTAll, 2, 'b', 'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    xFit = linspace(min(sebFactorAll), max(sebFactorAll), 200);
    plot(xFit, polyval(pFit, xFit), 'r-', 'LineWidth', 2);

    if ~isnan(options.flightPathLength)
        % Also show the free-fit for comparison
        pFree = polyfit(sebFactorAll, sebTAll, 1);
        plot(xFit, polyval(pFree, xFit), 'g--', 'LineWidth', 1.5);
        legend('Data', ...
            sprintf('Fixed L=%g: t0=%.1f', options.flightPathLength, t0), ...
            sprintf('Free fit: L=%.1f, t0=%.1f', pFree(1), pFree(2)), ...
            'Location', 'northwest');
    else
        legend('Data', sprintf('Fit: L=%.1f, t0=%.1f', L, t0), ...
            'Location', 'northwest');
    end

    xlabel('\surd(m_{ideal} \cdot amu / 2eV) \times 10^6');
    ylabel('Time-of-flight (ns)');
    title(sprintf('L = %.1f mm,  t_0 = %.1f ns  (R^2 = %.4f)', L, t0, rSquared));
    hold off;
end

end
