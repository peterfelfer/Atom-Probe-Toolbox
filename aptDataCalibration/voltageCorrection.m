function [pos, fitResult] = voltageCorrection(pos, peakRange, options)
% VOLTAGECORRECTION Apply voltage-dependent mass-to-charge correction.
%
% [pos, fitResult] = voltageCorrection(pos, peakRange)
% [pos, fitResult] = voltageCorrection(pos, peakRange, 'sampleSize', 2000)
%
% Fits a quadratic model f(V) = a + bV + cV^2 to the normalised peak
% position vs. voltage, then corrects: mc = mc ./ sqrt(f(V)).
%
% INPUT
%   pos       - table with mc (Da) and VDC (V) columns
%   peakRange - [lo hi] Da range of the reference peak
%
% OPTIONS
%   'sampleSize'   - ions per group (default 1000)
%   'mode'         - 'ionSequence' or 'voltage' grouping (default 'ionSequence')
%   'binWidth'     - histogram bin width in Da (default 0.01)
%   'sampleMethod' - 'histogram', 'mean', or 'median' (default 'histogram')
%   'showPlot'     - show diagnostic plot (default true)
%
% OUTPUT
%   pos       - table with corrected mc column
%   fitResult - struct with fields: poly, V_means, mc_norms
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    peakRange (1,2) double
    options.sampleSize (1,1) double = 1000
    options.mode (1,:) char {mustBeMember(options.mode, {'ionSequence','voltage'})} = 'ionSequence'
    options.binWidth (1,1) double = 0.01
    options.sampleMethod (1,:) char {mustBeMember(options.sampleMethod, {'histogram','mean','median'})} = 'histogram'
    options.showPlot (1,1) logical = true
end

mc = pos.mc;
hv = pos.VDC;

% 1. Mask reference peak
peakMask = (mc >= peakRange(1)) & (mc <= peakRange(2));
mcPeak   = mc(peakMask);
hvPeak   = hv(peakMask);

if sum(peakMask) < 100
    warning('voltageCorrection:tooFewIons', ...
        'Only %d ions in peak range. Skipping correction.', sum(peakMask));
    fitResult = [];
    return;
end

% 2. Reference peak location
refLoc = findPeakLocation(mcPeak, options.binWidth, options.sampleMethod);
if refLoc == 0
    warning('voltageCorrection:noReference', 'Could not locate reference peak.');
    fitResult = [];
    return;
end

% 3. Divide into groups, find normalised peak position per group
vMeans  = [];
mcNorms = [];

switch options.mode
    case 'ionSequence'
        nGroups = floor(numel(mcPeak) / options.sampleSize);
        for g = 1:nGroups
            idx = ((g-1)*options.sampleSize + 1):(g*options.sampleSize);
            loc = findPeakLocation(mcPeak(idx), options.binWidth, options.sampleMethod);
            if loc > 0
                vMeans(end+1)  = mean(hvPeak(idx)); %#ok<AGROW>
                mcNorms(end+1) = loc / refLoc;       %#ok<AGROW>
            end
        end

    case 'voltage'
        vEdges = min(hvPeak):options.sampleSize:max(hvPeak);
        for g = 1:numel(vEdges)-1
            maskG = (hvPeak >= vEdges(g)) & (hvPeak < vEdges(g+1));
            if sum(maskG) < 10, continue; end
            loc = findPeakLocation(mcPeak(maskG), options.binWidth, options.sampleMethod);
            if loc > 0
                vMeans(end+1)  = mean(hvPeak(maskG)); %#ok<AGROW>
                mcNorms(end+1) = loc / refLoc;          %#ok<AGROW>
            end
        end
end

if numel(vMeans) < 3
    warning('voltageCorrection:tooFewGroups', ...
        'Only %d groups — need at least 3 for quadratic fit.', numel(vMeans));
    fitResult = [];
    return;
end

% 4. Quadratic fit: polyfit returns [c, b, a] for c*x^2 + b*x + a
pCoeff = polyfit(vMeans, mcNorms, 2);
fitResult = struct('poly', pCoeff, 'V_means', vMeans, 'mc_norms', mcNorms);

% 5. Apply correction
fV = polyval(pCoeff, hv);
pos.mc = pos.mc ./ sqrt(fV);

% 6. Diagnostic plot
if options.showPlot
    figure('Name', 'Voltage Correction', 'NumberTitle', 'off');
    plot(vMeans, mcNorms, 'bo', 'MarkerSize', 4); hold on;
    vFit = linspace(min(vMeans), max(vMeans), 200);
    plot(vFit, polyval(pCoeff, vFit), 'r-', 'LineWidth', 1.5);
    xlabel('DC Voltage (V)');
    ylabel('Normalised peak position');
    title('Voltage correction: f(V) = a + bV + cV^2');
    legend('Data', 'Fit', 'Location', 'best');
    hold off;
end

end


function loc = findPeakLocation(mcSubset, binWidth, method)
    loc = 0;
    if numel(mcSubset) < 5, return; end
    switch method
        case 'histogram'
            edges = min(mcSubset):binWidth:max(mcSubset);
            if numel(edges) < 3, loc = mean(mcSubset); return; end
            counts  = histcounts(mcSubset, edges);
            centers = edges(1:end-1) + binWidth/2;
            [~, idx] = max(counts);
            loc = centers(idx);
        case 'mean'
            loc = mean(mcSubset);
        case 'median'
            loc = median(mcSubset);
    end
end
