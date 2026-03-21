function [pos, fitResult] = bowlCorrection(pos, peakRange, options)
% BOWLCORRECTION Apply spatial (bowl) mass-to-charge correction.
%
% [pos, fitResult] = bowlCorrection(pos, peakRange)
% [pos, fitResult] = bowlCorrection(pos, peakRange, 'gridSize', 3)
%
% Fits a 2-D quadratic surface f(x,y) = a + bx + cy + dx^2 + exy + fy^2
% to the normalised peak position across the detector, then corrects:
% mc = mc ./ f(x,y).
%
% INPUT
%   pos       - table with mc (Da), detx (mm), dety (mm) columns
%   peakRange - [lo hi] Da range of the reference peak
%
% OPTIONS
%   'gridSize'     - detector grid cell size in mm (default 5)
%   'binWidth'     - histogram bin width in Da (default 0.01)
%   'sampleMethod' - 'histogram', 'mean', or 'median' (default 'histogram')
%   'showPlot'     - show diagnostic plot (default true)
%
% OUTPUT
%   pos       - table with corrected mc column
%   fitResult - struct with field 'params' (6-element vector [a b c d e f])
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    peakRange (1,2) double
    options.gridSize (1,1) double = 5
    options.binWidth (1,1) double = 0.01
    options.sampleMethod (1,:) char {mustBeMember(options.sampleMethod, {'histogram','mean','median'})} = 'histogram'
    options.showPlot (1,1) logical = true
end

mc   = pos.mc;
detx = pos.detx;   % already in mm
dety = pos.dety;

% 1. Mask reference peak
peakMask = (mc >= peakRange(1)) & (mc <= peakRange(2));
mcPeak   = mc(peakMask);
xPeak    = detx(peakMask);
yPeak    = dety(peakMask);

if sum(peakMask) < 100
    warning('bowlCorrection:tooFewIons', ...
        'Only %d ions in peak range. Skipping correction.', sum(peakMask));
    fitResult = [];
    return;
end

% 2. Reference peak location
refLoc = findPeakLocation(mcPeak, options.binWidth, options.sampleMethod);
if refLoc == 0
    warning('bowlCorrection:noReference', 'Could not locate reference peak.');
    fitResult = [];
    return;
end

% 3. Grid the detector
d = options.gridSize;
xLo = floor(min(xPeak));  xHi = ceil(max(xPeak));
yLo = floor(min(yPeak));  yHi = ceil(max(yPeak));

xs = [];  ys = [];  zs = [];

for yi = yLo:d:(yHi - d)
    for xi = xLo:d:(xHi - d)
        cellMask = (xPeak >= xi) & (xPeak < xi + d) & ...
                   (yPeak >= yi) & (yPeak < yi + d);
        if sum(cellMask) < 10, continue; end

        loc = findPeakLocation(mcPeak(cellMask), options.binWidth, options.sampleMethod);
        if loc == 0, continue; end

        xs(end+1) = median(xPeak(cellMask)); %#ok<AGROW>
        ys(end+1) = median(yPeak(cellMask)); %#ok<AGROW>
        zs(end+1) = loc / refLoc;             %#ok<AGROW>
    end
end

xs = xs(:);  ys = ys(:);  zs = zs(:);

if numel(xs) < 6
    warning('bowlCorrection:tooFewCells', ...
        'Only %d grid cells — need at least 6 for 2-D quadratic.', numel(xs));
    fitResult = [];
    return;
end

% 4. Fit: f(x,y) = a + bx + cy + dx^2 + exy + fy^2
A = [ones(numel(xs),1), xs, ys, xs.^2, xs.*ys, ys.^2];
params = A \ zs;
fitResult = struct('params', params, 'xs', xs, 'ys', ys, 'zs', zs);

% 5. Apply correction to all ions
Aall = [ones(numel(detx),1), detx, dety, detx.^2, detx.*dety, dety.^2];
fBowl = Aall * params;
pos.mc = pos.mc ./ fBowl;

% 6. Diagnostic plot
if options.showPlot
    figure('Name', 'Bowl Correction', 'NumberTitle', 'off');

    subplot(1,2,1);
    scatter(xs, ys, 30, zs, 'filled');
    colorbar; axis equal;
    xlabel('detx (mm)'); ylabel('dety (mm)');
    title('Sampled peak position (normalised)');

    subplot(1,2,2);
    xg = linspace(min(xs), max(xs), 40);
    yg = linspace(min(ys), max(ys), 40);
    [Xg, Yg] = meshgrid(xg, yg);
    Zg = params(1) + params(2)*Xg + params(3)*Yg + ...
         params(4)*Xg.^2 + params(5)*Xg.*Yg + params(6)*Yg.^2;
    surf(Xg, Yg, Zg, 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on;
    scatter3(xs, ys, zs, 30, 'r', 'filled');
    xlabel('detx (mm)'); ylabel('dety (mm)'); zlabel('Normalised peak');
    title('Quadratic surface fit');
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
