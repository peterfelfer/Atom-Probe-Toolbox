function pos = plotFieldDesorptionMap(pos, options)
% PLOTFIELDDESORPTIONMAP Show FDM and crop with a circular ROI.
%
% pos = plotFieldDesorptionMap(pos)
%
% Displays a 2-D histogram of detector (detx, dety) positions.
% The user draws a circle to select the region of interest.
% Ions outside the circle are removed.
%
% INPUT
%   pos - table with detx (mm) and dety (mm) columns
%
% OPTIONS
%   'numBins' - histogram bins per axis (default 256)
%
% OUTPUT
%   pos - cropped table with ionIdx recomputed
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.numBins (1,1) double = 256
end

x = pos.detx;
y = pos.dety;

xEdges = linspace(min(x), max(x), options.numBins + 1);
yEdges = linspace(min(y), max(y), options.numBins + 1);
counts = histcounts2(x, y, xEdges, yEdges);
counts(counts == 0) = NaN;

fig = figure('Name', 'Field Desorption Map — draw circle to crop', ...
             'NumberTitle', 'off');
ax = axes(fig);
imagesc(ax, xEdges(1:end-1), yEdges(1:end-1), log10(counts'));
set(ax, 'YDir', 'normal');
axis(ax, 'equal', 'tight');
xlabel(ax, 'detx (mm)');
ylabel(ax, 'dety (mm)');
colormap(ax, hot);
cb = colorbar(ax);
cb.Label.String = 'log_{10}(counts)';
title(ax, 'Draw circle to select ROI, then double-click inside it');

% Interactive selection
try
    roi = drawcircle(ax, 'Color', 'c', 'LineWidth', 2);
    l = addlistener(roi, 'ROIClicked', @(~,evt) roiCallback(evt));
    uiwait(fig);
    delete(l);
    cx = roi.Center(1);
    cy = roi.Center(2);
    r  = roi.Radius;
catch
    fprintf('Click centre, then edge of ROI.\n');
    [gx, gy] = ginput(2);
    cx = gx(1);  cy = gy(1);
    r  = sqrt((gx(2)-gx(1))^2 + (gy(2)-gy(1))^2);
end

nBefore = height(pos);
dist = sqrt((x - cx).^2 + (y - cy).^2);
pos  = pos(dist <= r, :);
pos.ionIdx = (1:height(pos))';

fprintf('Spatial crop: kept %d of %d ions (circle r=%.1f mm at [%.1f, %.1f])\n', ...
        height(pos), nBefore, r, cx, cy);
close(fig);

end


function roiCallback(evt)
    if strcmp(evt.SelectionType, 'double')
        uiresume(ancestor(evt.Source.Parent, 'figure'));
    end
end
