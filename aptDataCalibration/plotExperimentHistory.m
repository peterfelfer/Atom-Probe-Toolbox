function pos = plotExperimentHistory(pos, options)
% PLOTEXPERIMENTHISTORY Show experiment history and crop temporally.
%
% pos = plotExperimentHistory(pos)
%
% Displays a 2-D histogram (ion index vs TOF) with DC voltage overlay.
% The user draws a rectangle to select a stable evaporation region.
% The dataset is cropped to that index range.
%
% INPUT
%   pos - table with tof (ns) and VDC (V) columns
%
% OPTIONS
%   'maxTof'   - y-axis limit in ns (default: auto from data)
%   'numBinsX' - histogram bins along ion index (default 500)
%   'numBinsY' - histogram bins along TOF (default 512)
%
% OUTPUT
%   pos - cropped table with ionIdx recomputed
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.maxTof (1,1) double = NaN
    options.numBinsX (1,1) double = 500
    options.numBinsY (1,1) double = 512
end

nIons = height(pos);
idx   = (1:nIons)';
maxT  = options.maxTof;
if isnan(maxT), maxT = max(pos.tof); end

% 2-D histogram
xEdges = linspace(1, nIons, min(options.numBinsX, nIons) + 1);
yEdges = linspace(0, maxT, options.numBinsY + 1);
counts = histcounts2(idx, pos.tof, xEdges, yEdges);
counts(counts == 0) = NaN;

fig = figure('Name', 'Experiment History — draw rectangle to crop', ...
             'NumberTitle', 'off');
ax = axes(fig);
imagesc(ax, xEdges(1:end-1), yEdges(1:end-1), log10(counts'));
set(ax, 'YDir', 'normal');
xlabel(ax, 'Ion index');
ylabel(ax, 'Time-of-flight (ns)');
colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = 'log_{10}(counts)';

% Voltage overlay
dec = max(1, floor(nIons / 5000));
yyaxis(ax, 'right');
plot(ax, idx(1:dec:end), pos.VDC(1:dec:end), 'r-', 'LineWidth', 1);
ylabel(ax, 'DC voltage (V)');
set(ax, 'YColor', 'r');
yyaxis(ax, 'left');

title(ax, 'Draw rectangle to crop, then double-click inside it');

% Interactive selection
try
    roi = drawrectangle(ax, 'Color', 'g', 'LineWidth', 2);
    l = addlistener(roi, 'ROIClicked', @(~,evt) roiCallback(evt));
    uiwait(fig);
    delete(l);
    p = roi.Position;
    x1 = max(1, round(p(1)));
    x2 = min(nIons, round(p(1) + p(3)));
catch
    fprintf('Click two points for crop range (left, right).\n');
    [gx, ~] = ginput(2);
    x1 = max(1, round(min(gx)));
    x2 = min(nIons, round(max(gx)));
end

fprintf('Temporal crop: ions %d–%d (%d of %d kept)\n', x1, x2, x2-x1+1, nIons);
pos = pos(x1:x2, :);
pos.ionIdx = (1:height(pos))';
close(fig);

end


function roiCallback(evt)
    if strcmp(evt.SelectionType, 'double')
        uiresume(ancestor(evt.Source.Parent, 'figure'));
    end
end
