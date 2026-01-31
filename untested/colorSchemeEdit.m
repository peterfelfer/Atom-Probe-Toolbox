function colorScheme = colorSchemeEdit(colorScheme, options)
% colorSchemeEdit opens an interactive HSV color wheel editor.
%
% colorScheme = colorSchemeEdit(colorScheme)
% colorScheme = colorSchemeEdit(colorScheme, options)
%
% INPUT
% colorScheme: table with columns ion and color
% options: struct with fields
%   title:      figure title (default: "Color Scheme Editor")
%   showLabels: initial label visibility (default: true)
%
% OUTPUT
% colorScheme: updated color scheme table
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    colorScheme table
    options.title (1,1) string = "Color Scheme Editor"
    options.showLabels (1,1) logical = true
end

if ~all(ismember({'ion','color'}, colorScheme.Properties.VariableNames))
    error('colorSchemeEdit:invalidColorScheme', ...
        'colorScheme must be a table with columns ion and color.');
end

if isempty(colorScheme)
    return;
end

rgb = colorScheme.color;
if size(rgb, 2) ~= 3
    error('colorSchemeEdit:invalidColorScheme', ...
        'colorScheme.color must be Nx3 RGB values.');
end

ionNames = string(colorScheme.ion);
hsv = rgb2hsv(rgb);

[x, y] = hsvToXY(hsv);

fig = figure('Name', char(options.title), ...
    'NumberTitle', 'off', ...
    'Color', 'w', ...
    'Units', 'pixels', ...
    'Position', [100 100 600 520]);
ax = axes(fig, 'Units', 'normalized', 'Position', [0.06 0.12 0.68 0.82]);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, [-1.05 1.05 -1.05 1.05]);
ax.XTick = [];
ax.YTick = [];
ax.Box = 'on';

drawColorWheel(ax);

scatterH = scatter(ax, x, y, 60, rgb, 'filled', ...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.5);

labelOffset = [0.02 0.02];
labels = gobjects(numel(ionNames), 1);
for i = 1:numel(ionNames)
    [labelColor, labelBg] = pickLabelColors(rgb(i, :));
    labels(i) = text(ax, x(i) + labelOffset(1), y(i) + labelOffset(2), ionNames(i), ...
        'FontSize', 9, 'Interpreter', 'none', 'Color', labelColor, ...
        'BackgroundColor', labelBg, 'Margin', 2);
end

labelToggle = uicontrol(fig, 'Style', 'checkbox', ...
    'String', 'Show labels', ...
    'Value', options.showLabels, ...
    'Units', 'normalized', ...
    'Position', [0.86 0.92 0.12 0.05], ...
    'Callback', @(src, ~) toggleLabels(src, labels));

doneBtn = uicontrol(fig, 'Style', 'pushbutton', ...
    'String', 'Done', ...
    'Units', 'normalized', ...
    'Position', [0.86 0.86 0.12 0.05], ...
    'Callback', @(~, ~) close(fig));

colorPanel = uipanel(fig, 'Units', 'normalized', ...
    'Position', [0.81 0.60 0.14 0.12], ...
    'BackgroundColor', rgb(1, :), ...
    'BorderType', 'line', 'HighlightColor', [0 0 0]);

valueLabel = uicontrol(fig, 'Style', 'text', ...
    'String', 'Value', ...
    'Units', 'normalized', ...
    'Position', [0.81 0.53 0.14 0.04], ...
    'BackgroundColor', 'w');

valueSlider = uicontrol(fig, 'Style', 'slider', ...
    'Min', 0, 'Max', 1, 'Value', hsv(1, 3), ...
    'Units', 'normalized', ...
    'Position', [0.865 0.18 0.03 0.30], ...
    'SliderStep', [0.01 0.1], ...
    'Callback', @(src, ~) updateValueFromSlider(fig, src));

toggleLabels(labelToggle, labels);

data = struct();
data.colorScheme = colorScheme;
data.hsv = hsv;
data.rgb = rgb;
data.x = x;
data.y = y;
data.scatter = scatterH;
data.labels = labels;
data.labelOffset = labelOffset;
data.dragIdx = [];
data.selectedIdx = 1;
data.colorPanel = colorPanel;
data.valueSlider = valueSlider;
data.valueLabel = valueLabel;
setappdata(fig, 'colorSchemeEditData', data);

fig.WindowButtonDownFcn = @(~, ~) startDrag(fig, ax);
fig.WindowButtonMotionFcn = @(~, ~) dragMove(fig, ax);
fig.WindowButtonUpFcn = @(~, ~) stopDrag(fig);
fig.CloseRequestFcn = @(~, ~) onClose(fig);

uiwait(fig);

data = getappdata(fig, 'colorSchemeEditData');
if isempty(data)
    return;
end
colorScheme = data.colorScheme;
delete(fig);
end

function drawColorWheel(ax)
    n = 300;
    lim = linspace(-1, 1, n);
    [X, Y] = meshgrid(lim, lim);
    [theta, r] = cart2pol(X, Y);
    hue = mod(theta / (2 * pi), 1);
    sat = min(r, 1);
    val = ones(size(hue));
    hsv = cat(3, hue, sat, val);
    rgb = hsv2rgb(hsv);
    alpha = r <= 1;
    image(ax, lim, lim, rgb, 'AlphaData', alpha);
    plot(ax, cos(linspace(0, 2*pi, 256)), sin(linspace(0, 2*pi, 256)), ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    ax.YDir = 'normal';
end

function toggleLabels(src, labels)
    if src.Value
        set(labels, 'Visible', 'on');
    else
        set(labels, 'Visible', 'off');
    end
end

function startDrag(fig, ax)
    data = getappdata(fig, 'colorSchemeEditData');
    if isempty(data)
        return;
    end
    cp = ax.CurrentPoint;
    x0 = cp(1,1);
    y0 = cp(1,2);
    dist = hypot(data.x - x0, data.y - y0);
    [minDist, idx] = min(dist);
    if minDist <= 0.08
        data.dragIdx = idx;
        data.selectedIdx = idx;
        updateSelectionUI(data);
        setappdata(fig, 'colorSchemeEditData', data);
    end
end

function dragMove(fig, ax)
    data = getappdata(fig, 'colorSchemeEditData');
    if isempty(data) || isempty(data.dragIdx)
        return;
    end

    cp = ax.CurrentPoint;
    x0 = cp(1,1);
    y0 = cp(1,2);
    r = hypot(x0, y0);
    if r > 1
        x0 = x0 / r;
        y0 = y0 / r;
    end

    idx = data.dragIdx;
    data.x(idx) = x0;
    data.y(idx) = y0;

    [h, s] = xyToHSV(x0, y0);
    data.hsv(idx, 1) = h;
    data.hsv(idx, 2) = s;
    data.rgb(idx, :) = hsv2rgb(data.hsv(idx, :));
    data.colorScheme.color(idx, :) = data.rgb(idx, :);

    data.scatter.XData = data.x;
    data.scatter.YData = data.y;
    data.scatter.CData = data.rgb;

    data.labels(idx).Position(1:2) = [x0 + data.labelOffset(1), y0 + data.labelOffset(2)];
    [labelColor, labelBg] = pickLabelColors(data.rgb(idx, :));
    data.labels(idx).Color = labelColor;
    data.labels(idx).BackgroundColor = labelBg;

    if data.selectedIdx == idx
        updateSelectionUI(data);
    end

    setappdata(fig, 'colorSchemeEditData', data);
end

function stopDrag(fig)
    data = getappdata(fig, 'colorSchemeEditData');
    if isempty(data)
        return;
    end
    data.dragIdx = [];
    setappdata(fig, 'colorSchemeEditData', data);
end

function onClose(fig)
    uiresume(fig);
end

function updateValueFromSlider(fig, slider)
    data = getappdata(fig, 'colorSchemeEditData');
    if isempty(data) || isempty(data.selectedIdx)
        return;
    end
    idx = data.selectedIdx;
    val = slider.Value;
    data.hsv(idx, 3) = val;
    data.rgb(idx, :) = hsv2rgb(data.hsv(idx, :));
    data.colorScheme.color(idx, :) = data.rgb(idx, :);

    data.scatter.CData = data.rgb;
    [labelColor, labelBg] = pickLabelColors(data.rgb(idx, :));
    data.labels(idx).Color = labelColor;
    data.labels(idx).BackgroundColor = labelBg;

    updateSelectionUI(data);
    setappdata(fig, 'colorSchemeEditData', data);
end

function updateSelectionUI(data)
    if isempty(data.selectedIdx)
        return;
    end
    idx = data.selectedIdx;
    if isgraphics(data.colorPanel)
        data.colorPanel.BackgroundColor = data.rgb(idx, :);
    end
    if isgraphics(data.valueSlider)
        data.valueSlider.Value = data.hsv(idx, 3);
    end
end

function [x, y] = hsvToXY(hsv)
    hue = hsv(:,1);
    sat = hsv(:,2);
    ang = 2 * pi * hue;
    x = sat .* cos(ang);
    y = sat .* sin(ang);
end

function [h, s] = xyToHSV(x, y)
    [ang, s] = cart2pol(x, y);
    h = mod(ang / (2 * pi), 1);
end

function [textColor, bgColor] = pickLabelColors(pointColor)
    lum = 0.2126 * pointColor(1) + 0.7152 * pointColor(2) + 0.0722 * pointColor(3);
    if lum > 0.5
        textColor = [0 0 0];
        bgColor = [1 1 1];
    else
        textColor = [1 1 1];
        bgColor = [0 0 0];
    end
end
