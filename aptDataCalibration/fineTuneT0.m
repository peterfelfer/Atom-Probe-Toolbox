function [pos, t0Final] = fineTuneT0(pos, flightPathLength, t0Initial, options)
% FINETUNET0 Interactively fine-tune the propagation delay t0.
%
% [pos, t0] = fineTuneT0(pos, flightPathLength, t0Initial)
% [pos, t0] = fineTuneT0(pos, 110, 38, 'mode', 'voltage', 'VP', pos.VP)
%
% Opens a figure with a mass spectrum histogram and a slider for t0.
% Adjust the slider until a known reference peak aligns with its expected
% mass. Press "Accept" or close the figure to confirm.
%
% INPUT
%   pos              - table with tof (ns), VDC (V), detx (mm), dety (mm)
%   flightPathLength - nominal flight path (mm)
%   t0Initial        - starting value for t0 (ns)
%
% OPTIONS
%   'mode'     - 'laser' (default) or 'voltage'
%   'VP'       - pulse voltage vector (V), required for voltage mode
%   'maxMc'    - x-axis limit in Da (default 200)
%   'binWidth' - histogram bin width in Da (default 0.1)
%   't0Range'  - slider range: t0Initial +/- t0Range (default 200)
%
% OUTPUT
%   pos     - table with mc column added/updated
%   t0Final - accepted t0 value (ns)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    flightPathLength (1,1) double
    t0Initial (1,1) double
    options.mode (1,:) char {mustBeMember(options.mode, {'laser','voltage'})} = 'laser'
    options.VP (:,1) double = []
    options.maxMc (1,1) double = 200
    options.binWidth (1,1) double = 0.1
    options.t0Range (1,1) double = 200
end

t0Min = max(0, t0Initial - options.t0Range);
t0Max = t0Initial + options.t0Range;

% Store data for callbacks
tof  = pos.tof;
VDC  = pos.VDC;
detx = pos.detx;
dety = pos.dety;

% Compute initial mc
mcArgs = {'mode', options.mode};
if strcmp(options.mode, 'voltage')
    mcArgs = [mcArgs, {'VP', options.VP}];
end
mc = tofToMassToCharge(tof, VDC, detx, dety, flightPathLength, t0Initial, mcArgs{:});
mc(mc < 0) = 0;

% Figure
edges = 0:options.binWidth:options.maxMc;
fig = figure('Name', 'Fine-tune t0', 'NumberTitle', 'off', ...
             'Position', [100 100 900 550]);
ax = axes(fig, 'Position', [0.08 0.22 0.88 0.72]);
counts = histcounts(mc, edges);
bar(ax, edges(1:end-1), counts, 'histc', 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
set(ax, 'YScale', 'log');
xlim(ax, [0 options.maxMc]);
xlabel(ax, 'Mass-to-charge (Da)');
ylabel(ax, 'Counts');
title(ax, sprintf('t_0 = %.2f ns', t0Initial));

% Slider
uicontrol(fig, 'Style', 'text', 'String', 't0 (ns):', ...
    'Units', 'normalized', 'Position', [0.02 0.04 0.07 0.04], ...
    'HorizontalAlignment', 'right');

slider = uicontrol(fig, 'Style', 'slider', ...
    'Min', t0Min, 'Max', t0Max, 'Value', t0Initial, ...
    'Units', 'normalized', 'Position', [0.10 0.04 0.70 0.04], ...
    'SliderStep', [1/(t0Max-t0Min), 10/(t0Max-t0Min)]);

valLabel = uicontrol(fig, 'Style', 'text', ...
    'String', sprintf('%.2f', t0Initial), ...
    'Units', 'normalized', 'Position', [0.81 0.04 0.08 0.04]);

uicontrol(fig, 'Style', 'pushbutton', 'String', 'Accept', ...
    'Units', 'normalized', 'Position', [0.90 0.03 0.08 0.06], ...
    'Callback', @(~,~) uiresume(fig));

% Store state for callback
setappdata(fig, 'ax', ax);
setappdata(fig, 'tof', tof);
setappdata(fig, 'VDC', VDC);
setappdata(fig, 'detx', detx);
setappdata(fig, 'dety', dety);
setappdata(fig, 'edges', edges);
setappdata(fig, 'valLabel', valLabel);
setappdata(fig, 'mcArgs', {mcArgs});
setappdata(fig, 'flightPathLength', flightPathLength);

addlistener(slider, 'Value', 'PostSet', @(~,~) updatePlot(fig, slider));

uiwait(fig);

if isvalid(fig)
    t0Final = slider.Value;
    close(fig);
else
    t0Final = t0Initial;
    fprintf('Figure closed — keeping t0 = %.2f ns\n', t0Initial);
end

% Apply final mc
pos.mc = tofToMassToCharge(tof, VDC, detx, dety, flightPathLength, t0Final, mcArgs{:});
pos.mc(pos.mc < 0) = 0;

if ~ismember('mc', pos.Properties.VariableNames)
    pos.Properties.VariableUnits{end} = 'Da';
end

fprintf('Accepted t0 = %.2f ns\n', t0Final);

end


function updatePlot(fig, slider)
    if ~isvalid(fig), return; end

    t0New = slider.Value;
    ax    = getappdata(fig, 'ax');
    edges = getappdata(fig, 'edges');
    tof   = getappdata(fig, 'tof');
    VDC   = getappdata(fig, 'VDC');
    detx  = getappdata(fig, 'detx');
    dety  = getappdata(fig, 'dety');
    L     = getappdata(fig, 'flightPathLength');
    args  = getappdata(fig, 'mcArgs');
    label = getappdata(fig, 'valLabel');

    mc = tofToMassToCharge(tof, VDC, detx, dety, L, t0New, args{1}{:});
    mc(mc < 0) = 0;
    counts = histcounts(mc, edges);

    cla(ax);
    bar(ax, edges(1:end-1), counts, 'histc', 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
    set(ax, 'YScale', 'log');
    xlim(ax, [0 max(edges)]);
    xlabel(ax, 'Mass-to-charge (Da)');
    ylabel(ax, 'Counts');
    title(ax, sprintf('t_0 = %.2f ns', t0New));
    label.String = sprintf('%.2f', t0New);
    drawnow limitrate;
end
