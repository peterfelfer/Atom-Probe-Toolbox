function [peakTable, peakInfo, params] = tunePeakDetection(posOrMc, options)
% TUNEPEAKDETECTION Interactive tuning of peak detection parameters.
%
% [peakTable, peakInfo, params] = tunePeakDetection(pos)
% [peakTable, peakInfo, params] = tunePeakDetection(mc)
%
% Opens a figure with the mass spectrum, detected peaks, baseline, and
% noise floor. Sliders let you adjust the detection parameters in real
% time. Press "Accept" when satisfied.
%
% INPUT
%   posOrMc - pos table with mc column, or numeric mc vector (Da)
%
% OPTIONS
%   'binWidth'            - initial histogram bin width (default: 0.01)
%   'minProminenceFactor' - initial prominence factor (default: 6)
%   'minDistanceDa'       - initial min peak distance (default: 0.3)
%   'baselineQuantile'    - initial baseline quantile (default: 0.1)
%   'smoothSpan'          - initial smoothing span (default: 0.1)
%   'minProminence'       - initial min prominence for display (default: 500)
%   'mcRange'             - [min max] x-axis range (default: [0 200])
%
% OUTPUT
%   peakTable - peak table from massSpecFindPeaks with final parameters
%   peakInfo  - info struct from massSpecFindPeaks
%   params    - struct with the accepted parameter values
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    posOrMc
    options.binWidth (1,1) double = 0.01
    options.minProminenceFactor (1,1) double = 6
    options.minDistanceDa (1,1) double = 0.3
    options.baselineQuantile (1,1) double = 0.1
    options.smoothSpan (1,1) double = 0.1
    options.minProminence (1,1) double = 500
    options.mcRange (1,2) double = [0 200]
end

% Resolve mc
if istable(posOrMc)
    mc = posOrMc.mc;
else
    mc = posOrMc;
end
mc = mc(:);
mc = mc(isfinite(mc) & mc > 0);

% Store state
state = struct();
state.mc = mc;
state.binWidth = options.binWidth;
state.promFactor = options.minProminenceFactor;
state.minDist = options.minDistanceDa;
state.blQuantile = options.baselineQuantile;
state.smoothSpan = options.smoothSpan;
state.minProm = options.minProminence;
state.mcRange = options.mcRange;

% Initial detection
[peakTable, peakInfo] = runDetection(state);

% --- Build figure ---
fig = figure('Name', 'Tune Peak Detection', 'NumberTitle', 'off', ...
             'Position', [80 80 1300 700]);

ax = axes(fig, 'Position', [0.06 0.32 0.90 0.63]);
plotSpectrum(ax, peakInfo, peakTable, state);

% --- Slider panel ---
sliderY = 0.22;
sliderH = 0.03;
labelW  = 0.12;
sliderW = 0.28;
valW    = 0.06;
gap     = 0.005;

% Column 1: prominence factor, min distance, min prominence display
col1x = 0.02;
[state, fig] = addSlider(fig, state, 'promFactor', 'Prominence factor', ...
    col1x, sliderY, labelW, sliderW, valW, sliderH, 1, 20, @onParamChange);
[state, fig] = addSlider(fig, state, 'minDist', 'Min distance (Da)', ...
    col1x, sliderY - sliderH - gap*4, labelW, sliderW, valW, sliderH, 0.05, 2, @onParamChange);
[state, fig] = addSlider(fig, state, 'minProm', 'Min prominence', ...
    col1x, sliderY - 2*(sliderH + gap*4), labelW, sliderW, valW, sliderH, 10, 10000, @onParamChange);

% Column 2: baseline quantile, smooth span
col2x = 0.52;
[state, fig] = addSlider(fig, state, 'blQuantile', 'Baseline quantile', ...
    col2x, sliderY, labelW, sliderW, valW, sliderH, 0.01, 0.5, @onParamChange);
[state, fig] = addSlider(fig, state, 'smoothSpan', 'Smooth span (Da)', ...
    col2x, sliderY - sliderH - gap*4, labelW, sliderW, valW, sliderH, 0.01, 1.0, @onParamChange);

% Accept button
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Accept', ...
    'Units', 'normalized', 'Position', [col2x, sliderY - 2*(sliderH+gap*4), 0.08, 0.04], ...
    'FontSize', 11, 'FontWeight', 'bold', ...
    'Callback', @(~,~) uiresume(fig));

% Peak count label
state.peakCountLabel = uicontrol(fig, 'Style', 'text', ...
    'String', sprintf('Peaks: %d', height(peakTable(peakTable.prominence > state.minProm, :))), ...
    'Units', 'normalized', 'Position', [col2x + 0.10, sliderY - 2*(sliderH+gap*4), 0.15, 0.04], ...
    'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');

setappdata(fig, 'state', state);
setappdata(fig, 'ax', ax);

uiwait(fig);

% Return final results
if isvalid(fig)
    state = getappdata(fig, 'state');
    close(fig);
else
    % Figure was closed — use last state
end

[peakTable, peakInfo] = runDetection(state);
params = struct( ...
    'binWidth', state.binWidth, ...
    'minProminenceFactor', state.promFactor, ...
    'minDistanceDa', state.minDist, ...
    'baselineQuantile', state.blQuantile, ...
    'smoothSpan', state.smoothSpan, ...
    'minProminence', state.minProm);

fprintf('Accepted parameters: promFactor=%.1f, minDist=%.2f, blQuantile=%.2f, smoothSpan=%.2f, minProm=%.0f\n', ...
    params.minProminenceFactor, params.minDistanceDa, params.baselineQuantile, ...
    params.smoothSpan, params.minProminence);
fprintf('Peaks above threshold: %d\n', height(peakTable(peakTable.prominence > params.minProminence, :)));

end


% ===== Helper functions ==================================================

function [peakTable, peakInfo] = runDetection(state)
    [peakTable, peakInfo] = massSpecFindPeaks(state.mc, ...
        'binWidth', state.binWidth, ...
        'minProminenceFactor', state.promFactor, ...
        'minDistanceDa', state.minDist, ...
        'baselineQuantile', state.blQuantile, ...
        'smoothSpan', state.smoothSpan);
end


function plotSpectrum(ax, peakInfo, peakTable, state)
    cla(ax);
    sigPeaks = peakTable(peakTable.prominence > state.minProm, :);

    semilogy(ax, peakInfo.centers, peakInfo.counts, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    hold(ax, 'on');
    semilogy(ax, peakInfo.centers, peakInfo.baseline, 'b-', 'LineWidth', 1.2);
    semilogy(ax, peakInfo.centers, peakInfo.baseline + peakInfo.noise, 'g--', 'LineWidth', 0.8);

    if ~isempty(sigPeaks)
        semilogy(ax, sigPeaks.mc, sigPeaks.height, 'rv', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
        for i = 1:height(sigPeaks)
            if sigPeaks.height(i) > 3000
                text(ax, sigPeaks.mc(i), sigPeaks.height(i)*1.8, sprintf('%.1f', sigPeaks.mc(i)), ...
                    'FontSize', 7, 'HorizontalAlignment', 'center', 'Color', 'r');
            end
        end
    end

    xlim(ax, state.mcRange);
    ylim(ax, [1 max(peakInfo.counts)*3]);
    xlabel(ax, 'Mass-to-charge (Da)');
    ylabel(ax, 'Counts per bin');
    title(ax, sprintf('Peak detection: %d peaks above threshold', height(sigPeaks)));
    hold(ax, 'off');
end


function onParamChange(fig)
    if ~isvalid(fig), return; end
    state = getappdata(fig, 'state');
    ax    = getappdata(fig, 'ax');

    [peakTable, peakInfo] = runDetection(state);
    plotSpectrum(ax, peakInfo, peakTable, state);

    nSig = height(peakTable(peakTable.prominence > state.minProm, :));
    state.peakCountLabel.String = sprintf('Peaks: %d', nSig);

    setappdata(fig, 'state', state);
    drawnow limitrate;
end


function [state, fig] = addSlider(fig, state, fieldName, label, x, y, labelW, sliderW, valW, h, minVal, maxVal, callback)
    curVal = state.(fieldName);

    uicontrol(fig, 'Style', 'text', 'String', [label ':'], ...
        'Units', 'normalized', 'Position', [x, y, labelW, h], ...
        'HorizontalAlignment', 'right', 'FontSize', 9);

    sl = uicontrol(fig, 'Style', 'slider', ...
        'Min', minVal, 'Max', maxVal, 'Value', curVal, ...
        'Units', 'normalized', 'Position', [x + labelW + 0.005, y, sliderW, h], ...
        'SliderStep', [0.01, 0.1]);

    valLabel = uicontrol(fig, 'Style', 'text', ...
        'String', formatVal(curVal), ...
        'Units', 'normalized', 'Position', [x + labelW + sliderW + 0.01, y, valW, h], ...
        'HorizontalAlignment', 'left', 'FontSize', 9);

    addlistener(sl, 'Value', 'PostSet', @(~,~) sliderMoved(fig, sl, valLabel, fieldName, callback));
end


function sliderMoved(fig, sl, valLabel, fieldName, callback)
    if ~isvalid(fig), return; end
    val = sl.Value;
    state = getappdata(fig, 'state');
    state.(fieldName) = val;
    setappdata(fig, 'state', state);
    valLabel.String = formatVal(val);
    callback(fig);
end


function s = formatVal(v)
    if v >= 100
        s = sprintf('%.0f', v);
    elseif v >= 1
        s = sprintf('%.2f', v);
    else
        s = sprintf('%.3f', v);
    end
end
