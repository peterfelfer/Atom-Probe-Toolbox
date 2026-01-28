function rangeExtendToNeighbor(spec)
% rangeExtendToNeighbor extends a range to touch its neighboring range
%
% rangeExtendToNeighbor(spec)
%
% Click on a range to extend it. The side of the range that is closer to
% where you clicked will be extended to touch the neighboring range on
% that side.
%
% INPUT
% spec:     area plot that displays the mass spectrum, or axes containing
%           the mass spectrum
%
% OUTPUT
% (none - modifies the range plot in place)
%
% USAGE
% Click on a range. If you click closer to the left edge, the range will
% extend left to meet the previous range. If you click closer to the right
% edge, the range will extend right to meet the next range.
%
% Right-click to cancel.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% Resolve axes
if isgraphics(spec, 'axes')
    ax = spec;
    % Find the mass spectrum plot
    plots = ax.Children;
    spec = [];
    for i = 1:numel(plots)
        try
            if plots(i).UserData.plotType == "massSpectrum"
                spec = plots(i);
                break;
            end
        catch
        end
    end
    if isempty(spec)
        % Fall back to first area plot
        for i = 1:numel(plots)
            if isa(plots(i), 'matlab.graphics.primitive.Area')
                spec = plots(i);
                break;
            end
        end
    end
    if isempty(spec)
        error('rangeExtendToNeighbor:noSpec', 'No mass spectrum found in axes.');
    end
else
    ax = spec.Parent;
end

%% Get all range plots sorted by position
plots = ax.Children;
rangePlots = [];
rangeBegins = [];
rangeEnds = [];

for i = 1:numel(plots)
    try
        if plots(i).UserData.plotType == "range"
            rangePlots = [rangePlots; plots(i)];
            rangeBegins = [rangeBegins; plots(i).XData(1)];
            rangeEnds = [rangeEnds; plots(i).XData(end)];
        end
    catch
    end
end

if isempty(rangePlots)
    error('rangeExtendToNeighbor:noRanges', 'No ranges found in the mass spectrum.');
end

% Sort ranges by position
[rangeBegins, sortIdx] = sort(rangeBegins);
rangeEnds = rangeEnds(sortIdx);
rangePlots = rangePlots(sortIdx);
nRanges = numel(rangePlots);

%% Get user click
[xClick, ~, button] = ginput(1);

% Right-click to cancel
if isempty(button) || button == 3
    return;
end

%% Find which range was clicked
clickedIdx = find(xClick >= rangeBegins & xClick <= rangeEnds, 1);

if isempty(clickedIdx)
    warning('rangeExtendToNeighbor:noRangeClicked', 'Click was not on a range.');
    return;
end

%% Determine which side of the range was clicked (closer to left or right boundary)
clickedBegin = rangeBegins(clickedIdx);
clickedEnd = rangeEnds(clickedIdx);
distToLeft = abs(xClick - clickedBegin);
distToRight = abs(xClick - clickedEnd);

if distToLeft < distToRight
    % Extend left side
    extendDirection = 'left';
else
    % Extend right side
    extendDirection = 'right';
end

%% Find the neighboring range and extend
clickedPlot = rangePlots(clickedIdx);

if strcmp(extendDirection, 'left')
    % Extend left to touch previous range
    if clickedIdx == 1
        warning('rangeExtendToNeighbor:noLeftNeighbor', 'No range to the left to extend to.');
        return;
    end

    % Get the end of the previous range
    prevRangeEnd = rangeEnds(clickedIdx - 1);
    newBegin = prevRangeEnd;

    % Update the range plot
    updateRangeLeft(clickedPlot, spec, newBegin);

else
    % Extend right to touch next range
    if clickedIdx == nRanges
        warning('rangeExtendToNeighbor:noRightNeighbor', 'No range to the right to extend to.');
        return;
    end

    % Get the begin of the next range
    nextRangeBegin = rangeBegins(clickedIdx + 1);
    newEnd = nextRangeBegin;

    % Update the range plot
    updateRangeRight(clickedPlot, spec, newEnd);

end

end

%% Helper functions

function updateRangeLeft(rangePlot, spec, newBegin)
    % Extend the left side of a range to newBegin
    oldBegin = rangePlot.XData(1);
    oldEnd = rangePlot.XData(end);

    if newBegin >= oldBegin
        warning('rangeExtendToNeighbor:alreadyTouching', 'Ranges are already touching or overlapping.');
        return;
    end

    % Get spectrum data for the new extended range
    isIn = (spec.XData >= newBegin) & (spec.XData <= oldEnd);

    if ~any(isIn)
        warning('rangeExtendToNeighbor:noData', 'No spectrum data in the extended range.');
        return;
    end

    % Update the range plot
    rangePlot.XData = spec.XData(isIn);
    rangePlot.YData = spec.YData(isIn);

    % Update associated text position if exists
    updateRangeText(rangePlot);
end

function updateRangeRight(rangePlot, spec, newEnd)
    % Extend the right side of a range to newEnd
    oldBegin = rangePlot.XData(1);
    oldEnd = rangePlot.XData(end);

    if newEnd <= oldEnd
        warning('rangeExtendToNeighbor:alreadyTouching', 'Ranges are already touching or overlapping.');
        return;
    end

    % Get spectrum data for the new extended range
    isIn = (spec.XData >= oldBegin) & (spec.XData <= newEnd);

    if ~any(isIn)
        warning('rangeExtendToNeighbor:noData', 'No spectrum data in the extended range.');
        return;
    end

    % Update the range plot
    rangePlot.XData = spec.XData(isIn);
    rangePlot.YData = spec.YData(isIn);

    % Update associated text position if exists
    updateRangeText(rangePlot);
end

function updateRangeText(rangePlot)
    % Find and update the text label associated with this range
    ax = rangePlot.Parent;
    plots = ax.Children;

    for i = 1:numel(plots)
        try
            if plots(i).UserData.plotType == "text"
                % Check if this text belongs to the range by matching DisplayName
                if isprop(plots(i), 'DisplayName') && isprop(rangePlot, 'DisplayName')
                    if strcmp(plots(i).DisplayName, rangePlot.DisplayName)
                        % Update text position to new range start
                        plots(i).Position(1) = rangePlot.XData(1);
                        plots(i).Position(2) = max(rangePlot.YData) * 1.4;
                        return;
                    end
                end
            end
        catch
        end
    end
end
