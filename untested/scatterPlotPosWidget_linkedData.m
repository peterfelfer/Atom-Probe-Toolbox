function out = scatterPlotPosWidget_linkedData(linkKey, speciesIndex, component, posVar)
% SCATTERPLOTPOSWIDGET_LINKEDDATA Helper for linkdata-driven updates.
%
% This function is evaluated via X/Y/ZDataSource strings to trigger
% refresh when the linked base workspace variable changes.

out = [];
if nargin < 3
    return;
end
if nargin < 4
    posVar = [];
end

if ~(ischar(linkKey) || isstring(linkKey))
    return;
end
key = char(linkKey);
if isempty(key) || ~isappdata(0, key)
    return;
end

linkState = getappdata(0, key);
if ~isstruct(linkState) || ~isfield(linkState, 'controlFig')
    return;
end

controlFig = linkState.controlFig;
if isempty(controlFig) || ~isgraphics(controlFig)
    return;
end

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data)
    return;
end

if ~istable(posVar)
    if isfield(linkState, 'linkedVar') && ischar(linkState.linkedVar)
        try
            posVar = evalin('base', linkState.linkedVar);
        catch
            posVar = [];
        end
    end
end

if istable(posVar)
    if ~isfield(linkState, 'updating') || ~linkState.updating
        if isfield(linkState, 'sigFcn') && isa(linkState.sigFcn, 'function_handle')
            newSig = linkState.sigFcn(posVar);
        else
            newSig = [];
        end
        isSame = ~isempty(newSig) && isfield(linkState, 'lastSig') && isequal(newSig, linkState.lastSig);
        if ~isSame
            if isfield(linkState, 'refreshFcn') && isa(linkState.refreshFcn, 'function_handle')
                linkState.updating = true;
                setappdata(0, key, linkState);
                data = getappdata(controlFig, 'scatterPlotPosWidget');
                if ~isempty(data)
                    data.linkUpdating = true;
                    setappdata(controlFig, 'scatterPlotPosWidget', data);
                end
                try
                    linkState.refreshFcn(posVar);
                catch
                end
                data = getappdata(controlFig, 'scatterPlotPosWidget');
                if ~isempty(data)
                    data.linkUpdating = false;
                    setappdata(controlFig, 'scatterPlotPosWidget', data);
                end
                linkState = getappdata(0, key);
                linkState.lastSig = newSig;
                linkState.updating = false;
                setappdata(0, key, linkState);
            else
                linkState.lastSig = newSig;
                setappdata(0, key, linkState);
            end
        end
    end
end

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data) || ~isfield(data, 'scatterHandles')
    return;
end
if speciesIndex < 1 || speciesIndex > numel(data.scatterHandles)
    return;
end

h = data.scatterHandles(speciesIndex);
if ~isgraphics(h)
    return;
end

switch lower(component)
    case 'x'
        out = h.XData;
    case 'y'
        out = h.YData;
    case 'z'
        out = h.ZData;
    case 'c'
        out = h.CData;
    otherwise
        out = h.XData;
end
end
