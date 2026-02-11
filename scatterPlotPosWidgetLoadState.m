function [controlFig, loaded, sourceName] = scatterPlotPosWidgetLoadState(controlFig, source)
% SCATTERPLOTPOSWIDGETLOADSTATE Load and apply state/profile to a running widget.
%
% [controlFig, loaded, sourceName] = scatterPlotPosWidgetLoadState(controlFig)
% [controlFig, loaded, sourceName] = scatterPlotPosWidgetLoadState(controlFig, source)
%
% INPUT:
%   controlFig - active scatterPlotPosWidget control figure (default: gcf)
%   source     - workspace variable name (char/string) or struct payload
%
% OUTPUT:
%   loaded     - true if a state/profile was applied
%   sourceName - workspace variable name, or "<struct>" when source is struct

if nargin < 1 || isempty(controlFig)
    controlFig = gcf;
end
if isgraphics(controlFig, 'axes') || isgraphics(controlFig, 'patch')
    controlFig = ancestor(controlFig, 'figure');
end

if ~isgraphics(controlFig, 'figure')
    error('scatterPlotPosWidgetLoadState:invalidHandle', ...
        'controlFig must be a valid figure handle.');
end

loaded = false;
sourceName = "";

if nargin < 2 || isempty(source)
    [sourceName, isSelected] = selectWorkspaceStructVariable();
    if ~isSelected
        return;
    end
    source = sourceName;
end

if ischar(source) || (isstring(source) && isscalar(source))
    sourceName = string(source);
    sourceName = strtrim(sourceName);
    if strlength(sourceName) == 0
        error('scatterPlotPosWidgetLoadState:emptyVariableName', ...
            'Workspace variable name must not be empty.');
    end
    if ~isvarname(char(sourceName))
        error('scatterPlotPosWidgetLoadState:invalidVariableName', ...
            'Invalid workspace variable name "%s".', char(sourceName));
    end
    if ~evalin('base', sprintf('exist(''%s'', ''var'')', char(sourceName)))
        error('scatterPlotPosWidgetLoadState:missingVariable', ...
            'Workspace variable "%s" does not exist.', char(sourceName));
    end
    payload = evalin('base', char(sourceName));
elseif isstruct(source)
    payload = source;
    sourceName = "<struct>";
else
    error('scatterPlotPosWidgetLoadState:invalidSource', ...
        'source must be a workspace variable name or a scalar struct.');
end

if ~isstruct(payload) || ~isscalar(payload)
    error('scatterPlotPosWidgetLoadState:invalidPayload', ...
        'Source payload must be a scalar struct.');
end

state = visualisationProfileToWidgetState(payload);
scatterPlotPosWidgetApplyState(controlFig, state);
loaded = true;

end

function [sourceName, isSelected] = selectWorkspaceStructVariable()
sourceName = "";
isSelected = false;

wsInfo = evalin('base', 'whos');
if isempty(wsInfo)
    errordlg('Workspace is empty. No profile/state variable available.', ...
        'Load From Workspace', 'modal');
    return;
end

isStruct = strcmp({wsInfo.class}, 'struct');
names = string({wsInfo(isStruct).name});
if isempty(names)
    errordlg('No scalar struct variables found in workspace.', ...
        'Load From Workspace', 'modal');
    return;
end

% Keep only scalar struct variables (likely profile/state candidates).
isScalar = false(size(names));
for i = 1:numel(names)
    varName = names(i);
    try
        value = evalin('base', char(varName));
        isScalar(i) = isstruct(value) && isscalar(value);
    catch
        isScalar(i) = false;
    end
end
names = names(isScalar);
if isempty(names)
    errordlg('No scalar struct variables found in workspace.', ...
        'Load From Workspace', 'modal');
    return;
end

[~, order] = sort(lower(names));
names = names(order);

defaultIdx = find(names == "visualisationProfile", 1, 'first');
if isempty(defaultIdx)
    defaultIdx = 1;
end

[idx, tf] = listdlg( ...
    'ListString', cellstr(names), ...
    'SelectionMode', 'single', ...
    'InitialValue', defaultIdx, ...
    'PromptString', 'Select visualisation profile/state variable:', ...
    'Name', 'Load Profile From Workspace', ...
    'ListSize', [360 280]);

if ~tf || isempty(idx)
    return;
end

sourceName = names(idx);
isSelected = true;
end
