function [posOut, info, hullRec] = reconstructionWidget(pos, options)
% reconstructionWidget interactive exploration of APT reconstructions.
%
% [posOut, info, hullRec] = reconstructionWidget(pos)
% [posOut, info, hullRec] = reconstructionWidget(pos, options)
%
% options fields (optional):
%   title           figure title
%   sample          fraction (0..1) or count (>1) for scatter preview
%   liveUpdate      true/false (default: true)
%   volumeTable     ion volume table (ion, ionVolume)
%   isotopeTable    isotope table for ionVolumesFromAtomVolumes
%   hullOptions     options passed to scanHull
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.title (1,1) string = "Reconstruction Explorer"
    options.sample (1,1) double = 0.05
    options.liveUpdate (1,1) logical = true
    options.volumeTable table = table()
    options.isotopeTable = []
    options.hullOptions = struct()
end

validatePos(pos);

[posOut, info, hullRec] = deal(pos, struct(), []);

fig = figure('Name', char(options.title), 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'pixels', 'Position', [100 100 980 620]);

ax = axes(fig, 'Units', 'normalized', 'Position', [0.05 0.10 0.60 0.85]);
axis(ax, 'vis3d');
grid(ax, 'on');

ctrl = uipanel(fig, 'Units', 'normalized', 'Position', [0.68 0.06 0.30 0.90], ...
    'Title', 'Reconstruction Controls');

data = initData(pos, options, ax, ctrl);
setappdata(fig, 'reconstructionWidgetData', data);

buildControls(fig);
applyReconstruction(fig, true);
uiwait(fig);

data = getappdata(fig, 'reconstructionWidgetData');
if isempty(data)
    return;
end
posOut = data.posOut;
info = data.info;
hullRec = data.hullRec;
delete(fig);
end

function data = initData(pos, options, ax, ctrl)
    data = struct();
    data.pos = pos;
    data.posOut = pos;
    data.info = struct();
    data.hullRec = [];
    data.ax = ax;
    data.ctrl = ctrl;
    data.scatter = [];
    data.hullPatch = [];
    data.sample = options.sample;
    data.liveUpdate = options.liveUpdate;
    data.volumeTable = options.volumeTable;
    data.isotopeTable = options.isotopeTable;
    data.hullOptions = options.hullOptions;

    data.model = "voltage";
    data.voltageField = detectVoltageField(pos);
    data.kf = 10;
    data.ICF = 1.65;
    data.Fevap = 65;
    data.flightLength = 110;
    data.detEff = 0.82;
    data.shankAngle = 6;
    data.radius0 = NaN;
    data.customFcn = [];
    data.reconstructHull = false;
end

function buildControls(fig)
    data = getappdata(fig, 'reconstructionWidgetData');
    ctrl = data.ctrl;

    addText(ctrl, [0.05 0.92 0.9 0.04], 'Radius evolution');
    data.modelPopup = uicontrol(ctrl, 'Style', 'popupmenu', ...
        'String', {'voltage','shank','custom'}, ...
        'Value', 1, 'Units', 'normalized', ...
        'Position', [0.05 0.88 0.9 0.04], ...
        'Callback', @(src,~) onModelChange(fig, src));

    data.liveCheck = uicontrol(ctrl, 'Style', 'checkbox', ...
        'String', 'Live update', ...
        'Value', data.liveUpdate, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.83 0.9 0.04], ...
        'Callback', @(src,~) onLiveToggle(fig, src));

    data.updateBtn = uicontrol(ctrl, 'Style', 'pushbutton', ...
        'String', 'Update', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.78 0.9 0.05], ...
        'Callback', @(~,~) applyReconstruction(fig, false));

    addText(ctrl, [0.05 0.72 0.45 0.04], 'V field');
    data.voltagePopup = uicontrol(ctrl, 'Style', 'popupmenu', ...
        'String', voltageFieldList(data.pos), ...
        'Value', voltageFieldIndex(data.pos, data.voltageField), ...
        'Units', 'normalized', ...
        'Position', [0.50 0.72 0.45 0.04], ...
        'Callback', @(src,~) onVoltageField(fig, src));

    addText(ctrl, [0.05 0.67 0.45 0.04], 'kf');
    data.kfEdit = addEdit(ctrl, [0.50 0.67 0.45 0.04], data.kf, @(src,~) onParamEdit(fig, src, 'kf'));
    addText(ctrl, [0.05 0.62 0.45 0.04], 'ICF');
    data.icfEdit = addEdit(ctrl, [0.50 0.62 0.45 0.04], data.ICF, @(src,~) onParamEdit(fig, src, 'ICF'));
    addText(ctrl, [0.05 0.57 0.45 0.04], 'Fevap [V/nm]');
    data.fevapEdit = addEdit(ctrl, [0.50 0.57 0.45 0.04], data.Fevap, @(src,~) onParamEdit(fig, src, 'Fevap'));
    addText(ctrl, [0.05 0.52 0.45 0.04], 'flight [mm]');
    data.flightEdit = addEdit(ctrl, [0.50 0.52 0.45 0.04], data.flightLength, @(src,~) onParamEdit(fig, src, 'flightLength'));
    addText(ctrl, [0.05 0.47 0.45 0.04], 'detEff');
    data.detEffEdit = addEdit(ctrl, [0.50 0.47 0.45 0.04], data.detEff, @(src,~) onParamEdit(fig, src, 'detEff'));

    addText(ctrl, [0.05 0.41 0.45 0.04], 'shank [deg]');
    data.shankEdit = addEdit(ctrl, [0.50 0.41 0.45 0.04], data.shankAngle, @(src,~) onParamEdit(fig, src, 'shankAngle'));
    addText(ctrl, [0.05 0.36 0.45 0.04], 'R0 [nm]');
    data.r0Edit = addEdit(ctrl, [0.50 0.36 0.45 0.04], data.radius0, @(src,~) onParamEdit(fig, src, 'radius0'));

    addText(ctrl, [0.05 0.30 0.9 0.04], 'Custom radius function');
    data.customEdit = uicontrol(ctrl, 'Style', 'edit', ...
        'String', '', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.26 0.9 0.04], ...
        'Callback', @(src,~) onCustomFcn(fig, src));
    data.fitBtn = uicontrol(ctrl, 'Style', 'pushbutton', ...
        'String', 'Fit from image', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.22 0.9 0.04], ...
        'Callback', @(~,~) onFitFromImage(fig));

    data.hullCheck = uicontrol(ctrl, 'Style', 'checkbox', ...
        'String', 'Reconstruct hull', ...
        'Value', data.reconstructHull, ...
        'Units', 'normalized', ...
        'Position', [0.05 0.16 0.9 0.04], ...
        'Callback', @(src,~) onHullToggle(fig, src));

    addText(ctrl, [0.05 0.11 0.45 0.04], 'sample');
    data.sampleEdit = addEdit(ctrl, [0.50 0.11 0.45 0.04], data.sample, @(src,~) onParamEdit(fig, src, 'sample'));

    setappdata(fig, 'reconstructionWidgetData', data);
end

function addText(parent, pos, str)
    uicontrol(parent, 'Style', 'text', 'String', str, ...
        'Units', 'normalized', 'HorizontalAlignment', 'left', ...
        'Position', pos, 'BackgroundColor', get(parent, 'BackgroundColor'));
end

function h = addEdit(parent, pos, value, cb)
    h = uicontrol(parent, 'Style', 'edit', 'String', num2str(value), ...
        'Units', 'normalized', 'Position', pos, 'Callback', cb);
end

function onModelChange(fig, src)
    data = getappdata(fig, 'reconstructionWidgetData');
    val = src.Value;
    items = src.String;
    data.model = string(items{val});
    setappdata(fig, 'reconstructionWidgetData', data);
    maybeUpdate(fig);
end

function onLiveToggle(fig, src)
    data = getappdata(fig, 'reconstructionWidgetData');
    data.liveUpdate = logical(src.Value);
    setappdata(fig, 'reconstructionWidgetData', data);
end

function onVoltageField(fig, src)
    data = getappdata(fig, 'reconstructionWidgetData');
    items = src.String;
    data.voltageField = string(items{src.Value});
    setappdata(fig, 'reconstructionWidgetData', data);
    maybeUpdate(fig);
end

function onParamEdit(fig, src, field)
    data = getappdata(fig, 'reconstructionWidgetData');
    val = str2double(src.String);
    if ~isfinite(val)
        return;
    end
    data.(field) = val;
    setappdata(fig, 'reconstructionWidgetData', data);
    maybeUpdate(fig);
end

function onCustomFcn(fig, src)
    data = getappdata(fig, 'reconstructionWidgetData');
    str = strtrim(src.String);
    if isempty(str)
        data.customFcn = [];
    else
        try
            data.customFcn = str2func(str);
        catch
            warning('reconstructionWidget:invalidCustomFcn', 'Invalid function handle.');
        end
    end
    setappdata(fig, 'reconstructionWidgetData', data);
    maybeUpdate(fig);
end

function onFitFromImage(fig)
    data = getappdata(fig, 'reconstructionWidgetData');
    [fcn, fcnText] = fitRadiusFromImage();
    if ~isempty(fcn)
        data.customFcn = fcn;
        data.customEdit.String = fcnText;
        data.model = "custom";
        data.modelPopup.Value = 3;
        setappdata(fig, 'reconstructionWidgetData', data);
        maybeUpdate(fig);
    end
end

function onHullToggle(fig, src)
    data = getappdata(fig, 'reconstructionWidgetData');
    data.reconstructHull = logical(src.Value);
    setappdata(fig, 'reconstructionWidgetData', data);
    maybeUpdate(fig);
end

function maybeUpdate(fig)
    data = getappdata(fig, 'reconstructionWidgetData');
    if data.liveUpdate
        applyReconstruction(fig, false);
    end
end

function applyReconstruction(fig, initial)
    data = getappdata(fig, 'reconstructionWidgetData');
    opts = struct();
    opts.evolutionModel = data.model;
    opts.voltageField = data.voltageField;
    opts.kf = data.kf;
    opts.ICF = data.ICF;
    opts.Fevap = data.Fevap;
    opts.flightLength = data.flightLength;
    opts.detEff = data.detEff;
    opts.shankAngle = data.shankAngle;
    opts.radius0 = data.radius0;
    opts.customRadiusFcn = data.customFcn;
    opts.ionVolumeTable = data.volumeTable;
    opts.isotopeTable = data.isotopeTable;
    opts.reconstructHull = data.reconstructHull;
    opts.hullOptions = data.hullOptions;

    [posOut, info, hullRec] = posReconstructFromDetectorAdvanced(data.pos, opts);
    data.posOut = posOut;
    data.info = info;
    data.hullRec = hullRec;

    data = updatePlot(data, initial);
    setappdata(fig, 'reconstructionWidgetData', data);
end

function data = updatePlot(data, initial)
    ax = data.ax;
    if initial || isempty(data.scatter) || ~isgraphics(data.scatter)
        cla(ax);
        hold(ax, 'on');
        [x, y, z] = samplePoints(data.posOut, data.sample);
        data.scatter = scatter3(ax, x, y, z, 8, z, 'filled');
        axis(ax, 'equal');
        view(ax, 3);
    else
        [x, y, z] = samplePoints(data.posOut, data.sample);
        data.scatter.XData = x;
        data.scatter.YData = y;
        data.scatter.ZData = z;
        data.scatter.CData = z;
    end

    if data.reconstructHull
        if ~isempty(data.hullPatch)
            delete(data.hullPatch(isgraphics(data.hullPatch)));
        end
        if ~isempty(data.hullRec)
            data.hullPatch = patch(ax, data.hullRec(1), 'FaceColor', [0.2 0.8 1], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
        end
    elseif ~isempty(data.hullPatch)
        delete(data.hullPatch(isgraphics(data.hullPatch)));
        data.hullPatch = [];
    end
end

function [x, y, z] = samplePoints(pos, sample)
    n = height(pos);
    if sample <= 0
        sample = 0.05;
    end
    if sample < 1
        idx = randperm(n, max(1, round(n * sample)));
    else
        idx = randperm(n, min(n, round(sample)));
    end
    x = pos.x(idx);
    y = pos.y(idx);
    z = pos.z(idx);
end

function field = detectVoltageField(pos)
    if ismember('VDC', pos.Properties.VariableNames)
        field = "VDC";
    elseif ismember('VP', pos.Properties.VariableNames)
        field = "VP";
    else
        field = "";
    end
end

function list = voltageFieldList(pos)
    list = {};
    if ismember('VDC', pos.Properties.VariableNames)
        list{end+1} = 'VDC';
    end
    if ismember('VP', pos.Properties.VariableNames)
        list{end+1} = 'VP';
    end
    if isempty(list)
        list = {'VDC'};
    end
end

function idx = voltageFieldIndex(pos, field)
    list = voltageFieldList(pos);
    idx = find(strcmp(list, char(field)), 1);
    if isempty(idx)
        idx = 1;
    end
end

function validatePos(pos)
    required = {'detx','dety'};
    missing = setdiff(required, pos.Properties.VariableNames);
    if ~isempty(missing)
        error('reconstructionWidget:missingColumns', ...
            'pos table missing columns: %s', strjoin(missing, ', '));
    end
end

function [fcn, fcnText] = fitRadiusFromImage()
    fcn = [];
    fcnText = '';
    [file, path] = uigetfile({'*.png;*.jpg;*.tif;*.bmp', 'Images'}, 'Select correlative image');
    if isequal(file, 0)
        return;
    end
    img = imread(fullfile(path, file));
    fig = figure('Name', 'Select points (press Enter when done)');
    imshow(img);
    title('Click points along radius evolution, press Enter to finish');
    [x, y] = ginput;
    close(fig);
    if numel(x) < 2
        return;
    end
    prompt = {'Ion index min', 'Ion index max', 'Radius min (nm)', 'Radius max (nm)'};
    def = {'1', '100000', '10', '200'};
    answ = inputdlg(prompt, 'Axis calibration', 1, def);
    if isempty(answ)
        return;
    end
    xMin = str2double(answ{1});
    xMax = str2double(answ{2});
    yMin = str2double(answ{3});
    yMax = str2double(answ{4});
    if any(~isfinite([xMin xMax yMin yMax]))
        return;
    end

    imgSize = size(img);
    xScaled = xMin + (x - 1) / max(1, imgSize(2)-1) * (xMax - xMin);
    yScaled = yMax - (y - 1) / max(1, imgSize(1)-1) * (yMax - yMin);

    [xScaled, sortIdx] = sort(xScaled);
    yScaled = yScaled(sortIdx);
    fcn = @(i) interp1(xScaled, yScaled, i, 'pchip', 'extrap');
    fcnText = '@(i) interp1(x, y, i, ''pchip'', ''extrap'')';
end
