function [spec, controlFig, info] = massSpecWidget(source, options)
% MASSSPECWIDGET Interactive widget for mass spectrum ion/range workflows.
%
% [spec, controlFig, info] = massSpecWidget(spec)
% [spec, controlFig, info] = massSpecWidget(axOrFig)
% [spec, controlFig, info] = massSpecWidget(posOrMc)
% [spec, controlFig, info] = massSpecWidget(fileName)
%
% INPUT
%   source     mass spectrum handle/axes/figure, pos table, mc vector,
%              or file path (.fig/.pos/.epos/.apt/.h5)
%
% OPTIONS
%   'binWidth'        bin width for newly created spectrum (default: 0.01)
%   'mode'            spectrum mode for new spectrum (default: "count")
%   'controlTitle'    widget title (default: "Mass Spectrum Controls")
%   'colorScheme'     optional color scheme table (ion, color)
%   'isotopeTable'    optional isotope table override
%   'ionList'         optional ion list for ionFind
%   'searchRange'     ionFind search range (default: 0.3)
%   'rangeMargin'     rangeAddAll view margin in Da (default: 1)
%   'useMinForAddAll' use threshold line in rangeAddAll (default: false)
%   'presetVarName'   default workspace profile variable name
%
% OUTPUT
%   spec        mass spectrum area handle
%   controlFig  control widget figure handle
%   info        struct with source/colorScheme context
%
% NOTES
%   This first version focuses on integrating existing toolbox operations:
%   ionAdd, ionFind, rangeAdd, rangeAddAll, rangeExtendToNeighbor,
%   ionsExtractFromMassSpec, rangesExtractFromMassSpec, and concentration
%   calculations via posCalculateConcentrationSimple,
%   posCalculateConcentrationBackgroundRemoved, and
%   posCalculateConcentrationDeconvolved.
%   Keyboard shortcuts are available (press F1 or ? in the widget).
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    source = []
    options.binWidth (1,1) double {mustBePositive} = 0.01
    options.mode (1,1) string = "count"
    options.controlTitle (1,1) string = "Mass Spectrum Controls"
    options.colorScheme = []
    options.isotopeTable = []
    options.ionList = []
    options.searchRange (1,1) double {mustBePositive} = 0.3
    options.rangeMargin (1,1) double {mustBePositive} = 1
    options.useMinForAddAll (1,1) logical = false
    options.presetVarName (1,1) string = "massSpecProfile"
end

[spec, sourceInfo] = resolveMassSpectrumSource(source, options);
ax = ancestor(spec, 'axes');
fig = ancestor(spec, 'figure');

data = struct();
data.spec = spec;
data.ax = ax;
data.fig = fig;
data.sourceInfo = sourceInfo;
data.colorScheme = resolveColorScheme(options.colorScheme);
data.colorScheme = ensureBackgroundColor(data.colorScheme);
data.isotopeTable = resolveIsotopeTable(options.isotopeTable);
data.ionList = resolveIonList(options.ionList);
data.searchRange = options.searchRange;
data.rangeMargin = options.rangeMargin;
data.useMinForAddAll = options.useMinForAddAll;
data.presetVarName = options.presetVarName;
data.binWidth = sourceInfo.binWidth;
if ~isfinite(data.binWidth) || data.binWidth <= 0
    data.binWidth = options.binWidth;
end
data.sourceMode = string(sourceInfo.mode);
if strlength(data.sourceMode) == 0
    data.sourceMode = options.mode;
end
data.selectedIonRow = [];
data.selectedRangeRow = [];
data.ionRows = struct('handle', {}, 'name', {}, 'baseName', {}, 'chargeState', {}, ...
    'isTracer', {}, 'color', {}, 'x', {}, 'xAll', {});
data.rangeRows = struct('handle', {}, 'name', {}, 'mcbegin', {}, 'mcend', {}, 'chargeState', {}, 'color', {});
data.extractedElements = strings(0, 1);
data.candidateIons = table();
data.candidateDisplay = strings(0, 1);
data.selectedCandidateRows = [];
data.previewIonHandle = gobjects(0);
data.previewCandidate = struct('ionName', "", 'chargeState', NaN);
data.compositionMethod = "simple";
data.compositionMode = "ionic";
data.compositionDetEff = 0.8;
data.compositionAllocation = "raw";
data.compositionExcludeText = "unranged";
data.compositionVolumeName = "none";
data.compositionPosVarName = "pos";
data.compositionOutputVar = "conc";
data.compositionInfoVar = "concInfo";
data.compositionBgMethod = "linearBetweenPeaks";
data.compositionBackgroundMethod = "none";
data.compositionPlotBackground = true;
data.compositionPlotFits = false;
data.compositionExtraOptions = "";
data.compositionUseSourcePos = isfield(sourceInfo, 'pos') && istable(sourceInfo.pos) && ~isempty(sourceInfo.pos);

controlFig = createControlWindow(options.controlTitle);
setappdata(controlFig, 'massSpecWidget', data);
setupControls(controlFig);
refreshUi(controlFig);

if isgraphics(fig)
    addlistener(fig, 'ObjectBeingDestroyed', @(~, ~) safeDelete(controlFig));
end

info = struct();
info.source = sourceInfo;
info.figure = fig;
info.axes = ax;
info.colorScheme = data.colorScheme;
info.binWidth = data.binWidth;

end

function [spec, sourceInfo] = resolveMassSpectrumSource(source, options)
sourceInfo = struct('type', "", 'path', "", 'createdNewSpectrum', false, ...
    'mc', [], 'pos', [], 'binWidth', NaN, 'mode', string(options.mode));

if isempty(source)
    source = gcf;
end

if isgraphics(source, 'matlab.graphics.chart.primitive.Area')
    spec = source;
    sourceInfo.type = "areaHandle";
    sourceInfo.binWidth = estimateBinWidthFromSpec(spec);
    return;
end

if isgraphics(source, 'axes')
    spec = findMassSpectrumInAxes(source);
    if isempty(spec)
        error('massSpecWidget:noMassSpectrum', ...
            'No mass spectrum area plot found in provided axes.');
    end
    sourceInfo.type = "axesHandle";
    sourceInfo.binWidth = estimateBinWidthFromSpec(spec);
    return;
end

if isgraphics(source, 'figure')
    axesList = findobj(source, 'Type', 'axes');
    if isempty(axesList)
        error('massSpecWidget:noAxes', 'No axes found in provided figure.');
    end
    spec = [];
    for i = 1:numel(axesList)
        spec = findMassSpectrumInAxes(axesList(i));
        if ~isempty(spec)
            break;
        end
    end
    if isempty(spec)
        error('massSpecWidget:noMassSpectrum', ...
            'No mass spectrum area plot found in provided figure.');
    end
    sourceInfo.type = "figureHandle";
    sourceInfo.binWidth = estimateBinWidthFromSpec(spec);
    return;
end

if istable(source)
    if ~ismember('mc', source.Properties.VariableNames)
        error('massSpecWidget:missingMc', ...
            'Provided table must contain an mc column.');
    end
    spec = massSpecPlot(source.mc, options.binWidth, char(options.mode));
    sourceInfo.type = "posTable";
    sourceInfo.createdNewSpectrum = true;
    sourceInfo.mc = source.mc(:);
    sourceInfo.pos = source;
    sourceInfo.binWidth = options.binWidth;
    return;
end

if isnumeric(source)
    mc = source(:);
    spec = massSpecPlot(mc, options.binWidth, char(options.mode));
    sourceInfo.type = "mcVector";
    sourceInfo.createdNewSpectrum = true;
    sourceInfo.mc = mc;
    sourceInfo.binWidth = options.binWidth;
    return;
end

if ischar(source) || (isstring(source) && isscalar(source))
    fileName = char(source);
    if ~isfile(fileName)
        error('massSpecWidget:fileNotFound', 'File not found: %s', fileName);
    end
    [~, ~, ext] = fileparts(fileName);
    ext = lower(ext);
    sourceInfo.path = string(fileName);

    switch ext
        case '.fig'
            fig = openfig(fileName, 'new', 'visible');
            axesList = findobj(fig, 'Type', 'axes');
            spec = [];
            for i = 1:numel(axesList)
                spec = findMassSpectrumInAxes(axesList(i));
                if ~isempty(spec)
                    break;
                end
            end
            if isempty(spec)
                error('massSpecWidget:noMassSpectrumInFig', ...
                    'No mass spectrum area plot found in .fig file.');
            end
            sourceInfo.type = "figFile";
            sourceInfo.binWidth = estimateBinWidthFromSpec(spec);

        case {'.pos', '.epos', '.apt', '.h5', '.hdf5'}
            pos = posLoad(fileName);
            spec = massSpecPlot(pos.mc, options.binWidth, char(options.mode));
            sourceInfo.type = "dataFile";
            sourceInfo.createdNewSpectrum = true;
            sourceInfo.mc = pos.mc(:);
            sourceInfo.pos = pos;
            sourceInfo.binWidth = options.binWidth;

        otherwise
            error('massSpecWidget:unsupportedFile', ...
                'Unsupported file type ''%s''. Use .fig, .pos, .epos, .apt, .h5, .hdf5.', ext);
    end
    return;
end

error('massSpecWidget:unsupportedSource', ...
    'Unsupported source type. Use mass spectrum handle, axes/figure, pos table, mc vector, or file path.');
end

function spec = findMassSpectrumInAxes(ax)
spec = [];
plots = ax.Children;
for i = 1:numel(plots)
    try
        if isfield(plots(i).UserData, 'plotType') && plots(i).UserData.plotType == "massSpectrum"
            spec = plots(i);
            return;
        end
    catch
    end
end
for i = 1:numel(plots)
    if isa(plots(i), 'matlab.graphics.chart.primitive.Area')
        spec = plots(i);
        return;
    end
end
end

function controlFig = createControlWindow(titleText)
controlFig = figure('Name', char(titleText), ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Color', [0.94 0.95 0.96], ...
    'Units', 'pixels', ...
    'Position', [110 70 980 930], ...
    'Resize', 'on', ...
    'DefaultUicontrolFontName', 'Helvetica', ...
    'DefaultUicontrolFontSize', 10, ...
    'WindowKeyPressFcn', @onKeyPress);
end

function setupControls(controlFig)
data = getappdata(controlFig, 'massSpecWidget');

bg = get(controlFig, 'Color');
panelBg = [0.98 0.985 0.99];
btnBg = [1 1 1];
accent = [0.16 0.33 0.53];

uicontrol(controlFig, 'Style', 'text', ...
    'String', 'Mass Spectrum Operations', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.972 0.94 0.022], ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'bold', ...
    'ForegroundColor', accent, ...
    'BackgroundColor', bg);

% Source panel
sourcePanel = uipanel(controlFig, 'Title', 'Source', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.91 0.66 0.06], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

sourceText = uicontrol(sourcePanel, 'Style', 'text', ...
    'String', sourceSummaryText(data), ...
    'Units', 'normalized', ...
    'Position', [0.02 0.18 0.30 0.64], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

uicontrol(sourcePanel, 'Style', 'text', ...
    'String', 'Bin [Da]', ...
    'Units', 'normalized', ...
    'Position', [0.33 0.20 0.09 0.60], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

binWidthEdit = uicontrol(sourcePanel, 'Style', 'edit', ...
    'String', num2str(data.binWidth, '%.5g'), ...
    'Units', 'normalized', ...
    'Position', [0.42 0.18 0.11 0.64], ...
    'TooltipString', 'Mass-spectrum histogram bin width in Da (press Enter to rebin)', ...
    'Callback', @(~, ~) onRebin(controlFig));

refreshBtn = uicontrol(sourcePanel, 'Style', 'pushbutton', ...
    'String', 'Refresh', ...
    'Units', 'normalized', ...
    'Position', [0.56 0.14 0.16 0.72], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Refresh ion/range tables and status (Shortcut: R)', ...
    'Callback', @(~, ~) refreshUi(controlFig));

hotkeyBtn = uicontrol(sourcePanel, 'Style', 'pushbutton', ...
    'String', 'Help', ...
    'Units', 'normalized', ...
    'Position', [0.74 0.14 0.16 0.72], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Show keyboard shortcuts (Shortcut: ? or F1)', ...
    'Callback', @(~, ~) showKeyboardHelp());

% Ion panel
ionPanel = uipanel(controlFig, 'Title', 'Ions', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.56 0.66 0.35], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

uicontrol(ionPanel, 'Style', 'text', ...
    'String', 'Ion', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.89 0.07 0.09], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

ionEdit = uicontrol(ionPanel, 'Style', 'edit', ...
    'String', 'Fe', ...
    'Units', 'normalized', ...
    'Position', [0.09 0.892 0.16 0.088], ...
    'TooltipString', 'Ion name to add (e.g., Fe, FeO, 56Fe)');

uicontrol(ionPanel, 'Style', 'text', ...
    'String', 'Charge(s)', ...
    'Units', 'normalized', ...
    'Position', [0.26 0.89 0.11 0.09], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

chargeEdit = uicontrol(ionPanel, 'Style', 'edit', ...
    'String', '1', ...
    'Units', 'normalized', ...
    'Position', [0.37 0.892 0.08 0.088], ...
    'TooltipString', 'Charge states, e.g. 1 or 1 2 3 or 1:3');

addIonBtn = uicontrol(ionPanel, 'Style', 'pushbutton', ...
    'String', 'Add Ion', ...
    'Units', 'normalized', ...
    'Position', [0.46 0.892 0.11 0.088], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Add ion from fields (Shortcut: I)', ...
    'Callback', @(~, ~) onAddIon(controlFig));

findIonBtn = uicontrol(ionPanel, 'Style', 'pushbutton', ...
    'String', 'Find Ion', ...
    'Units', 'normalized', ...
    'Position', [0.58 0.892 0.11 0.088], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Launch ionFind for interactive peak-to-ion assignment (Shortcut: F)', ...
    'Callback', @(~, ~) onFindIon(controlFig));

importIonsBtn = uicontrol(ionPanel, 'Style', 'pushbutton', ...
    'String', 'Import', ...
    'Units', 'normalized', ...
    'Position', [0.70 0.892 0.16 0.088], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', ['Import ions from another mass spectrum (.fig or open figure) and ', ...
    'scale stem heights to current spectrum'], ...
    'Callback', @(~, ~) onImportIons(controlFig));

removeIonBtn = uicontrol(ionPanel, 'Style', 'pushbutton', ...
    'String', 'Remove', ...
    'Units', 'normalized', ...
    'Position', [0.87 0.892 0.11 0.088], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Remove selected ion stem plot (Shortcut: Delete when ion selected)', ...
    'Callback', @(~, ~) onRemoveIon(controlFig));

ionTable = uitable(ionPanel, ...
    'Data', cell(0, 4), ...
    'ColumnName', {'Ion', 'Tracer', 'Peak m/c', 'Color'}, ...
    'RowName', {}, ...
    'ColumnEditable', [true true false true], ...
    'ColumnWidth', {140, 60, 170, 90}, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.86], ...
    'CellSelectionCallback', @(~, evd) onIonSelection(controlFig, evd), ...
    'CellEditCallback', @(~, evd) onIonTableEdit(controlFig, evd));

% Range panel
rangePanel = uipanel(controlFig, 'Title', 'Ranges', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.31 0.66 0.24], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

addRangeBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Add', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.84 0.10 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Add one range interactively (Shortcut: A)', ...
    'Callback', @(~, ~) onAddRange(controlFig));

addAllRangeBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Add All', ...
    'Units', 'normalized', ...
    'Position', [0.13 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Guide through all unranged ion peaks (Shortcut: G)', ...
    'Callback', @(~, ~) onAddAllRanges(controlFig));

extendRangeBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Extend', ...
    'Units', 'normalized', ...
    'Position', [0.26 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Extend clicked range to neighbor (Shortcut: X)', ...
    'Callback', @(~, ~) onExtendRange(controlFig));

removeRangeBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Remove', ...
    'Units', 'normalized', ...
    'Position', [0.39 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Remove selected range (Shortcut: Delete when range selected)', ...
    'Callback', @(~, ~) onRemoveRange(controlFig));

legendBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Legend', ...
    'Units', 'normalized', ...
    'Position', [0.52 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Show mass spectrum legend (Shortcut: L)', ...
    'Callback', @(~, ~) onShowLegend(controlFig));

reorderBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Reorder', ...
    'Units', 'normalized', ...
    'Position', [0.61 0.84 0.11 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Reorder mass spectrum overlays (Shortcut: O)', ...
    'Callback', @(~, ~) onReorder(controlFig));

zoomBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Zoom Full', ...
    'Units', 'normalized', ...
    'Position', [0.73 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Reset mass spectrum x/y limits (Shortcut: Z)', ...
    'Callback', @(~, ~) onZoomFull(controlFig));

zoomRangedBtn = uicontrol(rangePanel, 'Style', 'pushbutton', ...
    'String', 'Zoom Ions', ...
    'Units', 'normalized', ...
    'Position', [0.86 0.84 0.12 0.14], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Zoom from 0 Da to max(ion peak, range end) + 5 Da (Shortcut: V)', ...
    'Callback', @(~, ~) onZoomRanged(controlFig));

rangeTable = uitable(rangePanel, ...
    'Data', cell(0, 3), ...
    'ColumnName', {'Name', 'Begin', 'End'}, ...
    'RowName', {}, ...
    'ColumnEditable', [true true true], ...
    'ColumnWidth', {140, 75, 75}, ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.80], ...
    'CellSelectionCallback', @(~, evd) onRangeSelection(controlFig, evd), ...
    'CellEditCallback', @(~, evd) onRangeTableEdit(controlFig, evd));

% Element parsing panel
elementPanel = uipanel(controlFig, 'Title', 'Element Input', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.18 0.66 0.12], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

elementsInput = uicontrol(elementPanel, 'Style', 'edit', ...
    'Max', 1, ...
    'Min', 0, ...
    'String', 'Fe Cr Ni Mo Ti Al Co C Mn Si P S B Cu Ta', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.72 0.73 0.22], ...
    'HorizontalAlignment', 'left', ...
    'TooltipString', 'Paste composition/text and extract element symbols');

extractElementsBtn = uicontrol(elementPanel, 'Style', 'pushbutton', ...
    'String', 'Extract', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.45 0.16 0.22], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Extract valid chemical element symbols from text (Shortcut: E)', ...
    'Callback', @(~, ~) onExtractElements(controlFig));

uicontrol(elementPanel, 'Style', 'text', ...
    'String', 'Chargestate', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.20 0.16 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

elementChargeEdit = uicontrol(elementPanel, 'Style', 'edit', ...
    'String', '+ ++', ...
    'Units', 'normalized', ...
    'Position', [0.18 0.18 0.20 0.22], ...
    'TooltipString', 'Charge states for selected elements, e.g. + ++ or 1 2');

addElementsBtn = uicontrol(elementPanel, 'Style', 'pushbutton', ...
    'String', 'Add Elem', ...
    'Units', 'normalized', ...
    'Position', [0.20 0.45 0.18 0.22], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Add mono-atomic ions for selected elements', ...
    'Callback', @(~, ~) onAddSelectedElements(controlFig));

addComplexBtn = uicontrol(elementPanel, 'Style', 'pushbutton', ...
    'String', 'Build Pot', ...
    'Units', 'normalized', ...
    'Position', [0.40 0.45 0.18 0.22], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Build potential ion list from selected elements', ...
    'Callback', @(~, ~) onAddComplexFromElements(controlFig));

uicontrol(elementPanel, 'Style', 'text', ...
    'String', 'Ion Complexity', ...
    'Units', 'normalized', ...
    'Position', [0.40 0.20 0.18 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

complexityEdit = uicontrol(elementPanel, 'Style', 'edit', ...
    'String', '1 2', ...
    'Units', 'normalized', ...
    'Position', [0.58 0.18 0.17 0.22], ...
    'TooltipString', 'Ion complexity orders, e.g. 1 2');

elementsList = uicontrol(elementPanel, 'Style', 'listbox', ...
    'String', {}, ...
    'Max', 2, ...
    'Min', 0, ...
    'Units', 'normalized', ...
    'Position', [0.78 0.06 0.20 0.88], ...
    'BackgroundColor', [1 1 1]);

% Potential ions panel (right column)
candidatePanel = uipanel(controlFig, 'Title', 'Potential Ions', ...
    'Units', 'normalized', ...
    'Position', [0.71 0.18 0.26 0.79], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

buildCandidatesBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Build', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.94 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Create potential ion list from selected elements/charges/complexity', ...
    'Callback', @(~, ~) onBuildCandidateIons(controlFig));

previewCandidateBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Preview', ...
    'Units', 'normalized', ...
    'Position', [0.35 0.94 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Preview selected potential ion in mass spectrum', ...
    'Callback', @(~, ~) onPreviewCandidate(controlFig));

addCandidateBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Add', ...
    'Units', 'normalized', ...
    'Position', [0.67 0.94 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Confirm selected potential ion and add to ion list', ...
    'Callback', @(~, ~) onAddCandidate(controlFig));

findCandidateBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Find', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.88 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Run ionFind for the selected potential ion', ...
    'Callback', @(~, ~) onFindCandidate(controlFig));

importIonsRightBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Import', ...
    'Units', 'normalized', ...
    'Position', [0.35 0.88 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Import ions from open figure or .fig and scale to this spectrum', ...
    'Callback', @(~, ~) onImportIons(controlFig));

clearPreviewBtn = uicontrol(candidatePanel, 'Style', 'pushbutton', ...
    'String', 'Clear', ...
    'Units', 'normalized', ...
    'Position', [0.67 0.88 0.30 0.05], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Remove preview ion from the plot', ...
    'Callback', @(~, ~) clearPreviewIon(controlFig));

candidateTable = uitable(candidatePanel, ...
    'Data', cell(0, 2), ...
    'ColumnName', {'Label', 'm/c'}, ...
    'RowName', {}, ...
    'ColumnEditable', [false false], ...
    'ColumnWidth', {100, 210}, ...
    'Units', 'normalized', ...
    'Position', [0.03 0.03 0.94 0.83], ...
    'CellSelectionCallback', @(~, evd) onCandidateSelection(controlFig, evd), ...
    'TooltipString', 'Potential ions from elements. Select row to preview.');

% Profile/preset + composition panel
profilePanel = uipanel(controlFig, 'Title', 'Profiles && Composition', ...
    'Units', 'normalized', ...
    'Position', [0.03 0.01 0.94 0.16], ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', accent, ...
    'FontWeight', 'bold');

saveWsBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Save -> WS', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.72 0.12 0.24], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Save current profile to workspace (Shortcut: Ctrl/Cmd+S)', ...
    'Callback', @(~, ~) onSaveProfileWorkspace(controlFig));

loadWsBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Load <- WS', ...
    'Units', 'normalized', ...
    'Position', [0.15 0.72 0.12 0.24], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Load profile from workspace (Shortcut: Ctrl/Cmd+L)', ...
    'Callback', @(~, ~) onLoadProfileWorkspace(controlFig));

exportBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Export File', ...
    'Units', 'normalized', ...
    'Position', [0.28 0.72 0.12 0.24], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Export profile to file (Shortcut: Ctrl/Cmd+E)', ...
    'Callback', @(~, ~) onExportProfileFile(controlFig));

importBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Import File', ...
    'Units', 'normalized', ...
    'Position', [0.41 0.72 0.12 0.24], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Import profile from file (Shortcut: Ctrl/Cmd+I)', ...
    'Callback', @(~, ~) onImportProfileFile(controlFig));

applyDataBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Apply -> Data', ...
    'Units', 'normalized', ...
    'Position', [0.54 0.72 0.14 0.24], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Apply current profile to another data file or mass spectrum figure (Shortcut: P)', ...
    'Callback', @(~, ~) onApplyProfileToDataset(controlFig));

clearExistingCb = uicontrol(profilePanel, 'Style', 'checkbox', ...
    'String', 'Clear existing before load', ...
    'Units', 'normalized', ...
    'Position', [0.70 0.74 0.28 0.20], ...
    'BackgroundColor', panelBg, ...
    'Value', true);

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Method', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.46 0.04 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compMethodPopup = uicontrol(profilePanel, 'Style', 'popupmenu', ...
    'String', {'simple', 'backgroundRemoved', 'deconvolved'}, ...
    'Value', 1, ...
    'Units', 'normalized', ...
    'Position', [0.06 0.44 0.15 0.18], ...
    'TooltipString', 'Concentration method', ...
    'Callback', @(~, ~) onCompositionMethodChanged(controlFig));

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Mode', ...
    'Units', 'normalized', ...
    'Position', [0.22 0.46 0.04 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compModePopup = uicontrol(profilePanel, 'Style', 'popupmenu', ...
    'String', {'ionic', 'isotopic', 'atomic'}, ...
    'Value', 1, ...
    'Units', 'normalized', ...
    'Position', [0.26 0.44 0.12 0.18], ...
    'TooltipString', 'Output grouping mode');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'DetEff', ...
    'Units', 'normalized', ...
    'Position', [0.39 0.46 0.05 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compDetEffEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', num2str(data.compositionDetEff, '%.4g'), ...
    'Units', 'normalized', ...
    'Position', [0.44 0.44 0.07 0.18], ...
    'TooltipString', 'Detector efficiency (fraction or %)');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Alloc', ...
    'Units', 'normalized', ...
    'Position', [0.52 0.46 0.04 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compAllocPopup = uicontrol(profilePanel, 'Style', 'popupmenu', ...
    'String', {'raw', 'decompose'}, ...
    'Value', 1, ...
    'Units', 'normalized', ...
    'Position', [0.56 0.44 0.09 0.18], ...
    'TooltipString', 'Allocation mode passed to posAllocateRange');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Background', ...
    'Units', 'normalized', ...
    'Position', [0.66 0.46 0.07 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compBgMethodPopup = uicontrol(profilePanel, 'Style', 'popupmenu', ...
    'String', {'linearBetweenPeaks', 'massSpecInvSqrt', 'none'}, ...
    'Value', 1, ...
    'Units', 'normalized', ...
    'Position', [0.73 0.44 0.12 0.18], ...
    'TooltipString', 'Background method (used by background/deconvolved methods)');

runConcBtn = uicontrol(profilePanel, 'Style', 'pushbutton', ...
    'String', 'Calc Conc', ...
    'Units', 'normalized', ...
    'Position', [0.86 0.44 0.12 0.18], ...
    'BackgroundColor', btnBg, ...
    'TooltipString', 'Calculate concentration with current options (Shortcut: C)', ...
    'Callback', @(~, ~) onComputeConcentration(controlFig));

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Pos Var', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.22 0.05 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compPosVarEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', char(data.compositionPosVarName), ...
    'Units', 'normalized', ...
    'Position', [0.07 0.20 0.12 0.18], ...
    'TooltipString', 'Workspace variable name for pos table');

compUseSourcePosCb = uicontrol(profilePanel, 'Style', 'checkbox', ...
    'String', 'Use source pos', ...
    'Units', 'normalized', ...
    'Position', [0.20 0.20 0.14 0.18], ...
    'BackgroundColor', panelBg, ...
    'Value', logical(data.compositionUseSourcePos), ...
    'TooltipString', 'Use pos loaded with the source (if available)');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Exclude', ...
    'Units', 'normalized', ...
    'Position', [0.35 0.22 0.06 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compExcludeEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', char(data.compositionExcludeText), ...
    'Units', 'normalized', ...
    'Position', [0.41 0.20 0.13 0.18], ...
    'TooltipString', 'Excluded categories, separated by spaces/commas (e.g. unranged)');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Volume', ...
    'Units', 'normalized', ...
    'Position', [0.55 0.22 0.05 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compVolumeEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', char(data.compositionVolumeName), ...
    'Units', 'normalized', ...
    'Position', [0.60 0.20 0.08 0.18], ...
    'TooltipString', 'Volume name written to concentration table');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Out', ...
    'Units', 'normalized', ...
    'Position', [0.69 0.22 0.03 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compOutVarEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', char(data.compositionOutputVar), ...
    'Units', 'normalized', ...
    'Position', [0.72 0.20 0.06 0.18], ...
    'TooltipString', 'Workspace output variable for concentration table');

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Info', ...
    'Units', 'normalized', ...
    'Position', [0.79 0.22 0.03 0.16], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compInfoVarEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', char(data.compositionInfoVar), ...
    'Units', 'normalized', ...
    'Position', [0.82 0.20 0.06 0.18], ...
    'TooltipString', 'Workspace output variable for info struct');

compPlotBgCb = uicontrol(profilePanel, 'Style', 'checkbox', ...
    'String', 'Plot BG', ...
    'Units', 'normalized', ...
    'Position', [0.89 0.20 0.09 0.18], ...
    'BackgroundColor', panelBg, ...
    'Value', false, ...
    'TooltipString', 'Plot background estimate onto mass spectrum');

compPlotFitsCb = uicontrol(profilePanel, 'Style', 'checkbox', ...
    'String', 'Plot Fits', ...
    'Units', 'normalized', ...
    'Position', [0.89 0.06 0.09 0.12], ...
    'BackgroundColor', panelBg, ...
    'Value', false, ...
    'TooltipString', 'Plot deconvolution fit diagnostics');

statusText = uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Ready', ...
    'Units', 'normalized', ...
    'Position', [0.02 0.02 0.40 0.12], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg, ...
    'ForegroundColor', [0.15 0.15 0.15]);

uicontrol(profilePanel, 'Style', 'text', ...
    'String', 'Extra NV', ...
    'Units', 'normalized', ...
    'Position', [0.43 0.02 0.06 0.12], ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', panelBg);

compExtraOptsEdit = uicontrol(profilePanel, 'Style', 'edit', ...
    'String', '', ...
    'Units', 'normalized', ...
    'Position', [0.49 0.02 0.39 0.12], ...
    'HorizontalAlignment', 'left', ...
    'TooltipString', ['Additional name-value arguments, e.g. ', ...
    '''minPeakDistance'',0.2,''fitLimits'',[5 20;40 60]']);

ctrls = struct();
ctrls.sourceText = sourceText;
ctrls.binWidthEdit = binWidthEdit;
ctrls.refreshBtn = refreshBtn;
ctrls.hotkeyBtn = hotkeyBtn;
ctrls.ionEdit = ionEdit;
ctrls.chargeEdit = chargeEdit;
ctrls.addIonBtn = addIonBtn;
ctrls.findIonBtn = findIonBtn;
ctrls.importIonsBtn = importIonsBtn;
ctrls.removeIonBtn = removeIonBtn;
ctrls.ionTable = ionTable;
ctrls.rangeTable = rangeTable;
ctrls.addRangeBtn = addRangeBtn;
ctrls.addAllRangeBtn = addAllRangeBtn;
ctrls.extendRangeBtn = extendRangeBtn;
ctrls.removeRangeBtn = removeRangeBtn;
ctrls.legendBtn = legendBtn;
ctrls.reorderBtn = reorderBtn;
ctrls.zoomBtn = zoomBtn;
ctrls.zoomRangedBtn = zoomRangedBtn;
ctrls.elementsInput = elementsInput;
ctrls.extractElementsBtn = extractElementsBtn;
ctrls.elementChargeEdit = elementChargeEdit;
ctrls.addElementsBtn = addElementsBtn;
ctrls.addComplexBtn = addComplexBtn;
ctrls.complexityEdit = complexityEdit;
ctrls.elementsList = elementsList;
ctrls.candidatePanel = candidatePanel;
ctrls.buildCandidatesBtn = buildCandidatesBtn;
ctrls.previewCandidateBtn = previewCandidateBtn;
ctrls.addCandidateBtn = addCandidateBtn;
ctrls.findCandidateBtn = findCandidateBtn;
ctrls.importIonsRightBtn = importIonsRightBtn;
ctrls.clearPreviewBtn = clearPreviewBtn;
ctrls.candidateTable = candidateTable;
ctrls.saveWsBtn = saveWsBtn;
ctrls.loadWsBtn = loadWsBtn;
ctrls.exportBtn = exportBtn;
ctrls.importBtn = importBtn;
ctrls.applyDataBtn = applyDataBtn;
ctrls.clearExistingCb = clearExistingCb;
ctrls.compMethodPopup = compMethodPopup;
ctrls.compModePopup = compModePopup;
ctrls.compDetEffEdit = compDetEffEdit;
ctrls.compAllocPopup = compAllocPopup;
ctrls.compBgMethodPopup = compBgMethodPopup;
ctrls.compPosVarEdit = compPosVarEdit;
ctrls.compUseSourcePosCb = compUseSourcePosCb;
ctrls.compExcludeEdit = compExcludeEdit;
ctrls.compVolumeEdit = compVolumeEdit;
ctrls.compOutVarEdit = compOutVarEdit;
ctrls.compInfoVarEdit = compInfoVarEdit;
ctrls.compPlotBgCb = compPlotBgCb;
ctrls.compPlotFitsCb = compPlotFitsCb;
ctrls.compExtraOptsEdit = compExtraOptsEdit;
ctrls.runConcBtn = runConcBtn;
ctrls.statusText = statusText;
ctrls.sourcePanel = sourcePanel;
ctrls.ionPanel = ionPanel;
ctrls.rangePanel = rangePanel;
ctrls.elementPanel = elementPanel;
ctrls.profilePanel = profilePanel;

setappdata(controlFig, 'massSpecWidgetControls', ctrls);
onCompositionMethodChanged(controlFig);
configureControlWindowLayout(controlFig);
end

function refreshUi(controlFig)
if ~isgraphics(controlFig)
    return;
end

data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
if ~isgraphics(data.spec)
    setStatus(controlFig, 'Mass spectrum handle is no longer valid.');
    return;
end

ctrls.sourceText.String = sourceSummaryText(data);
if isfield(ctrls, 'binWidthEdit') && isgraphics(ctrls.binWidthEdit)
    ctrls.binWidthEdit.String = num2str(data.binWidth, '%.5g');
end
syncCompositionControls(controlFig, data, ctrls);

data.ionRows = collectIonRows(data.spec);
data.rangeRows = collectRangeRows(data.spec);
setappdata(controlFig, 'massSpecWidget', data);

ionData = cell(numel(data.ionRows), 4);
for i = 1:numel(data.ionRows)
    ionData{i, 1} = char(data.ionRows(i).name);
    ionData{i, 2} = data.ionRows(i).isTracer;
    ionData{i, 3} = char(formatMcVector(data.ionRows(i).xAll));
    ionData{i, 4} = sprintf('[%.2f %.2f %.2f]', data.ionRows(i).color);
end
ctrls.ionTable.Data = ionData;

rangeData = cell(numel(data.rangeRows), 3);
for i = 1:numel(data.rangeRows)
    rangeData{i, 1} = char(data.rangeRows(i).name);
    rangeData{i, 2} = data.rangeRows(i).mcbegin;
    rangeData{i, 3} = data.rangeRows(i).mcend;
end
ctrls.rangeTable.Data = rangeData;

if isfield(ctrls, 'candidateTable') && isgraphics(ctrls.candidateTable)
    ctab = candidateTableData(data);
    ctrls.candidateTable.Data = ctab;
end

setStatus(controlFig, sprintf('Ions: %d | Ranges: %d', numel(data.ionRows), numel(data.rangeRows)));
end

function rows = collectIonRows(spec)
rows = struct('handle', {}, 'name', {}, 'baseName', {}, 'chargeState', {}, ...
    'isTracer', {}, 'color', {}, 'x', {}, 'xAll', {});
ax = ancestor(spec, 'axes');
plots = ax.Children;
for i = 1:numel(plots)
    h = plots(i);
    if ~isgraphics(h)
        continue;
    end
    try
        if ~isfield(h.UserData, 'plotType') || h.UserData.plotType ~= "ion"
            continue;
        end
    catch
        continue;
    end

    r = struct();
    r.handle = h;
    if isprop(h, 'DisplayName')
        r.name = string(h.DisplayName);
    else
        r.name = "ion";
    end
    r.baseName = "";
    try
        [baseName, ~, ~] = ionIdentityFromHandle(h);
        r.baseName = baseName;
    catch
    end
    if strlength(r.baseName) == 0
        r.baseName = canonicalIonBaseName(r.name);
    end
    r.chargeState = NaN;
    try
        cs = h.UserData.chargeState;
        if isnumeric(cs)
            r.chargeState = cs(1);
        end
    catch
    end
    r.isTracer = false;
    try
        if isfield(h.UserData, 'isTracer')
            r.isTracer = logical(h.UserData.isTracer);
        end
    catch
    end
    r.color = [0.5 0.5 0.5];
    try
        if isprop(h, 'Color')
            c = h.Color;
            if isnumeric(c) && numel(c) == 3
                r.color = c;
            end
        end
    catch
    end
    r.x = NaN;
    r.xAll = [];
    try
        if ~isempty(h.XData)
            xVals = double(h.XData(:)');
            xVals = xVals(isfinite(xVals));
            if ~isempty(xVals)
                xVals = sort(unique(xVals), 'ascend');
                r.xAll = xVals;
                r.x = xVals(1);
            end
        end
    catch
    end
    rows(end+1) = r; %#ok<AGROW>
end

if isempty(rows)
    return;
end
[~, idx] = sort([rows.x], 'ascend');
rows = rows(idx);
end

function rows = collectRangeRows(spec)
rows = struct('handle', {}, 'name', {}, 'mcbegin', {}, 'mcend', {}, 'chargeState', {}, 'color', {});
ax = ancestor(spec, 'axes');
plots = ax.Children;
for i = 1:numel(plots)
    h = plots(i);
    if ~isgraphics(h)
        continue;
    end
    try
        if ~isfield(h.UserData, 'plotType') || h.UserData.plotType ~= "range"
            continue;
        end
    catch
        continue;
    end

    r = struct();
    r.handle = h;
    r.name = "range";
    r.mcbegin = NaN;
    r.mcend = NaN;
    try
        if ~isempty(h.XData)
            r.mcbegin = h.XData(1);
            r.mcend = h.XData(end);
        end
    catch
    end
    r.chargeState = NaN;
    try
        if isfield(h.UserData, 'chargeState')
            r.chargeState = h.UserData.chargeState;
        end
    catch
    end
    % Prefer isotopic ion naming from range metadata.
    try
        ud = h.UserData;
        isTracer = false;
        if isstruct(ud) && isfield(ud, 'isTracer')
            isTracer = logical(ud.isTracer);
        end
        if isstruct(ud) && isfield(ud, 'ion')
            ionVal = ud.ion;
            if iscell(ionVal) && ~isempty(ionVal)
                ionVal = ionVal{1};
            end
            if (istable(ionVal) || iscategorical(ionVal)) && isfinite(r.chargeState)
                r.name = string(ionConvertName(ionVal, r.chargeState, 'plain', isTracer));
            elseif ischar(ionVal) || isstring(ionVal)
                r.name = string(ionVal);
            end
        end
    catch
    end
    if strlength(r.name) == 0 || r.name == "range"
        try
            if isprop(h, 'DisplayName')
                r.name = string(h.DisplayName);
            end
        catch
        end
    end
    r.color = [0.8 0.8 0.8];
    try
        if isprop(h, 'FaceColor') && isnumeric(h.FaceColor) && numel(h.FaceColor) == 3
            r.color = h.FaceColor;
        end
    catch
    end
    rows(end+1) = r; %#ok<AGROW>
end

if isempty(rows)
    return;
end
[~, idx] = sort([rows.mcbegin], 'ascend');
rows = rows(idx);
end

function onIonSelection(controlFig, evd)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
if isempty(evd.Indices)
    data.selectedIonRow = [];
else
    data.selectedIonRow = evd.Indices(1);
end
setappdata(controlFig, 'massSpecWidget', data);
end

function onRangeSelection(controlFig, evd)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
if isempty(evd.Indices)
    data.selectedRangeRow = [];
else
    data.selectedRangeRow = evd.Indices(1);
end
setappdata(controlFig, 'massSpecWidget', data);
end

function onIonTableEdit(controlFig, evd)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls) || isempty(evd.Indices)
    return;
end

row = evd.Indices(1);
col = evd.Indices(2);
if row < 1 || row > numel(data.ionRows)
    refreshUi(controlFig);
    return;
end

tbl = ctrls.ionTable.Data;
if size(tbl, 1) < row || size(tbl, 2) < 4
    refreshUi(controlFig);
    return;
end

switch col
    case {1, 2}
        [ok, msg, data] = replaceIonFromTableRow(data, tbl, row);
        setappdata(controlFig, 'massSpecWidget', data);
        refreshUi(controlFig);
        if ok
            setStatus(controlFig, msg);
        else
            setStatus(controlFig, "Ion edit failed: " + msg);
        end

    case 4
        [c, ok] = parseColorTriplet(tbl{row, 4});
        if ~ok
            refreshUi(controlFig);
            setStatus(controlFig, 'Invalid ion color format. Use [r g b] with values in [0,1].');
            return;
        end
        h = data.ionRows(row).handle;
        if ~isgraphics(h)
            refreshUi(controlFig);
            setStatus(controlFig, 'Selected ion handle is no longer valid.');
            return;
        end
        try
            h.Color = c;
        catch
        end

        ionName = canonicalIonBaseName(string(tbl{row, 1}));
        if strlength(ionName) > 0
            data.colorScheme = setColorSchemeIonColor(data.colorScheme, ionName, c);
        end
        setappdata(controlFig, 'massSpecWidget', data);
        refreshUi(controlFig);
        setStatus(controlFig, 'Ion color updated.');
end
end

function onRangeTableEdit(controlFig, evd)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls) || isempty(evd.Indices)
    return;
end

row = evd.Indices(1);
if row < 1 || row > numel(data.rangeRows)
    refreshUi(controlFig);
    return;
end

tbl = ctrls.rangeTable.Data;
if size(tbl, 1) < row || size(tbl, 2) < 3
    refreshUi(controlFig);
    return;
end

[ok, msg, data] = replaceRangeFromTableRow(data, tbl, row);
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
if ok
    setStatus(controlFig, msg);
else
    setStatus(controlFig, "Range edit failed: " + msg);
end
end

function onAddIon(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
ionName = strtrim(string(ctrls.ionEdit.String));
if strlength(ionName) == 0
    setStatus(controlFig, 'Ion name is empty.');
    return;
end
chargeStates = parseChargeStates(ctrls.chargeEdit.String, 1);

added = 0;
skipped = 0;
failed = 0;
for i = 1:numel(chargeStates)
    if ionExistsInSpec(data.spec, ionName, chargeStates(i))
        skipped = skipped + 1;
        continue;
    end
    try
        data.colorScheme = ensureIonColor(data.colorScheme, ionName);
        ionAdd(data.spec, char(ionName), chargeStates(i), data.isotopeTable, data.colorScheme, ...
            0, 0.01, 'most abundant', 0.1);
        added = added + 1;
    catch ME
        warning('massSpecWidget:addIonFailed', 'Add ion failed (%s): %s', ionName, ME.message);
        failed = failed + 1;
    end
end
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, sprintf('Added %d ion(s) for %s, skipped %d duplicate(s), failed %d.', ...
    added, ionName, skipped, failed));
end

function onFindIon(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
if exist('ionFind', 'file') == 0
    setStatus(controlFig, 'ionFind not found on MATLAB path.');
    return;
end
try
    ionFind(data.spec, data.ionList, data.isotopeTable, data.colorScheme, data.searchRange);
    refreshUi(controlFig);
    setStatus(controlFig, 'ionFind completed.');
catch ME
    warning('massSpecWidget:ionFindFailed', 'ionFind failed: %s', ME.message);
    setStatus(controlFig, 'ionFind failed. See command window warning.');
end
end

function onRemoveIon(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || isempty(data.selectedIonRow)
    setStatus(controlFig, 'Select an ion row first.');
    return;
end
idx = data.selectedIonRow;
if idx < 1 || idx > numel(data.ionRows)
    setStatus(controlFig, 'Selected ion row is out of range.');
    return;
end
h = data.ionRows(idx).handle;
if isgraphics(h)
    delete(h);
end
data.selectedIonRow = [];
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, 'Ion stem plot removed.');
end

function onAddRange(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    [hNew, ~] = rangeAdd(data.spec, data.colorScheme);
    if isgraphics(hNew)
        [rangeName, mcbegin, mcend, chargeState] = rangeIdentityFromHandle(hNew);
        if rangeExistsInSpec(data.spec, rangeName, mcbegin, mcend, chargeState, hNew)
            deleteRangeWithLabel(data.ax, hNew, rangeName);
            refreshUi(controlFig);
            setStatus(controlFig, 'Range already enlisted; duplicate ignored.');
            return;
        end
    end
    refreshUi(controlFig);
    setStatus(controlFig, 'Range added.');
catch ME
    warning('massSpecWidget:addRangeFailed', 'Add range failed: %s', ME.message);
    setStatus(controlFig, 'Range add canceled or failed.');
end
end

function onAddAllRanges(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    rangeAddAll(data.spec, data.colorScheme, data.rangeMargin, data.useMinForAddAll);
    refreshUi(controlFig);
    setStatus(controlFig, 'rangeAddAll completed.');
catch ME
    warning('massSpecWidget:addAllRangesFailed', 'rangeAddAll failed: %s', ME.message);
    setStatus(controlFig, 'rangeAddAll failed.');
end
end

function onExtendRange(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    rangeExtendToNeighbor(data.spec);
    refreshUi(controlFig);
    setStatus(controlFig, 'Range extension completed.');
catch ME
    warning('massSpecWidget:extendRangeFailed', 'Range extension failed: %s', ME.message);
    setStatus(controlFig, 'Range extension failed.');
end
end

function onRemoveRange(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || isempty(data.selectedRangeRow)
    setStatus(controlFig, 'Select a range row first.');
    return;
end
idx = data.selectedRangeRow;
if idx < 1 || idx > numel(data.rangeRows)
    setStatus(controlFig, 'Selected range row is out of range.');
    return;
end
h = data.rangeRows(idx).handle;
name = "";
if isgraphics(h)
    try
        name = string(h.DisplayName);
    catch
    end
end
deleteRangeWithLabel(data.ax, h, name);

data.selectedRangeRow = [];
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, 'Range removed.');
end

function onShowLegend(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    massSpecShowLegend(data.spec, ["ion", "range", "massSpectrum"]);
    setStatus(controlFig, 'Legend shown.');
catch ME
    warning('massSpecWidget:legendFailed', 'Show legend failed: %s', ME.message);
    setStatus(controlFig, 'Show legend failed.');
end
end

function onReorder(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    massSpecReorderPlot(data.spec);
    setStatus(controlFig, 'Plot order updated.');
catch ME
    warning('massSpecWidget:reorderFailed', 'Reorder failed: %s', ME.message);
    setStatus(controlFig, 'Reorder failed.');
end
end

function onZoomFull(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || ~isgraphics(data.ax)
    return;
end
try
    if isgraphics(data.spec)
        data.ax.XLim = [min(data.spec.XData), max(data.spec.XData)];
    end
    data.ax.YLimMode = 'auto';
    setStatus(controlFig, 'Zoom reset.');
catch
    setStatus(controlFig, 'Zoom reset failed.');
end
end

function onZoomRanged(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || ~isgraphics(data.ax)
    return;
end
rangeRows = collectRangeRows(data.spec);
ionRows = collectIonRows(data.spec);

mcEnd = [];
if ~isempty(rangeRows)
    mcEnd = [mcEnd, [rangeRows.mcend]]; %#ok<AGROW>
end
if ~isempty(ionRows)
    mcEnd = [mcEnd, [ionRows.x]]; %#ok<AGROW>
end
mcEnd = mcEnd(isfinite(mcEnd));
if isempty(mcEnd)
    setStatus(controlFig, 'No ions or ranges available for Zoom Ions.');
    return;
end
try
    data.ax.XLim = [0, max(mcEnd) + 5];
    data.ax.YLimMode = 'auto';
    setStatus(controlFig, sprintf('Zoomed to ion/range region [0, %.2f] Da.', max(mcEnd) + 5));
catch
    setStatus(controlFig, 'Zoom Ions failed.');
end
end

function onRebin(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
newBin = parseFiniteScalar(ctrls.binWidthEdit.String, data.binWidth);
if ~isfinite(newBin) || newBin <= 0
    setStatus(controlFig, 'Invalid bin width. Enter a positive numeric value.');
    return;
end
if isempty(data.sourceInfo.mc)
    setStatus(controlFig, 'Rebin not available: source has no raw m/c data.');
    return;
end

try
    profile = massSpecProfileFromSpec(data.spec, data.colorScheme);
catch ME
    warning('massSpecWidget:rebinProfileCaptureFailed', ...
        'Could not capture current profile before rebin: %s', ME.message);
    setStatus(controlFig, 'Could not capture current profile before rebin.');
    return;
end

tmpFig = [];
try
    specNew = massSpecPlot(data.sourceInfo.mc, newBin, char(data.sourceMode));
    tmpFig = ancestor(specNew, 'figure');
    xNew = specNew.XData;
    yNew = specNew.YData;

    if isgraphics(data.spec)
        data.spec.XData = xNew;
        data.spec.YData = yNew;
    else
        error('Current mass spectrum handle is no longer valid.');
    end

    data.ax = ancestor(data.spec, 'axes');
    data.fig = ancestor(data.spec, 'figure');
    data.sourceInfo.binWidth = newBin;
    data.binWidth = newBin;

    [~, colorSchemeOut] = massSpecProfileApply(data.spec, profile, ...
        'colorScheme', data.colorScheme, ...
        'isotopeTable', data.isotopeTable, ...
        'clearExisting', true);
    data.colorScheme = colorSchemeOut;

    setappdata(controlFig, 'massSpecWidget', data);
    refreshUi(controlFig);
    setStatus(controlFig, sprintf('Rebinned to %.5g Da and reapplied ions/ranges.', newBin));
catch ME
    warning('massSpecWidget:rebinFailed', 'Rebin failed: %s', ME.message);
    setStatus(controlFig, 'Rebin failed. See warning in command window.');
end

if ~isempty(tmpFig) && isgraphics(tmpFig)
    delete(tmpFig);
end
end

function onExtractElements(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
textIn = readMultilineEdit(ctrls.elementsInput);
elements = extractElementsRobust(textIn);
data.extractedElements = elements;
setappdata(controlFig, 'massSpecWidget', data);
if isempty(elements)
    ctrls.elementsList.String = {};
    ctrls.elementsList.Value = [];
    setStatus(controlFig, 'No element symbols detected.');
else
    ctrls.elementsList.String = cellstr(elements);
    ctrls.elementsList.Value = 1:min(1, numel(elements));
    setStatus(controlFig, sprintf('Extracted %d element(s).', numel(elements)));
end
end

function onBuildCandidateIons(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end

items = string(ctrls.elementsList.String);
sel = ctrls.elementsList.Value;
if isempty(items) || isempty(sel)
    setStatus(controlFig, 'Extract/select elements first.');
    return;
end
sel = sel(sel >= 1 & sel <= numel(items));
elements = unique(items(sel), 'stable');
chargeStates = parseChargeStates(ctrls.elementChargeEdit.String, [1 2]);
complexity = parsePositiveVector(ctrls.complexityEdit.String, [1 2]);

try
    ionList = ionsCreateComplex(cellstr(elements), complexity, data.isotopeTable, chargeStates);
catch ME
    warning('massSpecWidget:buildCandidateIonsFailed', ...
        'Could not build potential ion list: %s', ME.message);
    setStatus(controlFig, 'Could not build potential ion list.');
    return;
end

clearPreviewIon(controlFig);
data = getappdata(controlFig, 'massSpecWidget');
data.candidateIons = summarizeCandidateIons(ionList);
data.selectedCandidateRows = [];
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, sprintf('Built %d potential ions.', height(data.candidateIons)));
end

function onCandidateSelection(controlFig, evd)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
if isempty(evd.Indices)
    data.selectedCandidateRows = [];
    setappdata(controlFig, 'massSpecWidget', data);
    clearPreviewIon(controlFig);
    return;
end
rows = unique(evd.Indices(:,1), 'stable');
rows = rows(rows >= 1 & rows <= height(data.candidateIons));
if isempty(data.candidateIons) || isempty(rows)
    data.selectedCandidateRows = [];
    setappdata(controlFig, 'massSpecWidget', data);
    clearPreviewIon(controlFig);
    return;
end
data.selectedCandidateRows = rows(:)';
setappdata(controlFig, 'massSpecWidget', data);
onPreviewCandidate(controlFig);
end

function onPreviewCandidate(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || isempty(data.candidateIons) || isempty(data.selectedCandidateRows)
    clearPreviewIon(controlFig);
    setStatus(controlFig, 'Select one or more potential ion rows first.');
    return;
end
rows = data.selectedCandidateRows;
rows = rows(rows >= 1 & rows <= height(data.candidateIons));
if isempty(rows)
    clearPreviewIon(controlFig);
    return;
end

clearPreviewIon(controlFig);
previewHandles = gobjects(0);
previewLabels = strings(0, 1);
try
    for i = 1:numel(rows)
        row = rows(i);
        ionName = string(data.candidateIons.ionName(row));
        chargeState = double(data.candidateIons.chargeState(row));
        data.colorScheme = ensureIonColor(data.colorScheme, ionName);
        h = ionAdd(data.spec, char(ionName), chargeState, data.isotopeTable, data.colorScheme, ...
            0, 0.01, 'most abundant', 0.1);
        if isgraphics(h)
            h.LineStyle = ':';
            h.LineWidth = 1.5;
            h.Marker = 'none';
            ud = h.UserData;
            if ~isstruct(ud)
                ud = struct();
            end
            ud.plotType = "ionPreview";
            h.UserData = ud;
            previewHandles(end+1, 1) = h; %#ok<AGROW>
            previewLabels(end+1, 1) = ionName + repmat('+', 1, chargeState); %#ok<AGROW>
        end
    end
    data.previewIonHandle = previewHandles;
    if numel(previewLabels) == 1
        data.previewCandidate.ionName = previewLabels(1);
    else
        data.previewCandidate.ionName = strjoin(previewLabels, ", ");
    end
    data.previewCandidate.chargeState = NaN;
    setappdata(controlFig, 'massSpecWidget', data);
    setStatus(controlFig, sprintf('Previewing %d selected potential ion(s).', numel(previewLabels)));
catch ME
    warning('massSpecWidget:previewCandidateFailed', ...
        'Preview ion failed: %s', ME.message);
    setStatus(controlFig, 'Preview failed.');
end
end

function onAddCandidate(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || isempty(data.candidateIons) || isempty(data.selectedCandidateRows)
    setStatus(controlFig, 'Select one or more potential ion rows first.');
    return;
end
rows = data.selectedCandidateRows;
rows = rows(rows >= 1 & rows <= height(data.candidateIons));
if isempty(rows)
    return;
end
clearPreviewIon(controlFig);
data = getappdata(controlFig, 'massSpecWidget');
added = 0;
skipped = 0;
failed = 0;
try
    for i = 1:numel(rows)
        row = rows(i);
        ionName = string(data.candidateIons.ionName(row));
        chargeState = double(data.candidateIons.chargeState(row));
        if ionExistsInSpec(data.spec, ionName, chargeState)
            skipped = skipped + 1;
            continue;
        end
        data.colorScheme = ensureIonColor(data.colorScheme, ionName);
        try
            ionAdd(data.spec, char(ionName), chargeState, data.isotopeTable, data.colorScheme, ...
                0, 0.01, 'most abundant', 0.1);
            added = added + 1;
        catch MEi
            warning('massSpecWidget:addCandidateIonFailed', ...
                'Add candidate ion failed (%s, %d+): %s', ionName, chargeState, MEi.message);
            failed = failed + 1;
        end
    end
    setappdata(controlFig, 'massSpecWidget', data);
    refreshUi(controlFig);
    setStatus(controlFig, sprintf('Added %d, skipped %d duplicate(s), failed %d.', added, skipped, failed));
catch ME
    warning('massSpecWidget:addCandidateFailed', ...
        'Add candidate ion failed: %s', ME.message);
    setStatus(controlFig, 'Add candidate ion failed.');
end
end

function onFindCandidate(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data) || isempty(data.candidateIons) || isempty(data.selectedCandidateRows)
    setStatus(controlFig, 'Select at least one potential ion row first.');
    return;
end
if exist('ionFind', 'file') == 0
    setStatus(controlFig, 'ionFind not found on MATLAB path.');
    return;
end
row = data.selectedCandidateRows(1);
ionName = string(data.candidateIons.ionName(row));
chargeState = double(data.candidateIons.chargeState(row));
try
    subList = buildCandidateSubIonList(ionName, chargeState, data.isotopeTable);
    ionFind(data.spec, subList, data.isotopeTable, data.colorScheme, data.searchRange);
    refreshUi(controlFig);
    setStatus(controlFig, 'ionFind executed for selected potential ion.');
catch ME
    warning('massSpecWidget:findCandidateFailed', ...
        'ionFind for selected candidate failed: %s', ME.message);
    setStatus(controlFig, 'Find selected candidate failed.');
end
end

function onAddComplexFromElements(controlFig)
onBuildCandidateIons(controlFig);
end

function onAddSelectedElements(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
items = string(ctrls.elementsList.String);
if isempty(items)
    setStatus(controlFig, 'No extracted elements available.');
    return;
end
if isempty(ctrls.elementsList.Value)
    setStatus(controlFig, 'Select at least one element.');
    return;
end
selIdx = ctrls.elementsList.Value;
selIdx = selIdx(selIdx >= 1 & selIdx <= numel(items));
if isempty(selIdx)
    setStatus(controlFig, 'Select at least one element.');
    return;
end
chargeStates = parseChargeStates(ctrls.elementChargeEdit.String, [1 2]);

added = 0;
skipped = 0;
failed = 0;
for i = 1:numel(selIdx)
    el = items(selIdx(i));
    for cs = chargeStates
        if ionExistsInSpec(data.spec, el, cs)
            skipped = skipped + 1;
            continue;
        end
        try
            data.colorScheme = ensureIonColor(data.colorScheme, el);
            ionAdd(data.spec, char(el), cs, data.isotopeTable, data.colorScheme, ...
                0, 0.01, 'most abundant', 0.1);
            added = added + 1;
        catch ME
            warning('massSpecWidget:addElementIonFailed', ...
                'Add ion failed (%s, %d+): %s', el, cs, ME.message);
            failed = failed + 1;
        end
    end
end
setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, sprintf('Added %d ion(s), skipped %d duplicate(s), failed %d.', added, skipped, failed));
end

function onSaveProfileWorkspace(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    profile = massSpecProfileFromSpec(data.spec, data.colorScheme);
    varName = data.presetVarName;
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
        ts = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
        varName = data.presetVarName + "_" + ts;
    end
    assignin('base', char(varName), profile);
    setStatus(controlFig, sprintf('Saved profile to workspace variable "%s".', varName));
catch ME
    warning('massSpecWidget:saveProfileWsFailed', 'Save profile failed: %s', ME.message);
    setStatus(controlFig, 'Save profile to workspace failed.');
end
end

function onLoadProfileWorkspace(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end

varName = selectWorkspaceProfileVariable(data.presetVarName);
if strlength(varName) == 0
    return;
end

try
    profile = evalin('base', char(varName));
    clearExisting = logical(ctrls.clearExistingCb.Value);
    [report, colorSchemeOut] = massSpecProfileApply(data.spec, profile, ...
        'colorScheme', data.colorScheme, ...
        'isotopeTable', data.isotopeTable, ...
        'clearExisting', clearExisting);
    data.colorScheme = colorSchemeOut;
    setappdata(controlFig, 'massSpecWidget', data);
    refreshUi(controlFig);
    setStatus(controlFig, sprintf('Loaded profile "%s" (ions: %d, ranges: %d).', ...
        varName, report.ionsAdded, report.rangesAdded));
catch ME
    warning('massSpecWidget:loadProfileWsFailed', 'Load profile failed: %s', ME.message);
    setStatus(controlFig, 'Load profile from workspace failed.');
end
end

function onExportProfileFile(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
try
    profile = massSpecProfileFromSpec(data.spec, data.colorScheme);
    [file, path, filterIdx] = uiputfile( ...
        {'*.yaml', 'YAML (*.yaml)'; '*.yml', 'YAML (*.yml)'; '*.json', 'JSON (*.json)'; '*.mat', 'MAT (*.mat)'}, ...
        'Export mass spectrum profile', ...
        'massSpecProfile.yaml');
    if isequal(file, 0)
        return;
    end
    outFile = fullfile(path, file);
    if isempty(fileparts(outFile))
        switch filterIdx
            case 1
                outFile = [outFile '.yaml'];
            case 2
                outFile = [outFile '.yml'];
            case 3
                outFile = [outFile '.json'];
            otherwise
                outFile = [outFile '.mat'];
        end
    end
    configExport(profile, outFile);
    setStatus(controlFig, sprintf('Profile exported: %s', outFile));
catch ME
    warning('massSpecWidget:exportProfileFailed', 'Export profile failed: %s', ME.message);
    setStatus(controlFig, 'Profile export failed.');
end
end

function onImportProfileFile(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end
try
    [file, path] = uigetfile( ...
        {'*.yaml;*.yml;*.json;*.mat', 'Profile files (*.yaml,*.yml,*.json,*.mat)'}, ...
        'Import mass spectrum profile');
    if isequal(file, 0)
        return;
    end
    inFile = fullfile(path, file);
    profile = configImport(inFile);
    clearExisting = logical(ctrls.clearExistingCb.Value);
    [report, colorSchemeOut] = massSpecProfileApply(data.spec, profile, ...
        'colorScheme', data.colorScheme, ...
        'isotopeTable', data.isotopeTable, ...
        'clearExisting', clearExisting);
    data.colorScheme = colorSchemeOut;
    setappdata(controlFig, 'massSpecWidget', data);
    refreshUi(controlFig);
    setStatus(controlFig, sprintf('Profile imported (ions: %d, ranges: %d).', ...
        report.ionsAdded, report.rangesAdded));
catch ME
    warning('massSpecWidget:importProfileFailed', 'Import profile failed: %s', ME.message);
    setStatus(controlFig, 'Profile import failed.');
end
end

function onApplyProfileToDataset(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end

try
    profile = massSpecProfileFromSpec(data.spec, data.colorScheme);
catch ME
    warning('massSpecWidget:applyToDatasetProfileCreateFailed', ...
        'Could not create profile from current mass spectrum: %s', ME.message);
    setStatus(controlFig, 'Failed to capture current profile.');
    return;
end

[file, path] = uigetfile( ...
    {'*.fig;*.pos;*.epos;*.apt;*.h5;*.hdf5', 'Mass spectrum/data (*.fig,*.pos,*.epos,*.apt,*.h5,*.hdf5)'}, ...
    'Apply current profile to');
if isequal(file, 0)
    return;
end
targetFile = fullfile(path, file);

binWidth = estimateBinWidthFromSpec(data.spec);
opts = struct('binWidth', binWidth, 'mode', "count");
try
    [targetSpec, ~] = resolveMassSpectrumSource(targetFile, opts);
catch ME
    warning('massSpecWidget:applyToDatasetOpenFailed', ...
        'Could not open target data/figure: %s', ME.message);
    setStatus(controlFig, 'Failed to open selected target file.');
    return;
end

try
    [report, colorSchemeOut] = massSpecProfileApply(targetSpec, profile, ...
        'colorScheme', data.colorScheme, ...
        'isotopeTable', data.isotopeTable, ...
        'clearExisting', logical(ctrls.clearExistingCb.Value));
    data.colorScheme = colorSchemeOut;
    setappdata(controlFig, 'massSpecWidget', data);
    refreshUi(controlFig);
    [~, shortName, shortExt] = fileparts(targetFile);
    setStatus(controlFig, sprintf('Applied to %s%s (ions: %d, ranges: %d).', ...
        shortName, shortExt, report.ionsAdded, report.rangesAdded));
catch ME
    warning('massSpecWidget:applyToDatasetApplyFailed', ...
        'Could not apply profile to target: %s', ME.message);
    setStatus(controlFig, 'Profile apply to selected target failed.');
end
end

function onCompositionMethodChanged(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls) || ~isfield(ctrls, 'compMethodPopup')
    return;
end

method = lower(readPopupString(ctrls.compMethodPopup, data.compositionMethod));
if method == "backgroundremoved"
    bgChoice = lower(readPopupString(ctrls.compBgMethodPopup, data.compositionBackgroundMethod));
    if bgChoice == "none"
        setPopupString(ctrls.compBgMethodPopup, "linearBetweenPeaks");
    end
end
data.compositionMethod = method;
setappdata(controlFig, 'massSpecWidget', data);
applyCompositionMethodUiState(ctrls, method);
end

function onComputeConcentration(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(data) || isempty(ctrls)
    return;
end

data = readCompositionControlsIntoState(data, ctrls);
setappdata(controlFig, 'massSpecWidget', data);

if ~isgraphics(data.spec)
    setStatus(controlFig, 'Mass spectrum handle is not valid.');
    return;
end

try
    rangeTable = rangesExtractFromMassSpec(data.spec);
catch ME
    warning('massSpecWidget:extractRangesForConcFailed', ...
        'Could not extract ranges from mass spectrum: %s', ME.message);
    setStatus(controlFig, 'Could not extract ranges from mass spectrum.');
    return;
end
if isempty(rangeTable) || ~istable(rangeTable) || height(rangeTable) == 0
    setStatus(controlFig, 'No ranges found in mass spectrum. Add ranges first.');
    return;
end

[posRaw, posLabel, okPos, msgPos] = resolveCompositionPos(data);
if ~okPos
    setStatus(controlFig, msgPos);
    return;
end

try
    posForConc = posAllocateRange(posRaw, rangeTable, char(data.compositionAllocation));
catch ME
    warning('massSpecWidget:allocateForConcFailed', ...
        'posAllocateRange failed for concentration: %s', ME.message);
    setStatus(controlFig, 'Range allocation for concentration failed.');
    return;
end

excludeList = parseExcludeList(data.compositionExcludeText);
[extraOpts, okExtra, extraMsg] = parseExtraNameValueList(data.compositionExtraOptions);
if ~okExtra
    setStatus(controlFig, extraMsg);
    return;
end

mode = char(data.compositionMode);
volumeName = char(data.compositionVolumeName);
if isempty(strtrim(volumeName))
    volumeName = 'none';
end

method = lower(data.compositionMethod);
infoOut = [];
plotBackgroundFlag = logical(data.compositionPlotBackground);
if method == "backgroundremoved"
    plotBackgroundFlag = true;
elseif method == "deconvolved"
    bgForDeconvolution = lower(string(data.compositionBackgroundMethod));
    if bgForDeconvolution ~= "none" && bgForDeconvolution ~= ""
        plotBackgroundFlag = true;
    end
end
try
    switch method
        case "simple"
            conc = posCalculateConcentrationSimple(posForConc, data.compositionDetEff, ...
                excludeList, volumeName, 'mode', mode, extraOpts{:});

        case "backgroundremoved"
            conc = [];
            [conc, infoOut] = posCalculateConcentrationBackgroundRemoved(posForConc, ...
                data.compositionDetEff, excludeList, volumeName, ...
                'method', char(data.compositionBgMethod), ...
                'mode', mode, ...
                'massSpec', data.spec, ...
                'plotBackground', plotBackgroundFlag, ...
                extraOpts{:});

        case "deconvolved"
            conc = [];
            [conc, infoOut] = posCalculateConcentrationDeconvolved(posForConc, ...
                data.compositionDetEff, excludeList, volumeName, ...
                'mode', mode, ...
                'massSpec', data.spec, ...
                'plotBackground', plotBackgroundFlag, ...
                'plotFits', logical(data.compositionPlotFits), ...
                'backgroundMethod', char(data.compositionBackgroundMethod), ...
                extraOpts{:});

        otherwise
            setStatus(controlFig, sprintf('Unknown composition method: %s', method));
            return;
    end
catch ME
    warning('massSpecWidget:computeConcentrationFailed', ...
        'Concentration calculation failed (%s): %s', method, ME.message);
    setStatus(controlFig, sprintf('Concentration failed (%s). See warning.', method));
    return;
end

outVar = sanitizeWorkspaceVarName(data.compositionOutputVar, "conc");
assignin('base', char(outVar), conc);

infoVarAssigned = "";
if ~isempty(infoOut)
    infoVar = sanitizeWorkspaceVarName(data.compositionInfoVar, outVar + "Info");
    assignin('base', char(infoVar), infoOut);
    infoVarAssigned = infoVar;
end

numSpecies = max(0, width(conc) - 4);
if strlength(infoVarAssigned) > 0
    setStatus(controlFig, sprintf('Concentration (%s) -> %s, %s | %d species | %s', ...
        method, outVar, infoVarAssigned, numSpecies, posLabel));
else
    setStatus(controlFig, sprintf('Concentration (%s) -> %s | %d species | %s', ...
        method, outVar, numSpecies, posLabel));
end
end

function applyCompositionMethodUiState(ctrls, method)
if ~isfield(ctrls, 'compBgMethodPopup')
    return;
end
method = lower(string(method));
switch method
    case "simple"
        setEnable(ctrls.compBgMethodPopup, false);
        setEnable(ctrls.compPlotBgCb, false);
        setEnable(ctrls.compPlotFitsCb, false);
    case "backgroundremoved"
        setEnable(ctrls.compBgMethodPopup, true);
        setEnable(ctrls.compPlotBgCb, true);
        setEnable(ctrls.compPlotFitsCb, false);
    case "deconvolved"
        setEnable(ctrls.compBgMethodPopup, true);
        setEnable(ctrls.compPlotBgCb, true);
        setEnable(ctrls.compPlotFitsCb, true);
    otherwise
        setEnable(ctrls.compBgMethodPopup, true);
        setEnable(ctrls.compPlotBgCb, true);
        setEnable(ctrls.compPlotFitsCb, true);
end
end

function syncCompositionControls(controlFig, data, ctrls)
if isempty(ctrls) || ~isfield(ctrls, 'compMethodPopup')
    return;
end
setPopupString(ctrls.compMethodPopup, data.compositionMethod);
setPopupString(ctrls.compModePopup, data.compositionMode);
setPopupString(ctrls.compAllocPopup, data.compositionAllocation);
setPopupString(ctrls.compBgMethodPopup, data.compositionBackgroundMethod);

if isgraphics(ctrls.compDetEffEdit)
    ctrls.compDetEffEdit.String = num2str(data.compositionDetEff, '%.6g');
end
if isgraphics(ctrls.compPosVarEdit)
    ctrls.compPosVarEdit.String = char(data.compositionPosVarName);
end
if isgraphics(ctrls.compExcludeEdit)
    ctrls.compExcludeEdit.String = char(data.compositionExcludeText);
end
if isgraphics(ctrls.compVolumeEdit)
    ctrls.compVolumeEdit.String = char(data.compositionVolumeName);
end
if isgraphics(ctrls.compOutVarEdit)
    ctrls.compOutVarEdit.String = char(data.compositionOutputVar);
end
if isgraphics(ctrls.compInfoVarEdit)
    ctrls.compInfoVarEdit.String = char(data.compositionInfoVar);
end
if isgraphics(ctrls.compPlotBgCb)
    ctrls.compPlotBgCb.Value = logical(data.compositionPlotBackground);
end
if isgraphics(ctrls.compPlotFitsCb)
    ctrls.compPlotFitsCb.Value = logical(data.compositionPlotFits);
end
if isgraphics(ctrls.compExtraOptsEdit)
    ctrls.compExtraOptsEdit.String = char(data.compositionExtraOptions);
end

hasSourcePos = isfield(data.sourceInfo, 'pos') && istable(data.sourceInfo.pos) && ~isempty(data.sourceInfo.pos);
if isgraphics(ctrls.compUseSourcePosCb)
    ctrls.compUseSourcePosCb.Enable = onOff(hasSourcePos);
    ctrls.compUseSourcePosCb.Value = hasSourcePos && logical(data.compositionUseSourcePos);
end

applyCompositionMethodUiState(ctrls, data.compositionMethod);
end

function data = readCompositionControlsIntoState(data, ctrls)
data.compositionMethod = lower(readPopupString(ctrls.compMethodPopup, data.compositionMethod));
data.compositionMode = lower(readPopupString(ctrls.compModePopup, data.compositionMode));
data.compositionAllocation = lower(readPopupString(ctrls.compAllocPopup, data.compositionAllocation));

bgChoice = lower(readPopupString(ctrls.compBgMethodPopup, data.compositionBackgroundMethod));
if bgChoice == "none"
    data.compositionBgMethod = "linearBetweenPeaks";
else
    data.compositionBgMethod = bgChoice;
end
data.compositionBackgroundMethod = bgChoice;

data.compositionDetEff = parseFiniteScalar(ctrls.compDetEffEdit.String, data.compositionDetEff);
if ~isfinite(data.compositionDetEff) || data.compositionDetEff <= 0
    data.compositionDetEff = 0.8;
end

data.compositionPosVarName = strtrim(string(ctrls.compPosVarEdit.String));
if strlength(data.compositionPosVarName) == 0
    data.compositionPosVarName = "pos";
end
data.compositionUseSourcePos = logical(ctrls.compUseSourcePosCb.Value);
data.compositionExcludeText = strtrim(string(ctrls.compExcludeEdit.String));
data.compositionVolumeName = strtrim(string(ctrls.compVolumeEdit.String));
if strlength(data.compositionVolumeName) == 0
    data.compositionVolumeName = "none";
end

data.compositionOutputVar = strtrim(string(ctrls.compOutVarEdit.String));
if strlength(data.compositionOutputVar) == 0
    data.compositionOutputVar = "conc";
end
data.compositionInfoVar = strtrim(string(ctrls.compInfoVarEdit.String));
if strlength(data.compositionInfoVar) == 0
    data.compositionInfoVar = data.compositionOutputVar + "Info";
end

data.compositionPlotBackground = logical(ctrls.compPlotBgCb.Value);
data.compositionPlotFits = logical(ctrls.compPlotFitsCb.Value);
data.compositionExtraOptions = strtrim(string(ctrls.compExtraOptsEdit.String));
end

function [posRaw, sourceLabel, ok, msg] = resolveCompositionPos(data)
posRaw = table();
sourceLabel = "";
ok = false;
msg = "";

hasSourcePos = isfield(data.sourceInfo, 'pos') && istable(data.sourceInfo.pos) && ~isempty(data.sourceInfo.pos);
if logical(data.compositionUseSourcePos) && hasSourcePos
    posRaw = data.sourceInfo.pos;
    sourceLabel = "source pos";
    ok = true;
else
    varName = char(strtrim(data.compositionPosVarName));
    if isempty(varName)
        msg = 'No pos variable name configured.';
        return;
    end
    if ~isvarname(varName)
        msg = sprintf('Invalid MATLAB variable name: "%s".', varName);
        return;
    end
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName)) ~= 1
        msg = sprintf('Workspace variable "%s" not found.', varName);
        return;
    end
    posCandidate = evalin('base', varName);
    if ~istable(posCandidate) || ~ismember('mc', posCandidate.Properties.VariableNames)
        msg = sprintf('Variable "%s" must be a pos table with column mc.', varName);
        return;
    end
    posRaw = posCandidate;
    sourceLabel = "workspace:" + string(varName);
    ok = true;
end

if ~ismember('mc', posRaw.Properties.VariableNames)
    ok = false;
    msg = 'Input pos table must contain column mc.';
end
end

function list = parseExcludeList(textIn)
txt = strtrim(char(string(textIn)));
if isempty(txt)
    list = {};
    return;
end
tokens = regexp(txt, '[^,\s;]+', 'match');
if isempty(tokens)
    list = {};
    return;
end
list = unique(tokens, 'stable');
end

function [opts, ok, msg] = parseExtraNameValueList(textIn)
opts = {};
ok = true;
msg = "";
txt = strtrim(char(string(textIn)));
if isempty(txt)
    return;
end
try
    opts = eval(['{' txt '}']); %#ok<EVLDIR>
catch ME
    ok = false;
    msg = "Extra NV parse failed: " + string(ME.message);
    opts = {};
    return;
end
if ~iscell(opts) || mod(numel(opts), 2) ~= 0
    ok = false;
    msg = 'Extra NV must evaluate to an even-length name-value list.';
    opts = {};
end
end

function s = readPopupString(hPopup, fallback)
s = string(fallback);
if ~isgraphics(hPopup)
    return;
end
items = string(hPopup.String);
if isempty(items)
    return;
end
idx = hPopup.Value;
if ~isscalar(idx) || ~isfinite(idx) || idx < 1 || idx > numel(items)
    idx = 1;
end
s = string(items(idx));
end

function setPopupString(hPopup, valueIn)
if ~isgraphics(hPopup)
    return;
end
items = string(hPopup.String);
if isempty(items)
    return;
end
value = string(valueIn);
idx = find(lower(items) == lower(value), 1, 'first');
if isempty(idx)
    idx = 1;
end
hPopup.Value = idx;
end

function setEnable(h, tf)
if isgraphics(h)
    h.Enable = onOff(tf);
end
end

function out = onOff(tf)
if logical(tf)
    out = 'on';
else
    out = 'off';
end
end

function varName = sanitizeWorkspaceVarName(nameIn, fallback)
candidate = strtrim(char(string(nameIn)));
if isempty(candidate)
    candidate = char(fallback);
end
varName = string(matlab.lang.makeValidName(candidate));
if strlength(varName) == 0
    varName = string(matlab.lang.makeValidName(char(fallback)));
end
end

function onImportIons(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end

[sourceSpec, tempFig, sourceLabel] = resolveIonImportSource(controlFig, data);
if isempty(sourceSpec) || ~isgraphics(sourceSpec)
    return;
end

try
    ionTable = ionsExtractFromMassSpec(sourceSpec);
catch ME
    warning('massSpecWidget:importIonsExtractFailed', ...
        'Could not extract ions from source spectrum: %s', ME.message);
    if ~isempty(tempFig) && isgraphics(tempFig)
        delete(tempFig);
    end
    setStatus(controlFig, 'Could not extract ions from selected source.');
    return;
end

if isempty(ionTable) || ~istable(ionTable) || height(ionTable) == 0
    if ~isempty(tempFig) && isgraphics(tempFig)
        delete(tempFig);
    end
    setStatus(controlFig, 'No ions found in selected source.');
    return;
end

added = 0;
for i = 1:height(ionTable)
    try
        ionName = string(ionTable.ionName(i));
        cs = round(double(ionTable.chargeState(i)));
        if strlength(ionName) == 0 || ~isfinite(cs) || cs <= 0
            continue;
        end
        data.colorScheme = ensureIonColor(data.colorScheme, ionName);
        removeMatchingIon(data.spec, ionName, cs);
        ionAdd(data.spec, char(ionName), cs, data.isotopeTable, data.colorScheme, ...
            0, 0.01, 'most abundant', 0.1);
        added = added + 1;
    catch ME
        warning('massSpecWidget:importIonsAddFailed', ...
            'Could not import ion row %d: %s', i, ME.message);
    end
end

if ~isempty(tempFig) && isgraphics(tempFig)
    delete(tempFig);
end

setappdata(controlFig, 'massSpecWidget', data);
refreshUi(controlFig);
setStatus(controlFig, sprintf('Imported %d ions from %s.', added, sourceLabel));
end

function [sourceSpec, tempFig, sourceLabel] = resolveIonImportSource(controlFig, data)
sourceSpec = [];
tempFig = [];
sourceLabel = "source";

[choiceIdx, ok] = listdlg( ...
    'PromptString', 'Import ions from:', ...
    'ListString', {'Open .fig file', 'Open mass-spectrum figure'}, ...
    'SelectionMode', 'single', ...
    'InitialValue', 1, ...
    'ListSize', [260 120], ...
    'Name', 'Import Ions');
if ~ok || isempty(choiceIdx)
    return;
end

switch choiceIdx
    case 1
        [file, path] = uigetfile({'*.fig', 'MATLAB Figure (*.fig)'}, 'Select source mass spectrum figure');
        if isequal(file, 0)
            return;
        end
        sourceFile = fullfile(path, file);
        tempFig = openfig(sourceFile, 'new', 'invisible');
        axList = findobj(tempFig, 'Type', 'axes');
        for i = 1:numel(axList)
            sourceSpec = findMassSpectrumInAxes(axList(i));
            if ~isempty(sourceSpec)
                break;
            end
        end
        sourceLabel = string(file);

    case 2
        figs = findall(groot, 'Type', 'figure');
        if isempty(figs)
            setStatus(controlFig, 'No open figure found.');
            return;
        end
        names = strings(0, 1);
        validFig = gobjects(0);
        validSpec = gobjects(0);
        for i = 1:numel(figs)
            f = figs(i);
            if f == controlFig
                continue;
            end
            axesList = findobj(f, 'Type', 'axes');
            sp = [];
            for j = 1:numel(axesList)
                sp = findMassSpectrumInAxes(axesList(j));
                if ~isempty(sp)
                    break;
                end
            end
            if isempty(sp)
                continue;
            end
            validFig(end+1,1) = f; %#ok<AGROW>
            validSpec(end+1,1) = sp; %#ok<AGROW>
            figName = string(f.Name);
            if strlength(figName) == 0
                figName = "Figure " + string(double(f.Number));
            end
            names(end+1,1) = figName; %#ok<AGROW>
        end
        if isempty(validFig)
            setStatus(controlFig, 'No open figure with mass spectrum found.');
            return;
        end
        [idx, ok2] = listdlg( ...
            'PromptString', 'Select source figure:', ...
            'ListString', cellstr(names), ...
            'SelectionMode', 'single', ...
            'InitialValue', 1, ...
            'ListSize', [360 240], ...
            'Name', 'Import Ions');
        if ~ok2 || isempty(idx)
            return;
        end
        sourceSpec = validSpec(idx);
        sourceLabel = names(idx);
end

if isempty(sourceSpec) || ~isgraphics(sourceSpec)
    setStatus(controlFig, 'Could not resolve source mass spectrum.');
end
end

function onKeyPress(src, evd)
try
    controlFig = src;
    if ~isgraphics(controlFig)
        return;
    end

    activeObj = gco(controlFig);
    if shouldIgnoreHotkey(activeObj, evd)
        return;
    end

    isCtrl = any(strcmp(evd.Modifier, 'control')) || any(strcmp(evd.Modifier, 'command'));
    key = lower(string(evd.Key));

    if key == "slash" || key == "question" || key == "f1"
        showKeyboardHelp();
        return;
    end

    if isCtrl
        switch key
            case "s"
                onSaveProfileWorkspace(controlFig);
            case "l"
                onLoadProfileWorkspace(controlFig);
            case "e"
                onExportProfileFile(controlFig);
            case "i"
                onImportProfileFile(controlFig);
        end
        return;
    end

    switch key
        case "r"
            refreshUi(controlFig);
        case "i"
            onAddIon(controlFig);
        case "f"
            onFindIon(controlFig);
        case "a"
            onAddRange(controlFig);
        case "g"
            onAddAllRanges(controlFig);
        case "x"
            onExtendRange(controlFig);
        case "delete"
            onDeleteSelected(controlFig);
        case "backspace"
            onDeleteSelected(controlFig);
        case "l"
            onShowLegend(controlFig);
        case "o"
            onReorder(controlFig);
        case "z"
            onZoomFull(controlFig);
        case "v"
            onZoomRanged(controlFig);
        case "e"
            onExtractElements(controlFig);
        case "b"
            onBuildCandidateIons(controlFig);
        case "c"
            onComputeConcentration(controlFig);
        case "p"
            onApplyProfileToDataset(controlFig);
    end
catch ME
    warning('massSpecWidget:keyPressError', 'Key press error: %s', ME.message);
end
end

function onDeleteSelected(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end

if ~isempty(data.selectedRangeRow)
    onRemoveRange(controlFig);
    return;
end
if ~isempty(data.selectedIonRow)
    onRemoveIon(controlFig);
    return;
end
setStatus(controlFig, 'Select an ion or range row, then press Delete.');
end

function tf = shouldIgnoreHotkey(activeObj, evd)
tf = false;
if isempty(activeObj) || ~isgraphics(activeObj)
    return;
end

% Keep text entry/edit interactions untouched unless control/command is held.
isCtrl = any(strcmp(evd.Modifier, 'control')) || any(strcmp(evd.Modifier, 'command'));
if isCtrl
    return;
end

objType = lower(string(get(activeObj, 'Type')));
if objType == "uicontrol"
    try
        style = lower(string(get(activeObj, 'Style')));
    catch
        style = "";
    end
    if any(style == ["edit", "listbox", "popupmenu"])
        tf = true;
        return;
    end
end

if contains(class(activeObj), 'uitable', 'IgnoreCase', true) || objType == "uitable"
    tf = true;
end
end

function showKeyboardHelp()
helpText = sprintf([ ...
    'Keyboard Shortcuts:\n\n' ...
    'R            Refresh tables\n' ...
    'I            Add ion from Ion field\n' ...
    'F            Launch ionFind\n' ...
    'A            Add one range\n' ...
    'G            Run Add All ranges\n' ...
    'X            Extend range to neighbor\n' ...
    'Delete       Remove selected range (or ion)\n' ...
    'L            Show legend\n' ...
    'O            Reorder overlays\n' ...
    'Z            Zoom full\n' ...
    'V            Zoom ions/ranges\n' ...
    'E            Extract elements from text\n' ...
    'B            Build potential ion list from elements\n' ...
    'C            Calculate concentration\n' ...
    'P            Apply current profile to selected data/fig\n\n' ...
    'Ctrl/Cmd+S   Save profile to workspace\n' ...
    'Ctrl/Cmd+L   Load profile from workspace\n' ...
    'Ctrl/Cmd+E   Export profile to file\n' ...
    'Ctrl/Cmd+I   Import profile from file\n' ...
    'F1 or ?      Show this help\n']);
msgbox(helpText, 'Mass Spectrum Widget Shortcuts', 'help');
end

function configureControlWindowLayout(controlFig)
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(ctrls) || ~isgraphics(controlFig)
    return;
end

layout = struct();
layout.margins.outer = 20;
layout.margins.colGap = 12;
layout.margins.rowGap = 10;
layout.minFigW = 900;
layout.minFigH = 760;

set(controlFig, 'Units', 'pixels');

panels = [ctrls.sourcePanel, ctrls.ionPanel, ctrls.rangePanel, ctrls.elementPanel, ctrls.profilePanel, ctrls.candidatePanel];
for i = 1:numel(panels)
    if isgraphics(panels(i))
        set(panels(i), 'Units', 'pixels');
    end
end

if isgraphics(ctrls.sourcePanel)
    p = getpixelposition(ctrls.sourcePanel);
    layout.sourceHeight = p(4);
end
if isgraphics(ctrls.profilePanel)
    p = getpixelposition(ctrls.profilePanel);
    layout.profileHeight = p(4);
end
if isgraphics(ctrls.candidatePanel)
    p = getpixelposition(ctrls.candidatePanel);
    layout.rightWidth = p(3);
end
if ~isfield(layout, 'rightWidth') || ~isfinite(layout.rightWidth)
    layout.rightWidth = 240;
end
if ~isfield(layout, 'sourceHeight')
    layout.sourceHeight = 52;
end
if ~isfield(layout, 'profileHeight')
    layout.profileHeight = 92;
end

% Preserve relative vertical emphasis of middle left panels.
ionH = getpixelposition(ctrls.ionPanel);
rangeH = getpixelposition(ctrls.rangePanel);
elemH = getpixelposition(ctrls.elementPanel);
pref = [ionH(4), rangeH(4), elemH(4)];
if any(pref <= 0)
    pref = [300, 220, 160];
end
layout.leftMiddlePref = pref;
layout.leftMiddleMin = [180, 140, 120];

layout.panels = struct();
layout.panels.source = capturePanelLayout(ctrls.sourcePanel, gobjects(0));
layout.panels.ion = capturePanelLayout(ctrls.ionPanel, ctrls.ionTable);
layout.panels.range = capturePanelLayout(ctrls.rangePanel, ctrls.rangeTable);
layout.panels.element = capturePanelLayout(ctrls.elementPanel, ctrls.elementsList);
layout.panels.profile = capturePanelLayout(ctrls.profilePanel, gobjects(0));
layout.panels.candidate = capturePanelLayout(ctrls.candidatePanel, ctrls.candidateTable);

setappdata(controlFig, 'massSpecWidgetLayout', layout);
set(controlFig, 'SizeChangedFcn', @(src, ~) onControlWindowResize(src));
onControlWindowResize(controlFig);
end

function onControlWindowResize(controlFig)
layout = getappdata(controlFig, 'massSpecWidgetLayout');
if isempty(layout) || ~isgraphics(controlFig)
    return;
end

figPos = getpixelposition(controlFig);
figW = figPos(3);
figH = figPos(4);

if figW < layout.minFigW || figH < layout.minFigH
    newPos = figPos;
    newPos(3) = max(figW, layout.minFigW);
    newPos(4) = max(figH, layout.minFigH);
    if any(newPos(3:4) ~= figPos(3:4))
        setpixelposition(controlFig, newPos);
        return;
    end
end

outer = layout.margins.outer;
colGap = layout.margins.colGap;
rowGap = layout.margins.rowGap;
rightW = min(layout.rightWidth, max(240, floor(0.34 * figW)));
leftW = max(480, figW - 2*outer - colGap - rightW);
rightX = outer + leftW + colGap;
leftX = outer;

topY = figH - outer;
sourceH = layout.sourceHeight;
profileH = layout.profileHeight;
profileY = outer;
sourceY = topY - sourceH;

middleTop = sourceY - rowGap;
middleBottom = profileY + profileH + rowGap;
middleH = max(120, middleTop - middleBottom);

leftHeights = distributeHeights(middleH - 2*rowGap, layout.leftMiddlePref, layout.leftMiddleMin);
ionH = leftHeights(1);
rangeH = leftHeights(2);
elemH = leftHeights(3);

ionY = middleBottom + rangeH + elemH + 2*rowGap;
rangeY = middleBottom + elemH + rowGap;
elemY = middleBottom;

setpanel(layout.panels.source.panel, [leftX, sourceY, leftW, sourceH]);
setpanel(layout.panels.profile.panel, [outer, profileY, max(300, figW - 2*outer), profileH]);
setpanel(layout.panels.ion.panel, [leftX, ionY, leftW, ionH]);
setpanel(layout.panels.range.panel, [leftX, rangeY, leftW, rangeH]);
setpanel(layout.panels.element.panel, [leftX, elemY, leftW, elemH]);
setpanel(layout.panels.candidate.panel, [rightX, middleBottom, rightW, middleH + sourceH + rowGap]);

applyPanelLayout(layout.panels.source);
applyPanelLayout(layout.panels.profile);
applyPanelLayout(layout.panels.ion);
applyPanelLayout(layout.panels.range);
applyPanelLayout(layout.panels.element);
applyPanelLayout(layout.panels.candidate);
end

function cfg = capturePanelLayout(panelHandle, resizableHandle)
cfg = struct();
cfg.panel = panelHandle;
cfg.fixedHandles = gobjects(0);
cfg.fixedLeft = [];
cfg.fixedTop = [];
cfg.fixedWidth = [];
cfg.fixedHeight = [];
cfg.resizableHandles = gobjects(0);
cfg.resizableLeft = [];
cfg.resizableRight = [];
cfg.resizableTop = [];
cfg.resizableBottom = [];
cfg.resizableMinH = [];
cfg.resizableMinW = [];

if ~isgraphics(panelHandle)
    return;
end
set(panelHandle, 'Units', 'pixels');
pp = getpixelposition(panelHandle);

allChildren = panelHandle.Children;
resizable = resizableHandle(isgraphics(resizableHandle));
for i = 1:numel(allChildren)
    h = allChildren(i);
    if ~isgraphics(h)
        continue;
    end
    set(h, 'Units', 'pixels');
    hp = getpixelposition(h);
    if any(h == resizable)
        cfg.resizableHandles(end+1,1) = h; %#ok<AGROW>
        cfg.resizableLeft(end+1,1) = hp(1); %#ok<AGROW>
        cfg.resizableRight(end+1,1) = pp(3) - (hp(1) + hp(3)); %#ok<AGROW>
        cfg.resizableTop(end+1,1) = pp(4) - (hp(2) + hp(4)); %#ok<AGROW>
        cfg.resizableBottom(end+1,1) = hp(2); %#ok<AGROW>
        cfg.resizableMinH(end+1,1) = max(80, min(140, hp(4))); %#ok<AGROW>
        cfg.resizableMinW(end+1,1) = max(60, min(200, hp(3))); %#ok<AGROW>
    else
        cfg.fixedHandles(end+1,1) = h; %#ok<AGROW>
        cfg.fixedLeft(end+1,1) = hp(1); %#ok<AGROW>
        cfg.fixedTop(end+1,1) = pp(4) - (hp(2) + hp(4)); %#ok<AGROW>
        cfg.fixedWidth(end+1,1) = hp(3); %#ok<AGROW>
        cfg.fixedHeight(end+1,1) = hp(4); %#ok<AGROW>
    end
end
end

function applyPanelLayout(cfg)
if ~isgraphics(cfg.panel)
    return;
end
pp = getpixelposition(cfg.panel);

for i = 1:numel(cfg.fixedHandles)
    h = cfg.fixedHandles(i);
    if ~isgraphics(h)
        continue;
    end
    x = cfg.fixedLeft(i);
    y = pp(4) - cfg.fixedTop(i) - cfg.fixedHeight(i);
    setpixelposition(h, [x, y, cfg.fixedWidth(i), cfg.fixedHeight(i)]);
end

for i = 1:numel(cfg.resizableHandles)
    h = cfg.resizableHandles(i);
    if ~isgraphics(h)
        continue;
    end
    x = cfg.resizableLeft(i);
    y = cfg.resizableBottom(i);
    w = max(cfg.resizableMinW(i), pp(3) - cfg.resizableLeft(i) - cfg.resizableRight(i));
    hgt = max(cfg.resizableMinH(i), pp(4) - cfg.resizableTop(i) - cfg.resizableBottom(i));
    setpixelposition(h, [x, y, w, hgt]);
end
end

function setpanel(h, pos)
if isgraphics(h)
    setpixelposition(h, pos);
end
end

function h = distributeHeights(totalH, pref, mins)
pref = double(pref(:)');
mins = double(mins(:)');
if totalH <= 0
    h = mins;
    return;
end
sumMin = sum(mins);
if totalH <= sumMin
    h = mins .* (totalH / sumMin);
    return;
end
extra = totalH - sumMin;
pw = pref / max(sum(pref), eps);
h = mins + extra .* pw;
end

function dataOut = candidateTableData(data)
if isempty(data.candidateIons) || ~istable(data.candidateIons) || height(data.candidateIons) == 0
    dataOut = cell(0,2);
    return;
end
n = height(data.candidateIons);
dataOut = cell(n, 2);
for i = 1:n
    dataOut{i,1} = char(string(data.candidateIons.label(i)));
    if ismember('mcAll', data.candidateIons.Properties.VariableNames)
        dataOut{i,2} = char(formatMcVector(data.candidateIons.mcAll{i}));
    else
        dataOut{i,2} = char(formatMcVector(double(data.candidateIons.mc(i))));
    end
end
end

function txt = formatMcVector(vals)
txt = "[]";
if ischar(vals) || isstring(vals)
    txt = string(vals);
    return;
end
if isempty(vals)
    return;
end
v = double(vals(:)');
v = v(isfinite(v));
if isempty(v)
    return;
end
v = sort(unique(v), 'ascend');
parts = arrayfun(@(x) strtrim(num2str(x, '%.6g')), v, 'UniformOutput', false);
txt = "[" + strjoin(string(parts), " ") + "]";
end

function out = summarizeCandidateIons(ionList)
if isempty(ionList) || ~istable(ionList) || ...
        ~all(ismember({'ion','ionIsotopic','mc'}, ionList.Properties.VariableNames))
    out = table();
    return;
end

ionBase = string(ionList.ion);
ionIso = string(ionList.ionIsotopic);
charge = zeros(height(ionList),1);
for i = 1:height(ionList)
    charge(i) = countTrailingPlus(ionIso(i));
    if charge(i) <= 0
        charge(i) = 1;
    end
end

[G, baseNames, cs] = findgroups(ionBase, charge);
mcRef = nan(numel(baseNames), 1);
mcAll = cell(numel(baseNames), 1);
label = strings(numel(baseNames),1);
for i = 1:numel(baseNames)
    vals = sort(unique(double(ionList.mc(G == i))), 'ascend');
    vals = vals(isfinite(vals));
    mcAll{i} = vals(:)';
    if ~isempty(vals)
        mcRef(i) = vals(1);
    end
    label(i) = baseNames(i) + repmat('+', 1, cs(i));
end

out = table(string(baseNames), cs, mcRef, mcAll, label, ...
    'VariableNames', {'ionName', 'chargeState', 'mc', 'mcAll', 'label'});
out = sortrows(out, {'mc', 'ionName', 'chargeState'});
end

function subList = buildCandidateSubIonList(ionName, chargeState, isotopeTable)
[isoCombos, ~, weights] = ionsCreateIsotopeList(char(ionName), isotopeTable);
if isempty(weights)
    subList = table();
    return;
end

ionCol = repmat(categorical(string(ionName)), numel(weights), 1);
isoCol = strings(numel(weights), 1);
for i = 1:numel(weights)
    isoCol(i) = string(ionConvertName(isoCombos{i})) + repmat('+',1,chargeState);
end
mcCol = double(weights(:)) ./ double(chargeState);
subList = table(ionCol, categorical(isoCol), mcCol, ...
    'VariableNames', {'ion', 'ionIsotopic', 'mc'});
end

function clearPreviewIon(controlFig)
data = getappdata(controlFig, 'massSpecWidget');
if isempty(data)
    return;
end
if isfield(data, 'previewIonHandle') && ~isempty(data.previewIonHandle)
    validHandles = data.previewIonHandle(isgraphics(data.previewIonHandle));
    if ~isempty(validHandles)
        delete(validHandles);
    end
end
data.previewIonHandle = gobjects(0);
data.previewCandidate = struct('ionName', "", 'chargeState', NaN);
setappdata(controlFig, 'massSpecWidget', data);
end

function removeMatchingIon(spec, ionName, chargeState)
ax = ancestor(spec, 'axes');
plots = ax.Children;
for i = 1:numel(plots)
    h = plots(i);
    if ~isgraphics(h)
        continue;
    end
    try
        if ~isfield(h.UserData, 'plotType') || h.UserData.plotType ~= "ion"
            continue;
        end
    catch
        continue;
    end
    try
        base = string(h.UserData.ion{1}.element);
        cs = double(h.UserData.chargeState(1));
    catch
        continue;
    end
    if base == ionName && cs == chargeState
        delete(h);
    end
end
end

function tf = ionExistsInSpec(spec, ionName, chargeState, isTracer, ignoreHandle)
tf = false;
if nargin < 4 || isempty(isTracer)
    isTracer = false;
end
if nargin < 5
    ignoreHandle = [];
end
if isempty(spec) || ~isgraphics(spec)
    return;
end
rows = collectIonRows(spec);
if isempty(rows)
    return;
end
targetBase = canonicalIonBaseName(ionName);
if strlength(targetBase) == 0
    return;
end
names = string({rows.baseName});
charges = [rows.chargeState];
tracerFlags = logical([rows.isTracer]);

mask = names == targetBase & charges == double(chargeState) & tracerFlags == logical(isTracer);
if ~isempty(ignoreHandle) && isgraphics(ignoreHandle)
    handles = [rows.handle];
    mask = mask & handles ~= ignoreHandle;
end
tf = any(mask);
end

function baseName = canonicalIonBaseName(ionName)
baseName = "";
try
    [ionParsed, ~, ~] = ionConvertName(char(string(ionName)));
    if istable(ionParsed) && ~isempty(ionParsed) && ismember('element', ionParsed.Properties.VariableNames)
        baseName = string(ionConvertName(ionParsed.element));
    end
catch
end
if strlength(baseName) > 0
    baseName = strtrim(baseName);
    return;
end
txt = char(string(ionName));
txt = regexprep(txt, '\s*tracer\s*$', '', 'ignorecase');
txt = regexprep(txt, '[+-]+$', '');
baseName = strtrim(string(txt));
end

function [baseName, chargeState, isTracer] = ionIdentityFromHandle(hIon)
baseName = "";
chargeState = NaN;
isTracer = false;
if ~isgraphics(hIon)
    return;
end
try
    ud = hIon.UserData;
catch
    ud = struct();
end
if isstruct(ud)
    try
        if isfield(ud, 'chargeState') && ~isempty(ud.chargeState)
            chargeState = double(ud.chargeState(1));
        end
    catch
    end
    try
        if isfield(ud, 'isTracer')
            isTracer = logical(ud.isTracer(1));
        end
    catch
    end
    try
        if isfield(ud, 'ion') && ~isempty(ud.ion)
            ionData = ud.ion;
            if iscell(ionData)
                ionData = ionData{1};
            end
            if istable(ionData) && ismember('element', ionData.Properties.VariableNames)
                baseName = string(ionConvertName(ionData.element));
            end
        end
    catch
    end
end
if strlength(baseName) == 0
    try
        baseName = canonicalIonBaseName(string(hIon.DisplayName));
    catch
        baseName = "";
    end
end
end

function n = countTrailingPlus(textIn)
txt = char(string(textIn));
n = 0;
for i = length(txt):-1:1
    if txt(i) == '+'
        n = n + 1;
    else
        break;
    end
end
end

function vals = parsePositiveVector(textIn, defaultVal)
if nargin < 2
    defaultVal = 1;
end
if isnumeric(textIn)
    vals = round(double(textIn(:)'));
else
    vals = str2num(strtrim(char(string(textIn)))); %#ok<ST2NM>
end
if isempty(vals)
    vals = defaultVal;
end
vals = round(vals(:)');
vals = vals(isfinite(vals) & vals > 0);
if isempty(vals)
    vals = defaultVal;
end
vals = unique(vals, 'stable');
end

function textOut = sourceSummaryText(data)
srcType = string(data.sourceInfo.type);
if strlength(srcType) == 0
    srcType = "unknown";
end
textOut = "Source: " + srcType;
if isfield(data.sourceInfo, 'path') && strlength(string(data.sourceInfo.path)) > 0
    [~, baseName, ext] = fileparts(char(string(data.sourceInfo.path)));
    textOut = textOut + " | " + string(baseName + ext);
end
end

function setStatus(controlFig, msg)
ctrls = getappdata(controlFig, 'massSpecWidgetControls');
if isempty(ctrls) || ~isfield(ctrls, 'statusText') || ~isgraphics(ctrls.statusText)
    return;
end
ctrls.statusText.String = char(msg);
drawnow limitrate nocallbacks;
end

function textOut = readMultilineEdit(hEdit)
raw = hEdit.String;
if ischar(raw)
    textOut = string(raw);
elseif iscell(raw)
    textOut = strjoin(string(raw), " ");
else
    textOut = string(raw);
end
end

function values = parseChargeStates(textIn, defaultValues)
if nargin < 2 || isempty(defaultValues)
    defaultValues = 1;
end
defaults = round(double(defaultValues(:)'));
defaults = defaults(isfinite(defaults) & defaults > 0);
if isempty(defaults)
    defaults = 1;
end

if isnumeric(textIn)
    values = round(double(textIn(:)'));
else
    str = strtrim(char(string(textIn)));
    if isempty(str)
        values = defaults;
        return;
    end

    % Numeric expressions remain supported (e.g. "1 2", "1:3", "[1 2 3]").
    values = str2num(str); %#ok<ST2NM>
    if isempty(values)
        values = [];
        tokens = regexp(str, '[^,\s;]+', 'match');
        for i = 1:numel(tokens)
            tok = strtrim(tokens{i});
            if isempty(tok)
                continue;
            end

            % "+" and "++" style charge state syntax.
            if all(tok == '+')
                values(end+1) = numel(tok); %#ok<AGROW>
                continue;
            end

            nPlus = countTrailingPlus(tok);
            if nPlus > 0
                core = strtrim(tok(1:end-nPlus));
                if isempty(core) || strcmp(core, '+')
                    values(end+1) = nPlus; %#ok<AGROW>
                    continue;
                end
                coreVal = str2double(core);
                if isfinite(coreVal) && coreVal > 0
                    values(end+1) = round(coreVal); %#ok<AGROW>
                    continue;
                end
            end

            val = str2double(tok);
            if isfinite(val) && val > 0
                values(end+1) = round(val); %#ok<AGROW>
            end
        end
    end
end
values = round(values(:)');
values = values(isfinite(values) & values > 0);
if isempty(values)
    values = defaults;
end
values = unique(values, 'stable');
end

function elements = extractElementsRobust(textIn)
if strlength(textIn) == 0
    elements = strings(0, 1);
    return;
end

allElements = periodicElements();
elements = strings(0, 1);

try
    base = elementsExtractFromText(textIn);
    elements = [elements; string(base(:))]; %#ok<AGROW>
catch
end

txt = char(textIn);
tok = regexp(txt, '(?<![A-Za-z])([A-Z][a-z]?)(?![A-Za-z])', 'tokens');
if ~isempty(tok)
    tok = string([tok{:}]);
    elements = [elements; tok(:)]; %#ok<AGROW>
end

tokCaps = regexp(txt, '(?<![A-Za-z])([A-Z]{1,2})(?![A-Za-z])', 'tokens');
if ~isempty(tokCaps)
    raw = string([tokCaps{:}]);
    mapped = strings(size(raw));
    for i = 1:numel(raw)
        t = char(raw(i));
        if numel(t) == 1
            mapped(i) = string(t);
        else
            mapped(i) = string([t(1), lower(t(2))]);
        end
    end
    elements = [elements; mapped(:)]; %#ok<AGROW>
end

elements = strtrim(elements);
elements = elements(elements ~= "");
elements = elements(ismember(elements, allElements));
elements = unique(elements, 'stable');
end

function colorScheme = resolveColorScheme(colorSchemeOpt)
if ~isempty(colorSchemeOpt)
    colorScheme = normalizeColorScheme(colorSchemeOpt);
    return;
end

if evalin('base', 'exist(''colorScheme'', ''var'')')
    cs = evalin('base', 'colorScheme');
    if istable(cs)
        colorScheme = normalizeColorScheme(cs);
        return;
    end
end

if isfile('colorScheme.mat')
    s = load('colorScheme.mat');
    if isfield(s, 'colorScheme') && istable(s.colorScheme)
        colorScheme = normalizeColorScheme(s.colorScheme);
        return;
    end
end

colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
    'VariableNames', {'ion', 'color'});
end

function isotopeTable = resolveIsotopeTable(isotopeTableOpt)
if ~isempty(isotopeTableOpt)
    isotopeTable = isotopeTableOpt;
    return;
end
s = load('isotopeTable_naturalAbundances.mat', 'isotopeTable');
isotopeTable = s.isotopeTable;
end

function ionList = resolveIonList(ionListOpt)
if ~isempty(ionListOpt)
    ionList = ionListOpt;
    return;
end

ionList = [];
if evalin('base', 'exist(''ionList'', ''var'')')
    ionList = evalin('base', 'ionList');
    return;
end

if isfile('IonenlisteAPT.xlsx')
    try
        ionList = readtable('IonenlisteAPT.xlsx', 'VariableNamingRule', 'preserve');
    catch
        ionList = [];
    end
end
end

function colorScheme = normalizeColorScheme(colorSchemeIn)
if isempty(colorSchemeIn)
    colorScheme = table(categorical.empty(0,1), zeros(0,3), ...
        'VariableNames', {'ion', 'color'});
    return;
end
if ~istable(colorSchemeIn) || ...
        ~all(ismember({'ion', 'color'}, colorSchemeIn.Properties.VariableNames))
    error('massSpecWidget:invalidColorScheme', ...
        'colorScheme must be a table with columns ion and color.');
end
colorScheme = colorSchemeIn(:, {'ion', 'color'});
colorScheme.ion = categorical(string(colorScheme.ion));
colorScheme.color = double(colorScheme.color);
if size(colorScheme.color, 2) ~= 3
    error('massSpecWidget:invalidColorSchemeColor', ...
        'colorScheme.color must be Nx3.');
end
end

function colorScheme = ensureBackgroundColor(colorScheme)
if any(string(colorScheme.ion) == "background")
    return;
end
colorScheme(end+1, :) = {categorical("background"), [0.7 0.7 0.7]};
end

function colorScheme = ensureIonColor(colorScheme, ionName)
baseIon = normalizeIonForColorScheme(ionName);
if any(string(colorScheme.ion) == baseIon)
    return;
end
try
    colorScheme = colorSchemeIonAdd(colorScheme, char(baseIon), 'create');
catch
    colorScheme(end+1, :) = {categorical(baseIon), rand(1, 3)};
end
end

function ionOut = normalizeIonForColorScheme(ionIn)
ionOut = strtrim(string(ionIn));
try
    ionOut = string(ionConvertMode(categorical(ionOut), 'atomic'));
catch
end
if strlength(ionOut) == 0
    ionOut = "unknown";
end
end

function varName = selectWorkspaceProfileVariable(defaultName)
varName = "";
ws = evalin('base', 'whos');
if isempty(ws)
    errordlg('Workspace is empty. No profile variable available.', ...
        'Load Profile', 'modal');
    return;
end

isStruct = strcmp({ws.class}, 'struct');
names = string({ws(isStruct).name});
if isempty(names)
    errordlg('No struct variables found in workspace.', ...
        'Load Profile', 'modal');
    return;
end

[~, order] = sort(lower(names));
names = names(order);
idxDefault = find(names == defaultName, 1, 'first');
if isempty(idxDefault)
    idxDefault = 1;
end

[idx, ok] = listdlg( ...
    'PromptString', 'Select profile variable:', ...
    'SelectionMode', 'single', ...
    'ListString', cellstr(names), ...
    'InitialValue', idxDefault, ...
    'ListSize', [360 280], ...
    'Name', 'Load Mass Spectrum Profile');
if ~ok || isempty(idx)
    return;
end
varName = names(idx);
end

function elems = periodicElements()
elems = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar", ...
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr", ...
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe", ...
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf", ...
    "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th", ...
    "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs", ...
    "Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"]';
end

function safeDelete(h)
if ~isempty(h) && isgraphics(h)
    delete(h);
end
end

function [ok, msg, data] = replaceIonFromTableRow(data, tbl, row)
ok = false;
msg = "";
ionName = strtrim(string(tbl{row, 1}));
if strlength(ionName) == 0
    msg = "Ion name cannot be empty.";
    return;
end
isTracer = parseLogicalScalar(tbl{row, 2}, data.ionRows(row).isTracer);
oldHandle = data.ionRows(row).handle;

ionNameToParse = ionName;
if isTracer && ~contains(lower(ionNameToParse), 'tracer')
    ionNameToParse = ionNameToParse + " tracer";
end

try
    [ionParsed, chargeState, isTracerParsed] = ionConvertName(char(ionNameToParse));
catch ME
    msg = string(ME.message);
    return;
end
if ~isfinite(chargeState) || chargeState <= 0
    msg = "Ion name must include charge state (e.g. Fe+, Fe++).";
    return;
end
isTracer = logical(isTracerParsed);
ionBase = string(ionConvertName(ionParsed.element));
ionNameToAdd = ionBase;
if isTracer && ~contains(lower(ionNameToAdd), 'tracer')
    ionNameToAdd = ionNameToAdd + " tracer";
end

if ionExistsInSpec(data.spec, ionNameToAdd, chargeState, isTracer, oldHandle)
    ok = true;
    msg = "Ion already enlisted; edit ignored.";
    return;
end

beforeRows = collectIonRows(data.spec);
beforeHandles = [beforeRows.handle];

try
    data.colorScheme = ensureIonColor(data.colorScheme, ionNameToAdd);
    ionAdd(data.spec, char(ionNameToAdd), chargeState, data.isotopeTable, data.colorScheme, ...
        0, 0.01, 'most abundant', 0.1);
catch ME
    msg = string(ME.message);
    return;
end

afterRows = collectIonRows(data.spec);
afterHandles = [afterRows.handle];
newHandle = findAddedHandle(beforeHandles, afterHandles);
if isempty(newHandle)
    msg = "Could not resolve newly added ion handle.";
    return;
end

oldHandle = data.ionRows(row).handle;
if isgraphics(oldHandle)
    delete(oldHandle);
end

try
    if isfinite(chargeState)
        ud = newHandle.UserData;
        if ~isstruct(ud)
            ud = struct();
        end
        ud.chargeState = chargeState;
        ud.isTracer = logical(isTracer);
        newHandle.UserData = ud;
    end
catch
end

ok = true;
updatedName = string(ionConvertName(ionParsed.element, chargeState, 'plain', isTracer));
msg = sprintf('Ion updated: %s.', updatedName);
end

function [ok, msg, data] = replaceRangeFromTableRow(data, tbl, row)
ok = false;
msg = "";

rangeName = strtrim(string(tbl{row, 1}));
if strlength(rangeName) == 0
    rangeName = "background";
end

mcbegin = parseFiniteScalar(tbl{row, 2}, data.rangeRows(row).mcbegin);
mcend = parseFiniteScalar(tbl{row, 3}, data.rangeRows(row).mcend);
if ~isfinite(mcbegin) || ~isfinite(mcend) || mcend <= mcbegin
    msg = "Range limits must satisfy End > Begin.";
    return;
end
chargeState = data.rangeRows(row).chargeState;
try
    [~, csParsed, ~] = ionConvertName(char(rangeName));
    if isfinite(csParsed)
        chargeState = double(csParsed);
    end
catch
end
oldHandle = data.rangeRows(row).handle;

if rangeExistsInSpec(data.spec, rangeName, mcbegin, mcend, chargeState, oldHandle)
    ok = true;
    msg = "Range already enlisted; edit ignored.";
    return;
end

beforeRows = collectRangeRows(data.spec);
beforeHandles = [beforeRows.handle];

try
    if strcmpi(rangeName, 'not assigned') || strcmpi(rangeName, 'background')
        rangeAdd(data.spec, data.colorScheme, 'background', [mcbegin, mcend]);
    else
        data.colorScheme = ensureIonColor(data.colorScheme, rangeName);
        rangeAdd(data.spec, data.colorScheme, char(rangeName), [mcbegin, mcend]);
    end
catch ME
    msg = string(ME.message);
    return;
end

afterRows = collectRangeRows(data.spec);
afterHandles = [afterRows.handle];
newHandle = findAddedHandle(beforeHandles, afterHandles);
if isempty(newHandle)
    msg = "Could not resolve newly added range handle.";
    return;
end

try
    ud = newHandle.UserData;
    if ~isstruct(ud)
        ud = struct();
    end
    ud.chargeState = chargeState;
    newHandle.UserData = ud;
catch
end

oldName = data.rangeRows(row).name;
deleteRangeWithLabel(data.ax, oldHandle, oldName);

ok = true;
msg = sprintf('Range updated: %.3f to %.3f Da.', mcbegin, mcend);
end

function hAdded = findAddedHandle(oldHandles, newHandles)
hAdded = [];
if isempty(newHandles)
    return;
end
for i = 1:numel(newHandles)
    h = newHandles(i);
    isOld = false;
    for j = 1:numel(oldHandles)
        if h == oldHandles(j)
            isOld = true;
            break;
        end
    end
    if ~isOld
        hAdded = h;
        return;
    end
end
end

function colorScheme = setColorSchemeIonColor(colorScheme, ionName, rgb)
baseIon = normalizeIonForColorScheme(ionName);
idx = find(string(colorScheme.ion) == baseIon, 1, 'first');
if isempty(idx)
    colorScheme(end+1, :) = {categorical(baseIon), rgb};
else
    colorScheme.color(idx, :) = rgb;
end
end

function [c, ok] = parseColorTriplet(valueIn)
c = [0 0 0];
ok = false;
if isnumeric(valueIn)
    v = double(valueIn(:)');
elseif ischar(valueIn) || isstring(valueIn)
    txt = strtrim(char(string(valueIn)));
    if isempty(txt)
        return;
    end
    txt = strrep(txt, '[', '');
    txt = strrep(txt, ']', '');
    v = str2num(txt); %#ok<ST2NM>
else
    return;
end
if isempty(v)
    return;
end
v = double(v(:)');
if numel(v) ~= 3 || any(~isfinite(v)) || any(v < 0) || any(v > 1)
    return;
end
c = v;
ok = true;
end

function val = parseFiniteScalar(valueIn, fallback)
val = fallback;
if isnumeric(valueIn) && isscalar(valueIn) && isfinite(valueIn)
    val = double(valueIn);
    return;
end
if ischar(valueIn) || isstring(valueIn)
    x = str2double(strtrim(char(string(valueIn))));
    if isfinite(x)
        val = x;
    end
end
end

function val = parsePositiveScalar(valueIn, fallback)
val = parseFiniteScalar(valueIn, fallback);
if ~isfinite(val) || val <= 0
    if isfinite(fallback) && fallback > 0
        val = fallback;
    else
        val = 1;
    end
end
val = round(val);
end

function val = parseLogicalScalar(valueIn, fallback)
val = logical(fallback);
if islogical(valueIn) && isscalar(valueIn)
    val = valueIn;
    return;
end
if isnumeric(valueIn) && isscalar(valueIn)
    val = valueIn ~= 0;
    return;
end
if ischar(valueIn) || isstring(valueIn)
    txt = lower(strtrim(char(string(valueIn))));
    if any(strcmp(txt, {'1','true','yes','y','on'}))
        val = true;
    elseif any(strcmp(txt, {'0','false','no','n','off'}))
        val = false;
    end
end
end

function [rangeName, mcbegin, mcend, chargeState] = rangeIdentityFromHandle(hRange)
rangeName = "background";
mcbegin = NaN;
mcend = NaN;
chargeState = NaN;
if ~isgraphics(hRange)
    return;
end
try
    if isprop(hRange, 'DisplayName') && strlength(string(hRange.DisplayName)) > 0
        rangeName = string(hRange.DisplayName);
    end
catch
end
try
    if ~isempty(hRange.XData)
        mcbegin = double(hRange.XData(1));
        mcend = double(hRange.XData(end));
    end
catch
end
try
    ud = hRange.UserData;
    if isstruct(ud) && isfield(ud, 'chargeState')
        chargeState = double(ud.chargeState);
    end
catch
end
if ~isfinite(chargeState)
    chargeState = NaN;
end
end

function tf = rangeExistsInSpec(spec, rangeName, mcbegin, mcend, chargeState, ignoreHandle)
tf = false;
if nargin < 6
    ignoreHandle = [];
end
if isempty(spec) || ~isgraphics(spec)
    return;
end
rows = collectRangeRows(spec);
if isempty(rows)
    return;
end
targetName = lower(strtrim(char(string(rangeName))));
if isempty(targetName)
    targetName = 'background';
end
if ~isfinite(mcbegin) || ~isfinite(mcend)
    return;
end
tol = 1e-9;
for i = 1:numel(rows)
    if ~isempty(ignoreHandle) && isgraphics(ignoreHandle) && rows(i).handle == ignoreHandle
        continue;
    end
    rowName = lower(strtrim(char(string(rows(i).name))));
    if isempty(rowName)
        rowName = 'background';
    end
    if ~strcmp(rowName, targetName)
        continue;
    end
    if ~isfinite(rows(i).mcbegin) || ~isfinite(rows(i).mcend)
        continue;
    end
    if abs(rows(i).mcbegin - mcbegin) > tol || abs(rows(i).mcend - mcend) > tol
        continue;
    end
    rowCs = rows(i).chargeState;
    sameCharge = (isnan(rowCs) && isnan(chargeState)) || (~isnan(rowCs) && ~isnan(chargeState) && rowCs == chargeState);
    if sameCharge
        tf = true;
        return;
    end
end
end

function deleteRangeWithLabel(ax, rangeHandle, rangeName)
if isgraphics(rangeHandle)
    delete(rangeHandle);
end
if ~isgraphics(ax) || strlength(string(rangeName)) == 0
    return;
end
plots = ax.Children;
for i = 1:numel(plots)
    try
        if isfield(plots(i).UserData, 'plotType') && plots(i).UserData.plotType == "text"
            if isprop(plots(i), 'DisplayName') && string(plots(i).DisplayName) == string(rangeName)
                delete(plots(i));
            end
        end
    catch
    end
end
end

function bw = estimateBinWidthFromSpec(spec)
bw = 0.01;
if ~isgraphics(spec)
    return;
end
try
    x = spec.XData(:);
    if numel(x) < 2
        return;
    end
    dx = diff(x);
    dx = dx(isfinite(dx) & dx > 0);
    if ~isempty(dx)
        bw = median(dx);
    end
catch
end
end
