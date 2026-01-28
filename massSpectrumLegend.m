function legendFig = massSpectrumLegend(spec, varargin)
% massSpectrumLegend creates a publication-ready legend for a mass spectrum
%
% legendFig = massSpectrumLegend(spec)
% legendFig = massSpectrumLegend(spec, Name, Value, ...)
%
% Creates a legend in a separate figure window with:
% - Colored boxes for each ion/range in the mass spectrum
% - Line samples for other plot elements (background curves, etc.)
% - Ion names formatted in LaTeX with proper sub/superscripts
%
% INPUT
% spec:     area plot of the mass spectrum, or axes containing it
%
% Name-Value Options:
%   'columns'           Number of columns (default: 2)
%   'boxSize'           [width height] of color boxes in points (default: [20 12])
%   'fontSize'          Font size for labels (default: 11)
%   'fontName'          Font name (default: 'Helvetica')
%   'sortBy'            Sort order: 'mass', 'name', or 'none' (default: 'mass')
%   'includeTheoretical' Include theoretical ion markers (default: false)
%   'includeBackground' Include background curves (default: true)
%   'title'             Legend title (default: '' for no title)
%   'spacing'           Vertical spacing between entries in points (default: 6)
%   'margin'            Margin around legend in points (default: 15)
%   'useUnicode'        Use Unicode sub/superscripts instead of LaTeX (default: false)
%                       Set to true for editable text in SVG exports
%   'output'            File path to save directly (e.g., 'legend.svg')
%                       Supported formats: .svg, .pdf, .png, .eps
%                       If provided, saves and closes figure automatically
%
% OUTPUT
% legendFig:    Handle to the legend figure
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% Parse options
options = struct( ...
    'columns', 2, ...
    'boxSize', [20 12], ...
    'fontSize', 11, ...
    'fontName', 'Helvetica', ...
    'sortBy', 'mass', ...
    'includeTheoretical', false, ...
    'includeBackground', true, ...
    'title', '', ...
    'spacing', 6, ...
    'margin', 15, ...
    'useUnicode', false, ...  % true = editable text in SVG, false = LaTeX (prettier but text as paths)
    'output', '');            % if provided, save directly to this path (SVG/PDF)

options = parseOptions(options, varargin{:});

%% Resolve axes
if isgraphics(spec, 'axes')
    ax = spec;
else
    ax = spec.Parent;
end

%% Extract range information (ions)
plots = ax.Children;
ionEntries = struct('name', {}, 'latexName', {}, 'color', {}, 'mcbegin', {});

for i = 1:numel(plots)
    try
        plotType = plots(i).UserData.plotType;
    catch
        continue;
    end

    if plotType == "range"
        % Get ion information
        entry = struct();
        entry.mcbegin = plots(i).XData(1);
        entry.color = plots(i).FaceColor;

        % Get name and convert to LaTeX
        if isfield(plots(i).UserData, 'ion') && istable(plots(i).UserData.ion)
            if isfield(plots(i).UserData, 'chargeState')
                chargeState = plots(i).UserData.chargeState;
            else
                chargeState = NaN;
            end
            entry.name = ionConvertName(plots(i).UserData.ion, chargeState, 'plain');
            entry.latexName = ionConvertName(plots(i).UserData.ion, chargeState, 'LaTeX');
        elseif isfield(plots(i).UserData, 'ion')
            entry.name = string(plots(i).UserData.ion);
            entry.latexName = string(plots(i).UserData.ion);
        elseif isprop(plots(i), 'DisplayName') && ~isempty(plots(i).DisplayName)
            entry.name = plots(i).DisplayName;
            entry.latexName = plots(i).DisplayName;
        else
            entry.name = 'Unknown';
            entry.latexName = 'Unknown';
        end

        ionEntries(end+1) = entry;
    end
end

%% Extract theoretical ion markers (optional)
theoreticalEntries = struct('name', {}, 'latexName', {}, 'color', {}, 'mcbegin', {});

if options.includeTheoretical
    for i = 1:numel(plots)
        try
            plotType = plots(i).UserData.plotType;
        catch
            continue;
        end

        if plotType == "ion"
            entry = struct();
            entry.color = plots(i).Color;

            % Get first peak position for sorting
            entry.mcbegin = min(plots(i).XData);

            if isprop(plots(i), 'DisplayName') && ~isempty(plots(i).DisplayName)
                entry.name = plots(i).DisplayName;
                % Try to convert to LaTeX if it looks like an ion name
                try
                    [ionTable, cs] = ionConvertName(plots(i).DisplayName);
                    entry.latexName = ionConvertName(ionTable, cs, 'LaTeX');
                catch
                    entry.latexName = plots(i).DisplayName;
                end
            else
                entry.name = 'Ion';
                entry.latexName = 'Ion';
            end

            theoreticalEntries(end+1) = entry;
        end
    end
end

%% Extract line plots (background, etc.)
lineEntries = struct('name', {}, 'color', {}, 'lineStyle', {}, 'lineWidth', {});

if options.includeBackground
    for i = 1:numel(plots)
        try
            plotType = plots(i).UserData.plotType;
        catch
            continue;
        end

        if plotType == "backgroundEstimate" || plotType == "background"
            entry = struct();
            entry.color = plots(i).Color;
            entry.lineStyle = plots(i).LineStyle;
            entry.lineWidth = plots(i).LineWidth;

            if isprop(plots(i), 'DisplayName') && ~isempty(plots(i).DisplayName)
                entry.name = plots(i).DisplayName;
            else
                entry.name = 'Background';
            end

            lineEntries(end+1) = entry;
        end
    end
end

%% Consolidate ion entries by ion composition (ignoring isotopes and charge states)
if ~isempty(ionEntries)
    % Compute base ion name for each entry (element composition only)
    baseIonNames = cell(numel(ionEntries), 1);
    for i = 1:numel(ionEntries)
        try
            origPlot = findRangePlotByMcbegin(plots, ionEntries(i).mcbegin);
            if ~isempty(origPlot) && isfield(origPlot.UserData, 'ion')
                ionData = origPlot.UserData.ion;
                if istable(ionData) && ismember('element', ionData.Properties.VariableNames)
                    % Get base ion name (elements + counts, no isotopes, no charge)
                    baseIonNames{i} = ionConvertName(ionData.element, NaN, 'plain');
                else
                    baseIonNames{i} = ionEntries(i).name;
                end
            else
                baseIonNames{i} = ionEntries(i).name;
            end
        catch
            baseIonNames{i} = ionEntries(i).name;
        end
    end

    % Group by base ion name
    [uniqueBaseNames, uniqueIdx] = unique(baseIonNames, 'stable');

    % Create consolidated entries
    consolidatedEntries = struct('name', {}, 'latexName', {}, 'color', {}, 'mcbegin', {});
    for i = 1:numel(uniqueBaseNames)
        idx = uniqueIdx(i);
        entry = ionEntries(idx);

        % Get LaTeX version of the base ion name
        try
            origPlot = findRangePlotByMcbegin(plots, ionEntries(idx).mcbegin);
            if ~isempty(origPlot) && isfield(origPlot.UserData, 'ion')
                ionData = origPlot.UserData.ion;
                if istable(ionData) && ismember('element', ionData.Properties.VariableNames)
                    entry.name = strtrim(ionConvertName(ionData.element, NaN, 'plain'));
                    entry.latexName = strtrim(ionConvertName(ionData.element, NaN, 'LaTeX'));
                end
            end
        catch
            % Keep original names
        end

        % Use minimum mcbegin among all entries with this base ion name
        sameIonIdx = find(strcmp(baseIonNames, uniqueBaseNames{i}));
        entry.mcbegin = min([ionEntries(sameIonIdx).mcbegin]);

        consolidatedEntries(end+1) = entry;
    end
    ionEntries = consolidatedEntries;
end

if ~isempty(theoreticalEntries)
    [~, uniqueIdx] = unique({theoreticalEntries.name}, 'stable');
    theoreticalEntries = theoreticalEntries(uniqueIdx);
end

%% Sort entries
if ~isempty(ionEntries)
    switch lower(options.sortBy)
        case 'mass'
            [~, sortIdx] = sort([ionEntries.mcbegin]);
            ionEntries = ionEntries(sortIdx);
        case 'name'
            [~, sortIdx] = sort({ionEntries.name});
            ionEntries = ionEntries(sortIdx);
        % 'none' - keep original order
    end
end

if ~isempty(theoreticalEntries)
    switch lower(options.sortBy)
        case 'mass'
            [~, sortIdx] = sort([theoreticalEntries.mcbegin]);
            theoreticalEntries = theoreticalEntries(sortIdx);
        case 'name'
            [~, sortIdx] = sort({theoreticalEntries.name});
            theoreticalEntries = theoreticalEntries(sortIdx);
    end
end

%% Calculate layout
nIons = numel(ionEntries);
nTheoretical = numel(theoreticalEntries);
nLines = numel(lineEntries);
nTotal = nIons + nTheoretical + nLines;

if nTotal == 0
    warning('massSpectrumLegend:noEntries', 'No legend entries found in mass spectrum.');
    legendFig = [];
    return;
end

nCols = min(options.columns, nTotal);
nRows = ceil(nTotal / nCols);

% Dimensions in points
boxW = options.boxSize(1);
boxH = options.boxSize(2);
fontSize = options.fontSize;
spacing = options.spacing;
margin = options.margin;
textGap = 5;  % gap between box and text

% Estimate text width (rough estimate based on font size and typical name length)
maxNameLen = 0;
for i = 1:nIons
    maxNameLen = max(maxNameLen, strlength(ionEntries(i).name));
end
for i = 1:nTheoretical
    maxNameLen = max(maxNameLen, strlength(theoreticalEntries(i).name));
end
for i = 1:nLines
    maxNameLen = max(maxNameLen, strlength(lineEntries(i).name));
end
estTextWidth = maxNameLen * fontSize * 0.6;  % rough estimate

colWidth = boxW + textGap + estTextWidth + spacing;
rowHeight = max(boxH, fontSize * 1.2) + spacing;

% Figure size
figWidth = nCols * colWidth + 2 * margin;
figHeight = nRows * rowHeight + 2 * margin;
if ~isempty(options.title)
    figHeight = figHeight + fontSize * 1.5;
end

%% Create legend figure
legendFig = figure('Name', 'Mass Spectrum Legend', ...
    'NumberTitle', 'off', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none', ...
    'Units', 'points', ...
    'Position', [100 100 figWidth figHeight], ...
    'Color', 'white', ...
    'Resize', 'on');

% Create axes for drawing
legendAx = axes(legendFig, ...
    'Units', 'points', ...
    'Position', [0 0 figWidth figHeight], ...
    'XLim', [0 figWidth], ...
    'YLim', [0 figHeight], ...
    'Visible', 'off', ...
    'YDir', 'reverse');  % Top-to-bottom layout
hold(legendAx, 'on');

%% Draw title if specified
yOffset = margin;
if ~isempty(options.title)
    text(legendAx, figWidth/2, yOffset, options.title, ...
        'FontSize', fontSize + 2, ...
        'FontName', options.fontName, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
    yOffset = yOffset + fontSize * 1.5;
end

%% Draw ion entries (colored boxes)
entryIdx = 0;
for i = 1:nIons
    entryIdx = entryIdx + 1;
    col = mod(entryIdx - 1, nCols);
    row = floor((entryIdx - 1) / nCols);

    x = margin + col * colWidth;
    y = yOffset + row * rowHeight;

    % Draw colored box
    rectangle(legendAx, ...
        'Position', [x, y, boxW, boxH], ...
        'FaceColor', ionEntries(i).color, ...
        'EdgeColor', 'k', ...
        'LineWidth', 0.5);

    % Draw label
    if options.useUnicode
        labelText = latexToUnicode(ionEntries(i).latexName);
        text(legendAx, x + boxW + textGap, y + boxH/2, labelText, ...
            'FontSize', fontSize, ...
            'FontName', options.fontName, ...
            'Interpreter', 'none', ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left');
    else
        text(legendAx, x + boxW + textGap, y + boxH/2, ...
            ['$' ionEntries(i).latexName '$'], ...
            'FontSize', fontSize, ...
            'FontName', options.fontName, ...
            'Interpreter', 'latex', ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left');
    end
end

%% Draw theoretical ion entries (colored boxes with dashed border)
for i = 1:nTheoretical
    entryIdx = entryIdx + 1;
    col = mod(entryIdx - 1, nCols);
    row = floor((entryIdx - 1) / nCols);

    x = margin + col * colWidth;
    y = yOffset + row * rowHeight;

    % Draw colored box with dashed border
    rectangle(legendAx, ...
        'Position', [x, y, boxW, boxH], ...
        'FaceColor', theoreticalEntries(i).color, ...
        'EdgeColor', 'k', ...
        'LineStyle', '--', ...
        'LineWidth', 0.5);

    % Draw label
    if options.useUnicode
        labelText = latexToUnicode(theoreticalEntries(i).latexName);
        text(legendAx, x + boxW + textGap, y + boxH/2, labelText, ...
            'FontSize', fontSize, ...
            'FontName', options.fontName, ...
            'Interpreter', 'none', ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left');
    else
        text(legendAx, x + boxW + textGap, y + boxH/2, ...
            ['$' theoreticalEntries(i).latexName '$'], ...
            'FontSize', fontSize, ...
            'FontName', options.fontName, ...
            'Interpreter', 'latex', ...
            'VerticalAlignment', 'middle', ...
            'HorizontalAlignment', 'left');
    end
end

%% Draw line entries (line samples)
for i = 1:nLines
    entryIdx = entryIdx + 1;
    col = mod(entryIdx - 1, nCols);
    row = floor((entryIdx - 1) / nCols);

    x = margin + col * colWidth;
    y = yOffset + row * rowHeight;

    % Draw line sample
    line(legendAx, [x, x + boxW], [y + boxH/2, y + boxH/2], ...
        'Color', lineEntries(i).color, ...
        'LineStyle', lineEntries(i).lineStyle, ...
        'LineWidth', lineEntries(i).lineWidth);

    % Draw label
    text(legendAx, x + boxW + textGap, y + boxH/2, ...
        lineEntries(i).name, ...
        'FontSize', fontSize, ...
        'FontName', options.fontName, ...
        'Interpreter', 'none', ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left');
end

hold(legendAx, 'off');

%% If output path provided, save directly and optionally close figure
if ~isempty(options.output)
    outputPath = options.output;
    [~, ~, ext] = fileparts(outputPath);

    drawnow;  % Ensure figure is fully rendered

    try
        if strcmpi(ext, '.svg')
            exportgraphics(legendFig, outputPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        elseif strcmpi(ext, '.pdf')
            exportgraphics(legendFig, outputPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        elseif strcmpi(ext, '.png')
            exportgraphics(legendFig, outputPath, 'Resolution', 300, 'BackgroundColor', 'white');
        elseif strcmpi(ext, '.eps')
            exportgraphics(legendFig, outputPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        else
            % Default to SVG if no extension or unknown
            if isempty(ext)
                outputPath = [outputPath '.svg'];
            end
            exportgraphics(legendFig, outputPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        end
        fprintf('Legend saved to: %s\n', outputPath);
    catch ME
        % Fallback for older MATLAB
        try
            if strcmpi(ext, '.svg') || isempty(ext)
                print(legendFig, outputPath, '-dsvg');
            elseif strcmpi(ext, '.pdf')
                print(legendFig, outputPath, '-dpdf', '-bestfit');
            elseif strcmpi(ext, '.png')
                print(legendFig, outputPath, '-dpng', '-r300');
            elseif strcmpi(ext, '.eps')
                print(legendFig, outputPath, '-depsc');
            end
            fprintf('Legend saved to: %s\n', outputPath);
        catch
            warning('massSpectrumLegend:saveFailed', 'Could not save: %s', ME.message);
        end
    end

    % Close figure if output was saved (return handle anyway)
    close(legendFig);
    legendFig = [];
    return;
end

%% Add export buttons
btnHeight = 20;
btnWidth = 80;
btnMargin = 5;
btnGap = 10;

% Resize figure to make room for buttons
figPos = legendFig.Position;
figPos(4) = figPos(4) + btnHeight + btnMargin * 2;
legendFig.Position = figPos;

% Move axes up to make room for buttons at bottom
axPos = legendAx.Position;
axPos(4) = axPos(4) + btnHeight + btnMargin * 2;
legendAx.Position = axPos;
legendAx.YLim = [0 axPos(4)];

% Calculate button positions (centered)
totalBtnWidth = btnWidth * 3 + btnGap * 2;
btnStartX = (figPos(3) - totalBtnWidth) / 2;

% Add "Copy" button (EMF on Windows, PDF on macOS)
if ismac
    copyTooltip = 'Copy to clipboard as PDF vector';
else
    copyTooltip = 'Copy to clipboard as EMF vector';
end
copyBtn = uicontrol(legendFig, ...
    'Style', 'pushbutton', ...
    'String', 'Copy', ...
    'Units', 'points', ...
    'Position', [btnStartX, btnMargin, btnWidth, btnHeight], ...
    'Callback', @(~,~) copyAsVector(legendFig), ...
    'TooltipString', copyTooltip);

% Add "Save SVG" button (works with Inkscape)
svgBtn = uicontrol(legendFig, ...
    'Style', 'pushbutton', ...
    'String', 'Save SVG', ...
    'Units', 'points', ...
    'Position', [btnStartX + btnWidth + btnGap, btnMargin, btnWidth, btnHeight], ...
    'Callback', @(~,~) saveAsSVG(legendFig), ...
    'TooltipString', 'Save as SVG file (Inkscape, web)');

% Add "Save PDF" button (universal)
pdfBtn = uicontrol(legendFig, ...
    'Style', 'pushbutton', ...
    'String', 'Save PDF', ...
    'Units', 'points', ...
    'Position', [btnStartX + 2*(btnWidth + btnGap), btnMargin, btnWidth, btnHeight], ...
    'Callback', @(~,~) saveAsPDF(legendFig), ...
    'TooltipString', 'Save as PDF file (universal vector format)');

end

function copyAsVector(fig)
    % Copy figure to clipboard as vector graphic
    % macOS: copies as PDF, Windows: copies as EMF
    btns = findobj(fig, 'Type', 'uicontrol');
    set(btns, 'Visible', 'off');
    drawnow;

    try
        copygraphics(fig, 'ContentType', 'vector', 'BackgroundColor', 'white');
        if ismac
            disp('Legend copied to clipboard as PDF vector.');
        else
            disp('Legend copied to clipboard as EMF vector.');
        end
    catch ME
        if contains(ME.identifier, 'MATLAB:UndefinedFunction') || contains(ME.message, 'copygraphics')
            if ispc
                try
                    print(fig, '-clipboard', '-dmeta');
                    disp('Legend copied to clipboard as EMF.');
                catch
                    warning('massSpectrumLegend:copyFailed', ...
                        'Clipboard copy failed. Use Save SVG or Save PDF instead.');
                end
            else
                warning('massSpectrumLegend:copyFailed', ...
                    'Clipboard copy not available. Use Save SVG or Save PDF instead.');
            end
        else
            warning('massSpectrumLegend:copyFailed', 'Copy failed: %s', ME.message);
        end
    end

    set(btns, 'Visible', 'on');
end

function saveAsSVG(fig)
    % Save figure as SVG file
    [filename, filepath] = uiputfile('*.svg', 'Save Legend as SVG', 'legend.svg');
    if isequal(filename, 0)
        return;  % User cancelled
    end

    fullpath = fullfile(filepath, filename);

    btns = findobj(fig, 'Type', 'uicontrol');
    set(btns, 'Visible', 'off');
    drawnow;

    try
        % exportgraphics was introduced in R2020a
        exportgraphics(fig, fullpath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        fprintf('Legend saved to: %s\n', fullpath);
    catch ME
        if contains(ME.identifier, 'MATLAB:UndefinedFunction')
            % Fallback for older MATLAB
            try
                print(fig, fullpath, '-dsvg');
                fprintf('Legend saved to: %s\n', fullpath);
            catch
                warning('massSpectrumLegend:saveFailed', 'Could not save SVG: %s', ME.message);
            end
        else
            warning('massSpectrumLegend:saveFailed', 'Save failed: %s', ME.message);
        end
    end

    set(btns, 'Visible', 'on');
end

function saveAsPDF(fig)
    % Save figure as PDF file
    [filename, filepath] = uiputfile('*.pdf', 'Save Legend as PDF', 'legend.pdf');
    if isequal(filename, 0)
        return;  % User cancelled
    end

    fullpath = fullfile(filepath, filename);

    btns = findobj(fig, 'Type', 'uicontrol');
    set(btns, 'Visible', 'off');
    drawnow;

    try
        exportgraphics(fig, fullpath, 'ContentType', 'vector', 'BackgroundColor', 'white');
        fprintf('Legend saved to: %s\n', fullpath);
    catch ME
        if contains(ME.identifier, 'MATLAB:UndefinedFunction')
            try
                print(fig, fullpath, '-dpdf', '-bestfit');
                fprintf('Legend saved to: %s\n', fullpath);
            catch
                warning('massSpectrumLegend:saveFailed', 'Could not save PDF: %s', ME.message);
            end
        else
            warning('massSpectrumLegend:saveFailed', 'Save failed: %s', ME.message);
        end
    end

    set(btns, 'Visible', 'on');
end

%% Helper function to parse options
function options = parseOptions(options, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('massSpectrumLegend:invalidOptions', 'Options must be name-value pairs.');
    end

    for k = 1:2:numel(varargin)
        name = lower(string(varargin{k}));
        value = varargin{k+1};

        switch name
            case "columns"
                options.columns = value;
            case "boxsize"
                options.boxSize = value;
            case "fontsize"
                options.fontSize = value;
            case "fontname"
                options.fontName = value;
            case "sortby"
                options.sortBy = value;
            case "includetheoretical"
                options.includeTheoretical = logical(value);
            case "includebackground"
                options.includeBackground = logical(value);
            case "title"
                options.title = char(value);
            case "spacing"
                options.spacing = value;
            case "margin"
                options.margin = value;
            case "useunicode"
                options.useUnicode = logical(value);
            case "output"
                options.output = char(value);
            otherwise
                warning('massSpectrumLegend:unknownOption', 'Unknown option: %s', name);
        end
    end
end

function rangePlot = findRangePlotByMcbegin(plots, targetMcbegin)
    % Find a range plot with the specified mcbegin (start position)
    rangePlot = [];
    for i = 1:numel(plots)
        try
            if plots(i).UserData.plotType == "range"
                if abs(plots(i).XData(1) - targetMcbegin) < 1e-6
                    rangePlot = plots(i);
                    return;
                end
            end
        catch
        end
    end
end

function unicodeStr = latexToUnicode(latexStr)
    % Convert LaTeX-style ion name to Unicode with super/subscripts
    % e.g., "^{56}Fe_{2} ^{+}" -> "⁵⁶Fe₂⁺"

    % Unicode superscript digits
    superDigits = {'⁰','¹','²','³','⁴','⁵','⁶','⁷','⁸','⁹'};
    % Unicode subscript digits
    subDigits = {'₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'};

    unicodeStr = latexStr;

    % Convert superscripts: ^{...}
    pattern = '\^\{([^}]*)\}';
    matches = regexp(unicodeStr, pattern, 'tokens');
    for i = 1:numel(matches)
        content = matches{i}{1};
        unicodeContent = '';
        for j = 1:length(content)
            c = content(j);
            if c >= '0' && c <= '9'
                unicodeContent = [unicodeContent superDigits{str2double(c) + 1}];
            elseif c == '+'
                unicodeContent = [unicodeContent '⁺'];
            elseif c == '-'
                unicodeContent = [unicodeContent '⁻'];
            else
                unicodeContent = [unicodeContent c];
            end
        end
        unicodeStr = regexprep(unicodeStr, '\^\{[^}]*\}', unicodeContent, 'once');
    end

    % Convert subscripts: _{...}
    pattern = '_\{([^}]*)\}';
    matches = regexp(unicodeStr, pattern, 'tokens');
    for i = 1:numel(matches)
        content = matches{i}{1};
        unicodeContent = '';
        for j = 1:length(content)
            c = content(j);
            if c >= '0' && c <= '9'
                unicodeContent = [unicodeContent subDigits{str2double(c) + 1}];
            else
                unicodeContent = [unicodeContent c];
            end
        end
        unicodeStr = regexprep(unicodeStr, '_\{[^}]*\}', unicodeContent, 'once');
    end

    % Clean up any remaining whitespace
    unicodeStr = strtrim(regexprep(unicodeStr, '\s+', ' '));
end
