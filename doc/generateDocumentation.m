function generateDocumentation(options)
% GENERATEDOCUMENTATION Generate HTML documentation from function help text
%
% generateDocumentation()
% generateDocumentation('outputDir', 'doc/html')
%
% Scans all .m files in the toolbox and generates HTML documentation
% compatible with MATLAB's Help Browser.
%
% OPTIONS:
%   'outputDir'    - Output directory for HTML files (default: 'doc')
%   'includeSource'- Include source code in docs (default: false)
%   'verbose'      - Show progress (default: true)
%
% This function generates:
%   - Individual HTML pages for each function
%   - Category index pages
%   - Main landing page (GettingStarted.html)
%   - Updated helptoc.xml
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    options.outputDir (1,:) char = ''
    options.includeSource (1,1) logical = false
    options.verbose (1,1) logical = true
end

% Get toolbox root
toolboxRoot = fileparts(fileparts(mfilename('fullpath')));

if isempty(options.outputDir)
    options.outputDir = fullfile(toolboxRoot, 'doc');
end

% Create output directory if needed
if ~isfolder(options.outputDir)
    mkdir(options.outputDir);
end

%% Define function categories
categories = functionCategories();

%% Generate main landing page
if options.verbose
    fprintf('Generating main documentation page...\n');
end
generateMainPage(options.outputDir, toolboxRoot);

%% Generate function reference page
if options.verbose
    fprintf('Generating function reference...\n');
end
generateFunctionReference(options.outputDir, categories);

% Update API_Reference.md
if options.verbose
    fprintf('Generating API reference...\n');
end
generateApiReference('outputFile', fullfile(options.outputDir, 'API_Reference.md'), 'categories', categories, 'verbose', options.verbose);

%% Generate individual function pages
categoryNames = fieldnames(categories);
totalFunctions = 0;
generatedFunctions = 0;

for c = 1:length(categoryNames)
    catName = categoryNames{c};
    funcs = categories.(catName).functions;
    totalFunctions = totalFunctions + length(funcs);
end

if options.verbose
    fprintf('Generating %d function pages...\n', totalFunctions);
end

for c = 1:length(categoryNames)
    catName = categoryNames{c};
    cat = categories.(catName);

    % Generate category index page
    generateCategoryPage(options.outputDir, catName, cat, toolboxRoot);

    % Generate individual function pages
    for f = 1:length(cat.functions)
        funcName = cat.functions{f};
        success = generateFunctionPage(options.outputDir, funcName, toolboxRoot, options);
        if success
            generatedFunctions = generatedFunctions + 1;
        end
    end
end

if options.verbose
    fprintf('Generated %d of %d function pages.\n', generatedFunctions, totalFunctions);
    fprintf('Documentation generated in: %s\n', options.outputDir);
end

end

%% Helper Functions

function generateMainPage(outputDir, toolboxRoot)
    % Generate GettingStarted.html

    html = ['<!DOCTYPE html>\n' ...
            '<html>\n' ...
            '<head>\n' ...
            '  <meta charset="utf-8">\n' ...
            '  <title>Atom Probe Toolbox</title>\n' ...
            '  <link rel="stylesheet" href="matlab:helpwin(''helpwin'')">\n' ...
            '</head>\n' ...
            '<body>\n' ...
            '<h1>Atom Probe Toolbox</h1>\n' ...
            '<p>A comprehensive MATLAB toolbox for Atom Probe Tomography (APT) data analysis.</p>\n' ...
            '\n' ...
            '<h2>Getting Started</h2>\n' ...
            '<ul>\n' ...
            '  <li><a href="matlab:setupToolbox">Run setupToolbox</a> to initialize the toolbox</li>\n' ...
            '  <li><a href="matlab:checkDependencies">Run checkDependencies</a> to verify installation</li>\n' ...
            '  <li><a href="matlab:open(''Workflow_FirstSteps.mlx'')">Open First Steps Workflow</a></li>\n' ...
            '</ul>\n' ...
            '\n' ...
            '<h2>Quick Start</h2>\n' ...
            '<pre>\n' ...
            '%% Initialize\n' ...
            'setupToolbox();\n' ...
            '\n' ...
            '%% Load data\n' ...
            'pos = posLoad(''mydata.pos'');\n' ...
            '\n' ...
            '%% Create mass spectrum\n' ...
            'massSpecPlot(pos.mc);\n' ...
            '\n' ...
            '%% Define ions\n' ...
            'ions = ionAdd(''Fe'');\n' ...
            'ions = ionAdd(''Cr'', ions);\n' ...
            '</pre>\n' ...
            '\n' ...
            '<h2>Documentation</h2>\n' ...
            '<ul>\n' ...
            '  <li><a href="function_reference.html">Function Reference</a></li>\n' ...
            '  <li><a href="matlab:open(''doc/API_Reference.md'')">API Quick Reference</a></li>\n' ...
            '  <li><a href="matlab:open(''doc/TROUBLESHOOTING.md'')">Troubleshooting Guide</a></li>\n' ...
            '</ul>\n' ...
            '\n' ...
            '<h2>Workflows</h2>\n' ...
            '<ul>\n' ...
            '  <li><a href="matlab:open(''Workflow_FirstSteps.mlx'')">First Steps</a></li>\n' ...
            '  <li><a href="matlab:open(''Workflow_1DConcentrationProfile.mlx'')">1D Concentration Profile</a></li>\n' ...
            '  <li><a href="matlab:open(''Workflow_ClusterDetermination.mlx'')">Cluster Determination</a></li>\n' ...
            '  <li><a href="matlab:open(''Workflow_Proxigram.mlx'')">Proxigram Analysis</a></li>\n' ...
            '  <li><a href="matlab:open(''Workflow_HDF5_IO.mlx'')">HDF5 Database</a></li>\n' ...
            '</ul>\n' ...
            '\n' ...
            '<h2>External Resources</h2>\n' ...
            '<ul>\n' ...
            '  <li><a href="https://github.com/peterfelfer/Atom-Probe-Toolbox">GitHub Repository</a></li>\n' ...
            '  <li><a href="https://www.youtube.com/playlist?list=PLLr-VNeczShDzrm-wdIyz9ORjFL5KDeYT">Video Tutorials</a></li>\n' ...
            '</ul>\n' ...
            '\n' ...
            '<p><em>(c) Prof. Peter Felfer Group @ FAU Erlangen-Nurnberg</em></p>\n' ...
            '</body>\n' ...
            '</html>\n'];

    fid = fopen(fullfile(outputDir, 'GettingStarted.html'), 'w');
    fprintf(fid, '%s', html);
    fclose(fid);
end

function generateFunctionReference(outputDir, categories)
    % Generate function_reference.html

    html = ['<!DOCTYPE html>\n' ...
            '<html>\n' ...
            '<head>\n' ...
            '  <meta charset="utf-8">\n' ...
            '  <title>Function Reference - Atom Probe Toolbox</title>\n' ...
            '</head>\n' ...
            '<body>\n' ...
            '<h1>Function Reference</h1>\n' ...
            '<p>Functions organized by category.</p>\n\n'];

    categoryNames = fieldnames(categories);
    for c = 1:length(categoryNames)
        catName = categoryNames{c};
        cat = categories.(catName);

        html = [html, sprintf('<h2 id="%s">%s</h2>\n', catName, cat.name)];
        html = [html, '<table border="1" cellpadding="5">\n'];
        html = [html, '<tr><th>Function</th><th>Description</th></tr>\n'];

        for f = 1:length(cat.functions)
            funcName = cat.functions{f};
            desc = getFunctionDescription(funcName);
            html = [html, sprintf('<tr><td><a href="%s.html">%s</a></td><td>%s</td></tr>\n', ...
                funcName, funcName, desc)];
        end

        html = [html, '</table>\n\n'];
    end

    html = [html, '</body>\n</html>\n'];

    fid = fopen(fullfile(outputDir, 'function_reference.html'), 'w');
    fprintf(fid, '%s', html);
    fclose(fid);
end

function generateCategoryPage(outputDir, catName, cat, toolboxRoot)
    % Generate category index page

    html = ['<!DOCTYPE html>\n' ...
            '<html>\n' ...
            '<head>\n' ...
            '  <meta charset="utf-8">\n' ...
            sprintf('  <title>%s - Atom Probe Toolbox</title>\n', cat.name) ...
            '</head>\n' ...
            '<body>\n' ...
            sprintf('<h1>%s</h1>\n', cat.name) ...
            '<table border="1" cellpadding="5">\n' ...
            '<tr><th>Function</th><th>Description</th></tr>\n'];

    for f = 1:length(cat.functions)
        funcName = cat.functions{f};
        desc = getFunctionDescription(funcName);
        html = [html, sprintf('<tr><td><a href="%s.html">%s</a></td><td>%s</td></tr>\n', ...
            funcName, funcName, desc)];
    end

    html = [html, '</table>\n' ...
            '<p><a href="function_reference.html">Back to Function Reference</a></p>\n' ...
            '</body>\n</html>\n'];

    fid = fopen(fullfile(outputDir, ['functions_' catName '.html']), 'w');
    fprintf(fid, '%s', html);
    fclose(fid);
end

function success = generateFunctionPage(outputDir, funcName, toolboxRoot, options)
    % Generate HTML page for a single function

    success = false;

    % Try to get help text
    try
        helpText = help(funcName);
        if isempty(helpText)
            return;
        end
    catch
        return;
    end

    % Convert help text to HTML
    helpHtml = helpTextToHtml(helpText);

    html = ['<!DOCTYPE html>\n' ...
            '<html>\n' ...
            '<head>\n' ...
            '  <meta charset="utf-8">\n' ...
            sprintf('  <title>%s - Atom Probe Toolbox</title>\n', funcName) ...
            '</head>\n' ...
            '<body>\n' ...
            sprintf('<h1>%s</h1>\n', funcName) ...
            '<pre>\n' ...
            helpHtml ...
            '</pre>\n' ...
            '<p><a href="function_reference.html">Back to Function Reference</a></p>\n' ...
            '</body>\n</html>\n'];

    fid = fopen(fullfile(outputDir, [funcName '.html']), 'w');
    if fid == -1
        return;
    end
    fprintf(fid, '%s', html);
    fclose(fid);

    success = true;
end

function desc = getFunctionDescription(funcName)
    % Get the H1 line (first line of help) for a function

    desc = '';
    try
        helpText = help(funcName);
        if ~isempty(helpText)
            lines = strsplit(helpText, '\n');
            if ~isempty(lines)
                % First non-empty line after function name
                for i = 1:min(3, length(lines))
                    line = strtrim(lines{i});
                    if ~isempty(line) && ~startsWith(upper(line), upper(funcName))
                        desc = line;
                        break;
                    end
                end
            end
        end
    catch
        desc = '';
    end

    % Escape HTML characters
    desc = strrep(desc, '&', '&amp;');
    desc = strrep(desc, '<', '&lt;');
    desc = strrep(desc, '>', '&gt;');
end

function html = helpTextToHtml(helpText)
    % Convert MATLAB help text to HTML

    % Escape HTML characters
    html = strrep(helpText, '&', '&amp;');
    html = strrep(html, '<', '&lt;');
    html = strrep(html, '>', '&gt;');

    % Convert URLs to links
    html = regexprep(html, '(https?://[^\s]+)', '<a href="$1">$1</a>');

    % Convert function references to links
    html = regexprep(html, '\b([a-zA-Z]\w+)\s*-\s*', '<a href="matlab:doc $1">$1</a> - ');
end
