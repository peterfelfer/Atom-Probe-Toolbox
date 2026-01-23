function results = batchProcess(dataFiles, analysisFcn, options)
% BATCHPROCESS Apply analysis pipeline to multiple APT datasets
%
% results = batchProcess(dataFiles, analysisFcn)
% results = batchProcess(dataFiles, analysisFcn, 'parallel', true)
%
% Processes multiple APT datasets through a common analysis pipeline,
% collecting results and handling errors gracefully.
%
% INPUT:
%   dataFiles   - Cell array of file paths OR folder path with wildcard
%                 e.g., {'file1.pos', 'file2.pos'} or '/data/*.pos'
%   analysisFcn - Function handle that takes a file path and returns results
%                 Signature: result = analysisFcn(filePath)
%                 OR: result = analysisFcn(filePath, idx, total) for progress info
%
% OPTIONS:
%   'parallel'      - Use parallel processing (default: false)
%   'continueOnError' - Continue if one file fails (default: true)
%   'saveResults'   - Save results after each file (default: false)
%   'outputFolder'  - Folder for saving results (default: pwd)
%   'outputPrefix'  - Prefix for output files (default: 'batch_')
%   'showProgress'  - Show progress bar (default: true)
%   'logFile'       - Path to log file (default: '' = no logging)
%
% OUTPUT:
%   results - Structure array with fields:
%       .fileName   - Original file name
%       .filePath   - Full file path
%       .success    - Whether analysis succeeded
%       .error      - Error message if failed
%       .data       - Analysis results (from analysisFcn)
%       .duration   - Processing time in seconds
%
% EXAMPLES:
%   % Simple concentration analysis
%   files = {'sample1.pos', 'sample2.pos', 'sample3.pos'};
%   results = batchProcess(files, @(f) posCalculateConcentrationSimple(posLoad(f)));
%
%   % With custom analysis function
%   myAnalysis = @(f) analyzeDataset(f, 'binWidth', 0.1);
%   results = batchProcess('/data/*.pos', myAnalysis, 'parallel', true);
%
%   % Using index information for logging
%   myAnalysis = @(f, idx, total) fprintf('Processing %d/%d\n', idx, total);
%   results = batchProcess(files, myAnalysis);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    dataFiles
    analysisFcn function_handle
    options.parallel (1,1) logical = false
    options.continueOnError (1,1) logical = true
    options.saveResults (1,1) logical = false
    options.outputFolder (1,:) char = pwd
    options.outputPrefix (1,:) char = 'batch_'
    options.showProgress (1,1) logical = true
    options.logFile (1,:) char = ''
end

% Resolve file list
if ischar(dataFiles) || isstring(dataFiles)
    dataFiles = char(dataFiles);
    if contains(dataFiles, '*') || contains(dataFiles, '?')
        % Wildcard pattern
        [folder, pattern, ext] = fileparts(dataFiles);
        if isempty(folder)
            folder = '.';
        end
        listing = dir(fullfile(folder, [pattern, ext]));
        dataFiles = cellfun(@(f) fullfile(folder, f), {listing.name}, 'UniformOutput', false);
    else
        dataFiles = {dataFiles};
    end
end

nFiles = length(dataFiles);

if nFiles == 0
    warning('batchProcess:noFiles', 'No files found to process.');
    results = struct([]);
    return;
end

% Initialize log file
if ~isempty(options.logFile)
    logFid = fopen(options.logFile, 'w');
    writeLog(logFid, 'Batch processing started: %s', datestr(now));
    writeLog(logFid, 'Number of files: %d', nFiles);
else
    logFid = -1;
end

% Check if analysis function accepts index arguments
try
    nArgIn = nargin(analysisFcn);
catch
    nArgIn = 1;  % Assume single argument if nargin fails
end
useIndexArgs = nArgIn >= 3 || nArgIn < 0;  % -1 means varargin

% Initialize results structure
results(nFiles) = struct('fileName', '', 'filePath', '', 'success', false, ...
    'error', '', 'data', [], 'duration', 0);

% Process files
if options.parallel && ~isempty(ver('parallel'))
    % Parallel processing
    results = processParallel(dataFiles, analysisFcn, options, useIndexArgs, logFid);
else
    % Sequential processing
    results = processSequential(dataFiles, analysisFcn, options, useIndexArgs, logFid);
end

% Summary
nSuccess = sum([results.success]);
nFailed = nFiles - nSuccess;
totalTime = sum([results.duration]);

if options.showProgress
    fprintf('\n========================================\n');
    fprintf('Batch Processing Complete\n');
    fprintf('========================================\n');
    fprintf('Total files:   %d\n', nFiles);
    fprintf('Successful:    %d\n', nSuccess);
    fprintf('Failed:        %d\n', nFailed);
    fprintf('Total time:    %.1f s\n', totalTime);
    fprintf('========================================\n');
end

if logFid > 0
    writeLog(logFid, '');
    writeLog(logFid, 'Batch processing complete');
    writeLog(logFid, 'Successful: %d / %d', nSuccess, nFiles);
    writeLog(logFid, 'Total time: %.1f s', totalTime);
    fclose(logFid);
end

% Save final combined results
if options.saveResults
    outputFile = fullfile(options.outputFolder, [options.outputPrefix 'results.mat']);
    save(outputFile, 'results');
    fprintf('Results saved to: %s\n', outputFile);
end

end

%% Sequential Processing
function results = processSequential(dataFiles, analysisFcn, options, useIndexArgs, logFid)
    nFiles = length(dataFiles);

    % Initialize results
    results(nFiles) = struct('fileName', '', 'filePath', '', 'success', false, ...
        'error', '', 'data', [], 'duration', 0);

    if options.showProgress
        prog = ProgressTracker(nFiles, 'Processing files');
    end

    for i = 1:nFiles
        filePath = dataFiles{i};
        [~, fileName, ext] = fileparts(filePath);

        results(i).fileName = [fileName, ext];
        results(i).filePath = filePath;

        if logFid > 0
            writeLog(logFid, 'Processing [%d/%d]: %s', i, nFiles, results(i).fileName);
        end

        startTime = tic;

        try
            % Call analysis function
            if useIndexArgs
                results(i).data = analysisFcn(filePath, i, nFiles);
            else
                results(i).data = analysisFcn(filePath);
            end
            results(i).success = true;

            if logFid > 0
                writeLog(logFid, '  SUCCESS (%.2f s)', toc(startTime));
            end
        catch ME
            results(i).success = false;
            results(i).error = sprintf('%s: %s', ME.identifier, ME.message);

            if logFid > 0
                writeLog(logFid, '  FAILED: %s', results(i).error);
            end

            if ~options.continueOnError
                error('batchProcess:analysisError', ...
                    'Analysis failed for %s: %s', filePath, ME.message);
            end
        end

        results(i).duration = toc(startTime);

        % Save intermediate result
        if options.saveResults && results(i).success
            outputFile = fullfile(options.outputFolder, ...
                sprintf('%s%s.mat', options.outputPrefix, fileName));
            resultData = results(i).data;
            save(outputFile, 'resultData');
        end

        if options.showProgress
            prog.update(i);
        end
    end

    if options.showProgress
        prog.finish();
    end
end

%% Parallel Processing
function results = processParallel(dataFiles, analysisFcn, options, useIndexArgs, logFid)
    nFiles = length(dataFiles);

    % Initialize results as cell array for parfor
    resultCell = cell(nFiles, 1);
    fileNames = cell(nFiles, 1);
    filePaths = dataFiles;

    % Extract file names
    for i = 1:nFiles
        [~, name, ext] = fileparts(dataFiles{i});
        fileNames{i} = [name, ext];
    end

    if options.showProgress
        fprintf('Processing %d files in parallel...\n', nFiles);
    end

    % Create local copies for parfor
    continueOnError = options.continueOnError;
    saveResults = options.saveResults;
    outputFolder = options.outputFolder;
    outputPrefix = options.outputPrefix;

    parfor i = 1:nFiles
        result = struct('fileName', fileNames{i}, 'filePath', filePaths{i}, ...
            'success', false, 'error', '', 'data', [], 'duration', 0);

        startTime = tic;

        try
            if useIndexArgs
                result.data = analysisFcn(filePaths{i}, i, nFiles);
            else
                result.data = analysisFcn(filePaths{i});
            end
            result.success = true;
        catch ME
            result.success = false;
            result.error = sprintf('%s: %s', ME.identifier, ME.message);

            if ~continueOnError
                error('batchProcess:analysisError', ...
                    'Analysis failed for %s: %s', filePaths{i}, ME.message);
            end
        end

        result.duration = toc(startTime);

        % Save intermediate result
        if saveResults && result.success
            [~, name, ~] = fileparts(filePaths{i});
            outputFile = fullfile(outputFolder, sprintf('%s%s.mat', outputPrefix, name));
            resultData = result.data;
            parsave(outputFile, resultData);
        end

        resultCell{i} = result;
    end

    % Convert cell array back to structure array
    results = [resultCell{:}];

    % Log results (can't do this inside parfor)
    if logFid > 0
        for i = 1:nFiles
            if results(i).success
                writeLog(logFid, '[%d/%d] %s: SUCCESS (%.2f s)', ...
                    i, nFiles, results(i).fileName, results(i).duration);
            else
                writeLog(logFid, '[%d/%d] %s: FAILED - %s', ...
                    i, nFiles, results(i).fileName, results(i).error);
            end
        end
    end
end

%% Helper Functions
function writeLog(fid, format, varargin)
    if fid > 0
        timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        fprintf(fid, '[%s] %s\n', timestamp, sprintf(format, varargin{:}));
    end
end

function parsave(filename, data)
    % Save function that works inside parfor
    save(filename, 'data');
end
