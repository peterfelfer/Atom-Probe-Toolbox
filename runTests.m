function results = runTests(options)
% RUNTESTS Run all tests for the Atom Probe Toolbox
%
% results = runTests()
% results = runTests('pattern', 'test_*.m')
% results = runTests('verbose', true)
%
% Discovers and runs test files throughout the toolbox. Test files should
% be named with 'test_' prefix and contain functions that return true/false.
%
% OPTIONS:
%   'pattern'    - File pattern to match (default: 'test_*.m')
%   'folder'     - Specific folder to search (default: all)
%   'verbose'    - Show detailed output (default: true)
%   'stopOnFail' - Stop on first failure (default: false)
%
% OUTPUT:
%   results - Structure with test results:
%       .total    - Total tests run
%       .passed   - Number passed
%       .failed   - Number failed
%       .skipped  - Number skipped
%       .duration - Total time
%       .details  - Table with per-test results
%
% TEST FILE FORMAT:
%   Test files should define functions that return true (pass) or false (fail).
%   Functions starting with 'test' are automatically discovered.
%
%   Example test_myFeature.m:
%   ```
%   function success = test_basic()
%       x = myFunction(1, 2);
%       success = (x == 3);
%   end
%
%   function success = test_edge_case()
%       try
%           myFunction([], []);
%           success = false;  % Should have thrown error
%       catch
%           success = true;
%       end
%   end
%   ```
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    options.pattern (1,:) char = 'test_*.m'
    options.folder (1,:) char = ''
    options.verbose (1,1) logical = true
    options.stopOnFail (1,1) logical = false
end

startTime = tic;

toolboxRoot = fileparts(mfilename('fullpath'));

% Find test files
if isempty(options.folder)
    searchPath = toolboxRoot;
else
    searchPath = fullfile(toolboxRoot, options.folder);
end

testFiles = dir(fullfile(searchPath, '**', options.pattern));

if isempty(testFiles)
    fprintf('No test files found matching pattern: %s\n', options.pattern);
    results = struct('total', 0, 'passed', 0, 'failed', 0, 'skipped', 0);
    return;
end

if options.verbose
    fprintf('\n');
    fprintf('==========================================\n');
    fprintf('  Atom Probe Toolbox Test Runner\n');
    fprintf('==========================================\n');
    fprintf('  Found %d test files\n', length(testFiles));
    fprintf('==========================================\n\n');
end

% Initialize results
testResults = {};
totalTests = 0;
passedTests = 0;
failedTests = 0;
skippedTests = 0;

% Run each test file
for f = 1:length(testFiles)
    testFile = testFiles(f);
    testPath = fullfile(testFile.folder, testFile.name);
    [~, testName, ~] = fileparts(testFile.name);

    if options.verbose
        fprintf('Running: %s\n', testName);
    end

    % Get test functions from file
    testFunctions = getTestFunctions(testPath);

    if isempty(testFunctions)
        if options.verbose
            fprintf('  (no test functions found)\n');
        end
        continue;
    end

    % Run each test function
    for t = 1:length(testFunctions)
        funcName = testFunctions{t};
        totalTests = totalTests + 1;

        result = struct();
        result.file = testName;
        result.function = funcName;
        result.passed = false;
        result.skipped = false;
        result.error = '';
        result.duration = 0;

        funcStart = tic;

        try
            % Navigate to test directory and run
            oldDir = cd(testFile.folder);
            cleanup = onCleanup(@() cd(oldDir));

            % Call the test function
            funcHandle = str2func(funcName);
            testPassed = funcHandle();

            result.duration = toc(funcStart);

            if islogical(testPassed) && testPassed
                result.passed = true;
                passedTests = passedTests + 1;
                statusStr = 'PASS';
            elseif islogical(testPassed) && ~testPassed
                result.passed = false;
                failedTests = failedTests + 1;
                statusStr = 'FAIL';
            else
                % Non-boolean return - treat as skip
                result.skipped = true;
                skippedTests = skippedTests + 1;
                statusStr = 'SKIP';
            end

        catch ME
            result.duration = toc(funcStart);
            result.passed = false;
            result.error = sprintf('%s: %s', ME.identifier, ME.message);
            failedTests = failedTests + 1;
            statusStr = 'ERROR';
        end

        if options.verbose
            if result.passed
                fprintf('  [%s] %s (%.3fs)\n', statusStr, funcName, result.duration);
            else
                fprintf('  [%s] %s (%.3fs)\n', statusStr, funcName, result.duration);
                if ~isempty(result.error)
                    fprintf('        %s\n', result.error);
                end
            end
        end

        testResults{end+1} = result;

        if options.stopOnFail && ~result.passed && ~result.skipped
            fprintf('\nStopping on first failure.\n');
            break;
        end
    end

    if options.stopOnFail && failedTests > 0
        break;
    end

    if options.verbose
        fprintf('\n');
    end
end

totalDuration = toc(startTime);

% Build results structure
results = struct();
results.total = totalTests;
results.passed = passedTests;
results.failed = failedTests;
results.skipped = skippedTests;
results.duration = totalDuration;

if ~isempty(testResults)
    results.details = struct2table([testResults{:}]);
end

% Print summary
if options.verbose
    fprintf('==========================================\n');
    fprintf('  Test Summary\n');
    fprintf('==========================================\n');
    fprintf('  Total:   %d\n', totalTests);
    fprintf('  Passed:  %d\n', passedTests);
    fprintf('  Failed:  %d\n', failedTests);
    fprintf('  Skipped: %d\n', skippedTests);
    fprintf('  Time:    %.2f s\n', totalDuration);
    fprintf('==========================================\n');

    if failedTests == 0
        fprintf('  All tests passed!\n');
    else
        fprintf('  %d test(s) failed.\n', failedTests);
    end
    fprintf('==========================================\n\n');
end

end

function testFunctions = getTestFunctions(testFilePath)
    % Extract test function names from a test file

    testFunctions = {};

    try
        % Read file content
        fid = fopen(testFilePath, 'r');
        if fid == -1
            return;
        end
        content = fread(fid, '*char')';
        fclose(fid);

        % Find function definitions starting with 'test'
        pattern = 'function\s+\w+\s*=\s*(test\w*)\s*\(';
        matches = regexp(content, pattern, 'tokens');

        for i = 1:length(matches)
            testFunctions{end+1} = matches{i}{1};
        end

        % Also check for functions without return value
        pattern2 = 'function\s+(test\w*)\s*\(';
        matches2 = regexp(content, pattern2, 'tokens');

        for i = 1:length(matches2)
            if ~any(strcmp(testFunctions, matches2{i}{1}))
                testFunctions{end+1} = matches2{i}{1};
            end
        end

    catch
        % Failed to parse file
    end
end
