function [status, report] = checkDependencies(options)
% CHECKDEPENDENCIES Verify MATLAB version and required toolboxes
%
% [status, report] = checkDependencies()
% [status, report] = checkDependencies('verbose', true)
% [status, report] = checkDependencies('autoAddPath', true)
%
% Checks that all required MATLAB toolboxes and internal utilities are
% available for the Atom Probe Toolbox to function correctly.
%
% OPTIONS:
%   'verbose'     - Display detailed information (default: true)
%   'autoAddPath' - Automatically add toolbox paths (default: true)
%   'checkOptional' - Also check optional dependencies (default: false)
%
% OUTPUT:
%   status - true if all required dependencies are met
%   report - Structure containing detailed dependency information
%       .matlabVersion  - Current MATLAB version
%       .required       - Table of required dependencies and their status
%       .optional       - Table of optional dependencies and their status
%       .missing        - Cell array of missing required dependencies
%       .warnings       - Cell array of warning messages
%
% EXAMPLE:
%   checkDependencies();  % Quick check with console output
%   [ok, rep] = checkDependencies('verbose', false);  % Silent check
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    options.verbose (1,1) logical = true
    options.autoAddPath (1,1) logical = true
    options.checkOptional (1,1) logical = false
end

report = struct();
report.matlabVersion = version;
report.matlabRelease = version('-release');
report.warnings = {};
report.missing = {};

% Minimum MATLAB version (R2019b for arguments block support)
minVersion = '9.7';  % R2019b
currentVersion = version;
versionNum = regexp(currentVersion, '^\d+\.\d+', 'match', 'once');

if str2double(versionNum) < str2double(minVersion)
    report.warnings{end+1} = sprintf('MATLAB version %s or later recommended (current: %s)', ...
        minVersion, currentVersion);
end

%% Check Required Toolboxes
requiredToolboxes = {
    'MATLAB',                           'Core MATLAB',          true
    'Statistics and Machine Learning Toolbox', 'Statistics',    true
    'Image Processing Toolbox',         'Image Processing',     true
};

toolboxStatus = cell(size(requiredToolboxes, 1), 4);
v = ver;
installedToolboxes = {v.Name};

for i = 1:size(requiredToolboxes, 1)
    tbName = requiredToolboxes{i, 1};
    tbShort = requiredToolboxes{i, 2};
    tbRequired = requiredToolboxes{i, 3};

    isInstalled = any(strcmpi(installedToolboxes, tbName));
    toolboxStatus{i, 1} = tbShort;
    toolboxStatus{i, 2} = tbName;
    toolboxStatus{i, 3} = isInstalled;
    toolboxStatus{i, 4} = tbRequired;

    if ~isInstalled && tbRequired
        report.missing{end+1} = tbName;
    end
end

report.required = cell2table(toolboxStatus, ...
    'VariableNames', {'ShortName', 'FullName', 'Installed', 'Required'});

%% Check Optional Toolboxes
optionalToolboxes = {
    'Curve Fitting Toolbox',            'Curve Fitting',        false
    'Optimization Toolbox',             'Optimization',         false
    'Parallel Computing Toolbox',       'Parallel Computing',   false
    'Computer Vision Toolbox',          'Computer Vision',      false
};

optionalStatus = cell(size(optionalToolboxes, 1), 4);

for i = 1:size(optionalToolboxes, 1)
    tbName = optionalToolboxes{i, 1};
    tbShort = optionalToolboxes{i, 2};

    isInstalled = any(strcmpi(installedToolboxes, tbName));
    optionalStatus{i, 1} = tbShort;
    optionalStatus{i, 2} = tbName;
    optionalStatus{i, 3} = isInstalled;
    optionalStatus{i, 4} = false;  % Not required
end

report.optional = cell2table(optionalStatus, ...
    'VariableNames', {'ShortName', 'FullName', 'Installed', 'Required'});

%% Check Internal Utilities
toolboxRoot = fileparts(mfilename('fullpath'));
internalPaths = {
    'utilities',            'Core utilities'
    'utilities_Cluster',    'Cluster analysis'
    'utilities_dualMesh',   'Dual mesh generation'
    'utilities_geom2d',     'Geometry utilities'
    'utilities_IO',         'Input/Output'
    'analysis',             'Analysis modules'
    'correlative',          'Correlative imaging'
    'crystallography',      'Crystallography'
    'doc',                  'Documentation'
};

internalStatus = cell(size(internalPaths, 1), 3);

for i = 1:size(internalPaths, 1)
    pathName = internalPaths{i, 1};
    pathDesc = internalPaths{i, 2};
    fullPath = fullfile(toolboxRoot, pathName);

    exists = isfolder(fullPath);
    internalStatus{i, 1} = pathName;
    internalStatus{i, 2} = pathDesc;
    internalStatus{i, 3} = exists;

    if ~exists
        report.warnings{end+1} = sprintf('Internal folder missing: %s', pathName);
    elseif options.autoAddPath
        % Add to path if not already there
        if ~contains(path, fullPath)
            addpath(genpath(fullPath));
        end
    end
end

report.internal = cell2table(internalStatus, ...
    'VariableNames', {'Folder', 'Description', 'Exists'});

%% Check for GPU availability (informational)
report.gpuAvailable = false;
report.gpuInfo = '';
try
    if exist('gpuDeviceCount', 'file') && gpuDeviceCount > 0
        gpu = gpuDevice;
        report.gpuAvailable = true;
        report.gpuInfo = sprintf('%s (%.1f GB)', gpu.Name, gpu.TotalMemory/1e9);
    end
catch
    % GPU not available or error checking
end

%% Check parallel pool availability
report.parallelAvailable = ~isempty(ver('parallel'));
report.parallelWorkers = 0;
if report.parallelAvailable
    try
        pool = gcp('nocreate');
        if ~isempty(pool)
            report.parallelWorkers = pool.NumWorkers;
        end
    catch
        % Parallel pool not available
    end
end

%% Determine overall status
status = isempty(report.missing);

%% Display results if verbose
if options.verbose
    fprintf('\n');
    fprintf('========================================\n');
    fprintf('  Atom Probe Toolbox Dependency Check\n');
    fprintf('========================================\n\n');

    fprintf('MATLAB Version: %s (%s)\n', report.matlabVersion, report.matlabRelease);
    fprintf('\n');

    % Required toolboxes
    fprintf('REQUIRED TOOLBOXES:\n');
    fprintf('-------------------\n');
    for i = 1:height(report.required)
        if report.required.Installed(i)
            statusStr = '[OK]';
        else
            statusStr = '[MISSING]';
        end
        fprintf('  %s %s\n', statusStr, report.required.ShortName{i});
    end
    fprintf('\n');

    % Optional toolboxes
    if options.checkOptional
        fprintf('OPTIONAL TOOLBOXES:\n');
        fprintf('-------------------\n');
        for i = 1:height(report.optional)
            if report.optional.Installed(i)
                statusStr = '[OK]';
            else
                statusStr = '[--]';
            end
            fprintf('  %s %s\n', statusStr, report.optional.ShortName{i});
        end
        fprintf('\n');
    end

    % Internal modules
    fprintf('INTERNAL MODULES:\n');
    fprintf('-----------------\n');
    for i = 1:height(report.internal)
        if report.internal.Exists(i)
            statusStr = '[OK]';
        else
            statusStr = '[MISSING]';
        end
        fprintf('  %s %s\n', statusStr, report.internal.Folder{i});
    end
    fprintf('\n');

    % Hardware capabilities
    fprintf('HARDWARE CAPABILITIES:\n');
    fprintf('----------------------\n');
    if report.gpuAvailable
        fprintf('  [OK] GPU: %s\n', report.gpuInfo);
    else
        fprintf('  [--] GPU: Not available (CPU fallback will be used)\n');
    end
    if report.parallelAvailable
        if report.parallelWorkers > 0
            fprintf('  [OK] Parallel: %d workers active\n', report.parallelWorkers);
        else
            fprintf('  [OK] Parallel: Available (no active pool)\n');
        end
    else
        fprintf('  [--] Parallel: Not available\n');
    end
    fprintf('\n');

    % Warnings
    if ~isempty(report.warnings)
        fprintf('WARNINGS:\n');
        fprintf('---------\n');
        for i = 1:length(report.warnings)
            fprintf('  ! %s\n', report.warnings{i});
        end
        fprintf('\n');
    end

    % Summary
    fprintf('========================================\n');
    if status
        fprintf('  Status: READY\n');
        fprintf('  All required dependencies are met.\n');
    else
        fprintf('  Status: MISSING DEPENDENCIES\n');
        fprintf('  Please install: %s\n', strjoin(report.missing, ', '));
    end
    fprintf('========================================\n\n');
end

end
