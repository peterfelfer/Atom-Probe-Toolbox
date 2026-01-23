function setupToolbox(options)
% SETUPTOOLBOX Initialize the Atom Probe Toolbox
%
% setupToolbox()
% setupToolbox('checkDeps', true)
% setupToolbox('quiet', true)
%
% Adds all necessary paths and optionally checks dependencies.
% Run this once per MATLAB session before using the toolbox.
%
% OPTIONS:
%   'checkDeps'  - Run dependency check (default: true)
%   'quiet'      - Suppress output messages (default: false)
%   'permanent'  - Save paths permanently (default: false)
%
% EXAMPLES:
%   setupToolbox();                    % Standard setup
%   setupToolbox('quiet', true);       % Silent setup
%   setupToolbox('permanent', true);   % Add to MATLAB path permanently
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    options.checkDeps (1,1) logical = true
    options.quiet (1,1) logical = false
    options.permanent (1,1) logical = false
end

% Get toolbox root directory
toolboxRoot = fileparts(mfilename('fullpath'));

if ~options.quiet
    fprintf('\n');
    fprintf('==========================================\n');
    fprintf('  Atom Probe Toolbox Setup\n');
    fprintf('==========================================\n');
    fprintf('  Version: 1.0\n');
    fprintf('  Location: %s\n', toolboxRoot);
    fprintf('==========================================\n\n');
end

% Define directories to add
subfolders = {
    '',                     % Root
    'analysis'
    'correlative'
    'crystallography'
    'utilities'
    'utilities_Cluster'
    'utilities_dualMesh'
    'utilities_geom2d'
    'utilities_geom2d/geom2d'
    'utilities_geom2d/polygons2d'
    'utilities_IO'
    'utilities_IO/APT_reader_MPIE'
    'utilities_IO/EMIODist2'
    'doc'
};

% Add paths
nAdded = 0;
for i = 1:length(subfolders)
    folderPath = fullfile(toolboxRoot, subfolders{i});
    if isfolder(folderPath)
        if ~contains(path, folderPath)
            addpath(folderPath);
            nAdded = nAdded + 1;
        end
    else
        if ~options.quiet
            warning('setupToolbox:missingFolder', 'Folder not found: %s', subfolders{i});
        end
    end
end

if ~options.quiet
    fprintf('Added %d folders to MATLAB path.\n\n', nAdded);
end

% Save path permanently if requested
if options.permanent
    try
        savepath;
        if ~options.quiet
            fprintf('Path saved permanently.\n\n');
        end
    catch ME
        warning('setupToolbox:savePathFailed', ...
            'Could not save path permanently: %s', ME.message);
    end
end

% Check dependencies
if options.checkDeps
    if ~options.quiet
        checkDependencies('verbose', true);
    else
        [status, ~] = checkDependencies('verbose', false);
        if ~status
            warning('setupToolbox:missingDeps', ...
                'Some dependencies are missing. Run checkDependencies() for details.');
        end
    end
end

if ~options.quiet
    fprintf('\nAtom Probe Toolbox is ready.\n');
    fprintf('Run "doc GettingStarted" for documentation.\n\n');
end

end
