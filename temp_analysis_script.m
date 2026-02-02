%% Temporary Analysis Script for A718 EPOS data
% This script loads and analyzes the EPOS file

% Add toolbox to path
addpath(genpath(pwd));

% File path
eposFile = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/recons/recon-v01/default/R56_09162-v01.epos';

% Create output folder
[parentDir, ~, ~] = fileparts(eposFile);
outputDir = fullfile(parentDir, 'analysis_output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
fprintf('Output directory: %s\n', outputDir);

% Load the EPOS file
fprintf('Loading EPOS file...\n');
pos = posLoad(eposFile);
fprintf('Loaded %d atoms\n', height(pos));

% Show the columns available
fprintf('Columns in pos table:\n');
disp(pos.Properties.VariableNames);

% Get mass range
fprintf('Mass range: %.2f - %.2f Da\n', min(pos.mc), max(pos.mc));

% Save pos to output folder
save(fullfile(outputDir, 'pos.mat'), 'pos');
fprintf('Saved pos.mat\n');

fprintf('Step 1 complete!\n');
