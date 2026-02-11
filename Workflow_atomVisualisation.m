%[text] # **Atom Visualisation Workflow**
%[text] **Interactive and scriptable atom visualisation with** **`scatterPlotPosWidget`** **and visualisation profiles**
%[text] This workflow is designed as a practical bridge between interactive exploration and reproducible automation. In a typical APT analysis session, users first inspect the reconstruction visually, tune what is shown (species visibility, fractions, marker sizes, camera view, clipping), and then need to apply exactly the same visualisation settings to other datasets or regenerate figures later in a consistent way. This script demonstrates that full cycle with a clear progression from minimal function calls to profile-based scripting.
%[text] 
%[text] The workflow is intentionally structured in blocks that correspond to common analysis tasks. First, a position dataset is loaded and linked to a ranged mass spectrum figure. The ranges are extracted and allocated so the plotting table carries ion identity information. Next, `scatterPlotPosWidget` is launched with the simplest call so beginners can start quickly without a large options list. After interactive tuning, the current visualisation is captured as a profile struct and exported to multiple formats. That profile is then imported and applied in non-interactive mode to show how visualisation can become fully scriptable. Finally, still-image and turntable exports are created from the same profile, and figure state persistence is demonstrated by saving and re-opening a scatter figure.
%[text] 
%[text] Conceptually, this script uses four central function groups, and each section maps directly to one of them.
%[text] If you run this file section-by-section, you can stop after the interactive steps for exploratory work, or continue through the profile/export sections for fully reproducible workflows.
%[text] 
%[text] The toolbox is assumed to already be on the MATLAB path.
%%
%[text] ## 1) Load input data and ranged mass spectrum
%[text] Beginner version: one `posLoad(...)` call and one `openfig(...)` call.
%[text] If the default example files are not available, file dialogs are used as fallback.

% Root folder of this workflow file (used for outputs and documentation assets).
toolboxRoot = fileparts(mfilename('fullpath'));

% Folder where workflow outputs are written.
outputDir = fullfile(toolboxRoot, 'CodexTest_visualisation');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

% Example input files used for demonstration when present.
exampleEposFile = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/recons/recon-v01/default/R56_09162-v01.epos';
exampleMassSpecFigFile = '/Users/peterfelfer/Dropbox/research/01_research_Erlangen/05_publications/05_01_papers_unfinished_to_write/2025_OttMonGobFel_LaserAPTOxcart/Göbel25_MA_Hydrogen Laser-Oxcart/Atom Probe Data/LEAP/A718_9162/Massspectrum.fig';

posRaw = posLoad(exampleEposFile);
specFig = openfig(exampleMassSpecFigFile);

% Resolve the mass-spectrum graphics handle from the opened figure.
spec = resolveMassSpecHandle(specFig);

fprintf('Loaded raw ions: %d\n', height(posRaw));
%%
%[text] ## 2) Load colour scheme, extract ranges, and allocate ions
%[text] This converts raw positions into a plotting-ready `pos` table that includes ion/range assignments.
% Load default colour scheme used by plotting functions.
load colorScheme.mat

% Pull all defined ranges from the ranged mass spectrum figure.
rangeTable = rangesExtractFromMassSpec(spec);

% Allocate each ion hit to its corresponding range; 'raw' keeps one row per ion event.
pos = posAllocateRange(posRaw, rangeTable, 'raw');
%%
%[text] ## 3) Launch interactive scatter widget (minimal call)
%[text] Start with the default call so beginners can focus on interaction first.
%[text] You can tune visibility, fraction shown, marker sizes, and view directly in the widget.
% Minimal call with defaults.
[scatterHandles, ax, controlFig, infoInteractive] = scatterPlotPosWidget(pos, colorScheme);
%%
%[text] ## 
%%
%[text] ## 4) Capture and export a visualisation profile
%[text] The profile captures current widget state and enables reproducible re-use.
%[text] It is exported as `MAT` (canonical), `JSON`, and `YAML`.

% Capture profile from current interactive widget state.
visualisationProfile = visualisationProfileFromWidget(controlFig);

% Validate schema before writing to disk.
visualisationProfileValidate(visualisationProfile);

% Target files for the three config formats.
profileMatFile = fullfile(outputDir, 'visualisationProfile_example.mat');
profileJsonFile = fullfile(outputDir, 'visualisationProfile_example.json');
profileYamlFile = fullfile(outputDir, 'visualisationProfile_example.yaml');

% Export profile in all supported formats.
save(profileMatFile, 'visualisationProfile');
configExport(visualisationProfile, profileJsonFile);
configExport(visualisationProfile, profileYamlFile);

% Also expose profile in base workspace for ad-hoc scripting.
assignin('base', 'visualisationProfile', visualisationProfile);

fprintf('Exported profile files:\n');
fprintf('  %s\n  %s\n  %s\n', profileMatFile, profileJsonFile, profileYamlFile);
%%
%[text] ## 5) Import profile and apply scripted visualisation (minimal call)
%[text] This shows non-interactive re-use of the profile with the default apply call.

% Import profile from YAML (human-readable format).
profileFromYaml = configImport(profileYamlFile);

% Minimal scripted apply call; returns resolved profile and diagnostic info.
[~, axScripted, controlFigScripted, infoScripted, resolvedProfile] = ...
    visualisationProfileApply(pos, colorScheme, profileFromYaml);

% Close control figure if created/hidden by apply routine.
if isgraphics(controlFigScripted)
    delete(controlFigScripted);
end

fprintf('Scripted apply figure: %s\n', char(get(ancestor(axScripted, 'figure'), 'Name')));
fprintf('Matched/new/missing species: %d / %d / %d\n', ...
    numel(infoScripted.resolve.matchedSpecies), ...
    numel(infoScripted.resolve.newSpecies), ...
    numel(infoScripted.resolve.missingSpecies));
%%
%[text] ## 6) Scripted exports: images and turntable
%[text] Minimal calls are used; only output locations are specified explicitly.

% Output directory for still-image exports.
imageOutDir = fullfile(outputDir, 'exports_images');
if ~isfolder(imageOutDir)
    mkdir(imageOutDir);
end

% Minimal image-export call (defaults handle base name, per-species, combined export, etc.).
imageInfo = visualisationProfileExportImages(pos, colorScheme, resolvedProfile, ...
    'outputDir', imageOutDir);

% Minimal turntable export call with explicit output file path.
turntableFile = fullfile(outputDir, 'visualisation_turntable.mp4');
turntableInfo = visualisationProfileExportTurntable(pos, colorScheme, resolvedProfile, ...
    'outputFile', turntableFile);

fprintf('Exported still-image files: %d\n', numel(imageInfo.files));
fprintf('Turntable file:\n  %s\n', turntableInfo.outputFile);
%%
%[text] ## 7) Save figure and resume widget later
%[text] This demonstrates that widget state can persist with the saved axes/figure.

% Save current interactive scatter figure.
savedScatterFig = fullfile(outputDir, 'scatter_saved.fig');
savefig(ancestor(ax, 'figure'), savedScatterFig);

% Reopen the saved scatter figure and find its plotting axes.
figReloaded = openfig(savedScatterFig);
axReloaded = findobj(figReloaded, 'Type', 'axes', '-not', 'Tag', 'legend');
axReloaded = axReloaded(1);

% Reattach widget to existing axes (keeps beginner-level call, one required option).
[~, ~, controlFigReloaded] = scatterPlotPosWidget(pos, colorScheme, 'axes', axReloaded);

%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":25}
%---
