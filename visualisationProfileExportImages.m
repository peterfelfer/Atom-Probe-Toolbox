function exportInfo = visualisationProfileExportImages(pos, colorScheme, visualisationProfile, options)
% VISUALISATIONPROFILEEXPORTIMAGES Scriptable image export for profile-based visualisation.
%
% exportInfo = visualisationProfileExportImages(pos, colorScheme, visualisationProfile, ...)

arguments
    pos table
    colorScheme
    visualisationProfile (1,1) struct
    options.outputDir (1,:) char = pwd
    options.baseName (1,:) char = 'apt_visualisation'
    options.extension (1,:) char = ''
    options.perSpecies (1,1) logical = true
    options.includeCombined (1,1) logical = true
    options.axes = []
    options.speciesPolicy (1,1) string = ""
    options.showWidget (1,1) logical = false
    options.closeControlFig (1,1) logical = true
end

if ~isfolder(options.outputDir)
    mkdir(options.outputDir);
end

profile = visualisationProfileMigrate(visualisationProfile);
visualisationProfileValidate(profile);

ext = string(options.extension);
if strlength(ext) == 0
    ext = string(profile.export.imageExtension);
end
ext = lower(ext);
if ~startsWith(ext, '.')
    ext = '.' + ext;
end

[~, ax, controlFig, infoApply, resolvedProfile] = visualisationProfileApply( ...
    pos, colorScheme, profile, ...
    'axes', options.axes, ...
    'showWidget', options.showWidget, ...
    'speciesPolicy', options.speciesPolicy, ...
    'restoreStateFromAxis', false, ...
    'persistStateToAxis', true);

data = getappdata(controlFig, 'scatterPlotPosWidget');
if isempty(data) || ~isfield(data, 'fn')
    error('visualisationProfileExportImages:missingWidgetState', ...
        'Could not access scatterPlotPosWidget state for export.');
end

exportInfo = struct();
exportInfo.files = strings(0, 1);
exportInfo.species = strings(0, 1);
exportInfo.resolvedProfile = resolvedProfile;
exportInfo.applyInfo = infoApply;

isVector = ismember(ext, [".pdf", ".eps", ".svg"]);
baseName = string(options.baseName);

if options.includeCombined
    outFile = fullfile(options.outputDir, char(baseName + ext));
    exportOne(ax, outFile, isVector);
    exportInfo.files(end+1, 1) = string(outFile); %#ok<AGROW>
    exportInfo.species(end+1, 1) = "combined"; %#ok<AGROW>
end

if options.perSpecies
    visibleIdx = find(data.visible);
    originalVisible = data.visible;

    for k = 1:numel(visibleIdx)
        i = visibleIdx(k);
        data.visible(:) = false;
        data.visible(i) = true;
        setappdata(controlFig, 'scatterPlotPosWidget', data);
        data.fn.updateScatter(controlFig);
        drawnow;

        tag = sanitizeFileTagLocal(data.displayNames(i));
        outFile = fullfile(options.outputDir, char(baseName + "_" + tag + ext));
        exportOne(ax, outFile, isVector);

        exportInfo.files(end+1, 1) = string(outFile); %#ok<AGROW>
        exportInfo.species(end+1, 1) = string(data.displayNames(i)); %#ok<AGROW>
    end

    data.visible = originalVisible;
    setappdata(controlFig, 'scatterPlotPosWidget', data);
    data.fn.updateScatter(controlFig);
end

if options.closeControlFig && isgraphics(controlFig)
    close(controlFig);
end

end

function exportOne(ax, outFile, isVector)
if isVector
    exportgraphics(ax, outFile, 'ContentType', 'vector');
else
    exportgraphics(ax, outFile, 'Resolution', 300);
end
end

function tag = sanitizeFileTagLocal(nameIn)
tag = string(nameIn);
if strlength(tag) == 0
    tag = "species";
end
tag = regexprep(tag, '[^a-zA-Z0-9]+', '_');
tag = regexprep(tag, '^_+|_+$', '');
if strlength(tag) == 0
    tag = "species";
end
end
