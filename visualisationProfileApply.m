function [scatterHandles, ax, controlFig, info, resolvedProfile] = visualisationProfileApply( ...
    pos, colorScheme, visualisationProfile, options)
% VISUALISATIONPROFILEAPPLY Apply a visualisation profile to a dataset.
%
% [scatterHandles, ax, controlFig, info, resolvedProfile] = ...
%    visualisationProfileApply(pos, colorScheme, visualisationProfile)

arguments
    pos table
    colorScheme
    visualisationProfile (1,1) struct
    options.axes = []
    options.showWidget (1,1) logical = false
    options.controlTitle (1,1) string = "APT Scatter Controls"
    options.speciesPolicy (1,1) string = ""
    options.restoreStateFromAxis (1,1) logical = false
    options.persistStateToAxis (1,1) logical = true
end

[resolvedProfile, resolveInfo] = visualisationProfileResolve( ...
    pos, colorScheme, visualisationProfile, 'speciesPolicy', options.speciesPolicy);
settings = resolvedProfile.settings;

groupBy = "auto";
mode = lower(string(settings.mode));
if mode == "ionic"
    groupBy = "ion";
elseif mode == "atomic"
    groupBy = "atom";
end

[scatterHandles, ax, controlFig, widgetInfo] = scatterPlotPosWidget(pos, colorScheme, ...
    'axes', options.axes, ...
    'groupBy', groupBy, ...
    'sample', double(settings.sampleValue), ...
    'markerSizeGlobal', double(settings.markerSizeGlobal), ...
    'randomSeed', double(settings.randomSeed), ...
    'showUnranged', logical(settings.showUnranged), ...
    'controlTitle', options.controlTitle, ...
    'splitIsotope', logical(settings.splitIsotope), ...
    'splitCharge', logical(settings.splitCharge), ...
    'fixBoundingBox', logical(settings.fixBoundingBox), ...
    'clipToBoundingBox', logical(settings.clipToBoundingBox), ...
    'bboxUseSpan', logical(settings.bboxUseSpan), ...
    'restoreStateFromAxis', options.restoreStateFromAxis, ...
    'persistStateToAxis', options.persistStateToAxis);

state = stateFromResolvedProfile(resolvedProfile);
scatterPlotPosWidgetApplyState(controlFig, state);

if ~options.showWidget && isgraphics(controlFig)
    controlFig.Visible = 'off';
end

info = struct();
info.resolve = resolveInfo;
info.widget = widgetInfo;

end

function state = stateFromResolvedProfile(profile)
state = struct();
state.version = 3;

settings = profile.settings;
species = profile.species;

state.mode = settings.mode;
state.splitIsotope = settings.splitIsotope;
state.splitCharge = settings.splitCharge;
state.showUnranged = settings.showUnranged;
state.sampleValue = settings.sampleValue;
state.markerSize = settings.markerSize;
state.markerSizeGlobal = settings.markerSizeGlobal;
state.randomSeed = settings.randomSeed;
state.rotationAngle = settings.rotationAngle;
state.fixBoundingBox = settings.fixBoundingBox;
state.clipToBoundingBox = settings.clipToBoundingBox;
state.bboxUseSpan = settings.bboxUseSpan;
state.showAxisLabels = settings.showAxisLabels;
state.showAxisTicks = settings.showAxisTicks;
state.showScaleCube = settings.showScaleCube;
state.scaleCubeSize = settings.scaleCubeSize;
state.searchFilter = settings.searchFilter;
state.sortMode = settings.sortMode;

if isfield(settings, 'boundingBox')
    state.boundingBox = settings.boundingBox;
end
if isfield(settings, 'camera')
    state.camera = settings.camera;
end

state.species = struct();
state.species.name = string(species.name(:));
state.species.visible = logical(species.visible(:));
state.species.fraction = double(species.fraction(:));
state.species.markerSize = double(species.markerSize(:));
state.species.markerSizeLinked = logical(species.markerSizeLinked(:));
state.species.color = double(species.color);

if isfield(profile, 'meta') && isfield(profile.meta, 'posSignature')
    state.posSignature = profile.meta.posSignature;
end

end
