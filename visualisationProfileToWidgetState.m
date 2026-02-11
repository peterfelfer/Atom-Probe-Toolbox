function state = visualisationProfileToWidgetState(profileOrState)
% VISUALISATIONPROFILETOWIDGETSTATE Convert profile/state structs to widget state.
%
% state = visualisationProfileToWidgetState(profileOrState)
%
% Accepts either:
% - canonical visualisation profiles (schema: visualisationProfile), or
% - legacy scatterPlotPosWidget state structs.

arguments
    profileOrState (1,1) struct
end

if isfield(profileOrState, 'settings') && isfield(profileOrState, 'species')
    profile = visualisationProfileMigrate(profileOrState);
    visualisationProfileValidate(profile);
    state = stateFromProfile(profile);
else
    state = profileOrState;
    if ~isfield(state, 'version') || isempty(state.version)
        state.version = 3;
    end
end

end

function state = stateFromProfile(profile)
settings = profile.settings;
species = profile.species;

state = struct();
state.version = 3;

state.mode = string(settings.mode);
state.splitIsotope = logical(settings.splitIsotope);
state.splitCharge = logical(settings.splitCharge);
state.showUnranged = logical(settings.showUnranged);
state.sampleValue = double(settings.sampleValue);
state.markerSize = double(settings.markerSize);
state.markerSizeGlobal = double(settings.markerSizeGlobal);
state.randomSeed = double(settings.randomSeed);
state.rotationAngle = double(settings.rotationAngle);
state.fixBoundingBox = logical(settings.fixBoundingBox);
state.clipToBoundingBox = logical(settings.clipToBoundingBox);
state.bboxUseSpan = logical(settings.bboxUseSpan);
state.showAxisLabels = logical(settings.showAxisLabels);
state.showAxisTicks = logical(settings.showAxisTicks);
state.showScaleCube = logical(settings.showScaleCube);
state.scaleCubeSize = double(settings.scaleCubeSize);
state.searchFilter = string(settings.searchFilter);
state.sortMode = string(settings.sortMode);

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
