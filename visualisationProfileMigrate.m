function visualisationProfile = visualisationProfileMigrate(profileIn)
% VISUALISATIONPROFILEMIGRATE Migrate profile/state structs to canonical schema.
%
% visualisationProfile = visualisationProfileMigrate(profileIn)

if nargin < 1 || isempty(profileIn)
    profileIn = struct();
end
if ~isstruct(profileIn) || ~isscalar(profileIn)
    error('visualisationProfileMigrate:invalidInput', ...
        'Input must be a scalar struct.');
end

if isfield(profileIn, 'settings') && isfield(profileIn, 'species')
    visualisationProfile = profileIn;
else
    visualisationProfile = convertLegacyState(profileIn);
end

visualisationProfile = applyDefaults(visualisationProfile);

end

function profile = convertLegacyState(state)
profile = struct();
profile.version = 1;
profile.schema = 'visualisationProfile';
profile.settings = struct();
profile.species = struct();

settingFields = {'mode', 'splitIsotope', 'splitCharge', 'showUnranged', ...
    'sampleValue', 'markerSize', 'markerSizeGlobal', 'randomSeed', ...
    'rotationAngle', 'boundingBox', 'fixBoundingBox', 'clipToBoundingBox', ...
    'bboxUseSpan', 'liveUpdate', 'showAxisLabels', 'showAxisTicks', ...
    'showScaleCube', 'scaleCubeSize', 'searchFilter', 'sortMode', ...
    'camera'};

for i = 1:numel(settingFields)
    fieldName = settingFields{i};
    if isfield(state, fieldName)
        profile.settings.(fieldName) = state.(fieldName);
    end
end

if isfield(state, 'bbox') && ~isfield(profile.settings, 'boundingBox')
    profile.settings.boundingBox = state.bbox;
end
if isfield(state, 'camera')
    profile.settings.camera = state.camera;
end

if isfield(state, 'species')
    profile.species = state.species;
else
    species = struct();
    if isfield(state, 'speciesNames')
        species.name = state.speciesNames;
    else
        species.name = strings(0, 1);
    end
    if isfield(state, 'visible')
        species.visible = state.visible;
    else
        species.visible = true(numel(species.name), 1);
    end
    if isfield(state, 'sampleFractions')
        species.fraction = state.sampleFractions;
    else
        species.fraction = ones(numel(species.name), 1);
    end
    if isfield(state, 'markerSizes')
        species.markerSize = state.markerSizes;
    end
    if isfield(state, 'markerSizeLinked')
        species.markerSizeLinked = state.markerSizeLinked;
    end
    if isfield(state, 'colors')
        species.color = state.colors;
    end
    profile.species = species;
end

if isfield(state, 'posSignature')
    profile.meta.posSignature = state.posSignature;
end
end

function profile = applyDefaults(profile)
if ~isfield(profile, 'version') || isempty(profile.version)
    profile.version = 1;
end
if ~isfield(profile, 'schema') || isempty(profile.schema)
    profile.schema = 'visualisationProfile';
end
if ~isfield(profile, 'settings') || ~isstruct(profile.settings)
    profile.settings = struct();
end
if ~isfield(profile, 'species') || ~isstruct(profile.species)
    profile.species = struct();
end
if ~isfield(profile, 'policies') || ~isstruct(profile.policies)
    profile.policies = struct();
end
if ~isfield(profile, 'export') || ~isstruct(profile.export)
    profile.export = struct();
end
if ~isfield(profile, 'meta') || ~isstruct(profile.meta)
    profile.meta = struct();
end

settingsDefaults = struct( ...
    'mode', "ionic", ...
    'splitIsotope', false, ...
    'splitCharge', false, ...
    'showUnranged', true, ...
    'sampleValue', 1, ...
    'markerSizeGlobal', 15, ...
    'markerSize', 15, ...
    'randomSeed', 1, ...
    'rotationAngle', 0, ...
    'fixBoundingBox', true, ...
    'clipToBoundingBox', true, ...
    'bboxUseSpan', false, ...
    'showAxisLabels', true, ...
    'showAxisTicks', true, ...
    'showScaleCube', false, ...
    'scaleCubeSize', 10, ...
    'searchFilter', "", ...
    'sortMode', "none");
profile.settings = fillMissingFields(profile.settings, settingsDefaults);

if ~isfield(profile.settings, 'boundingBox')
    profile.settings.boundingBox = struct();
end
if ~isfield(profile.settings, 'camera')
    profile.settings.camera = struct();
end

policyDefaults = struct( ...
    'speciesPolicy', "expand", ...
    'bboxPolicy', "dataExtent", ...
    'cameraPolicy', "keep", ...
    'rotationPolicy', "keep", ...
    'seedPolicy', "fixed");
profile.policies = fillMissingFields(profile.policies, policyDefaults);

exportDefaults = struct( ...
    'imageExtension', '.png', ...
    'imageResolution', 300, ...
    'turntableStepDeg', 0.5, ...
    'turntableFrameRate', 20, ...
    'turntableTotalAngleDeg', 360);
profile.export = fillMissingFields(profile.export, exportDefaults);

if ~isfield(profile.meta, 'createdAt') || isempty(profile.meta.createdAt)
    profile.meta.createdAt = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

profile.species = normalizeSpecies(profile.species, profile.settings.markerSizeGlobal);
end

function species = normalizeSpecies(species, markerSizeGlobal)
if nargin < 2 || isempty(markerSizeGlobal)
    markerSizeGlobal = 15;
end
if ~isfield(species, 'name') || isempty(species.name)
    species.name = strings(0, 1);
else
    species.name = string(species.name(:));
end
n = numel(species.name);

if ~isfield(species, 'visible') || numel(species.visible) ~= n
    species.visible = true(n, 1);
else
    species.visible = logical(species.visible(:));
end

if ~isfield(species, 'fraction') || numel(species.fraction) ~= n
    species.fraction = ones(n, 1);
else
    species.fraction = double(species.fraction(:));
end

if ~isfield(species, 'markerSize') || numel(species.markerSize) ~= n
    species.markerSize = repmat(markerSizeGlobal, n, 1);
else
    species.markerSize = double(species.markerSize(:));
end

if ~isfield(species, 'markerSizeLinked') || numel(species.markerSizeLinked) ~= n
    tol = max(1e-9, 1e-6 * max(1, markerSizeGlobal));
    species.markerSizeLinked = abs(species.markerSize - markerSizeGlobal) <= tol;
else
    species.markerSizeLinked = logical(species.markerSizeLinked(:));
end

if ~isfield(species, 'color') || size(species.color, 1) ~= n || size(species.color, 2) ~= 3
    species.color = NaN(n, 3);
else
    species.color = double(species.color);
end
end

function outStruct = fillMissingFields(inStruct, defaults)
outStruct = inStruct;
fieldNames = fieldnames(defaults);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if ~isfield(outStruct, fieldName) || isempty(outStruct.(fieldName))
        outStruct.(fieldName) = defaults.(fieldName);
    end
end
end
