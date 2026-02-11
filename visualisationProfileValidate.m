function visualisationProfileValidate(visualisationProfile)
% VISUALISATIONPROFILEVALIDATE Validate visualisation profile schema.
%
% visualisationProfileValidate(visualisationProfile)

if ~isstruct(visualisationProfile) || ~isscalar(visualisationProfile)
    error('visualisationProfileValidate:invalidType', ...
        'visualisationProfile must be a scalar struct.');
end

requiredTopFields = {'version', 'schema', 'settings', 'species', 'policies', 'export'};
for i = 1:numel(requiredTopFields)
    fieldName = requiredTopFields{i};
    if ~isfield(visualisationProfile, fieldName)
        error('visualisationProfileValidate:missingField', ...
            'Missing required field ''%s''.', fieldName);
    end
end

if ~strcmp(char(string(visualisationProfile.schema)), 'visualisationProfile')
    error('visualisationProfileValidate:invalidSchema', ...
        'Invalid schema ''%s''.', char(string(visualisationProfile.schema)));
end

if ~isstruct(visualisationProfile.settings) || ~isscalar(visualisationProfile.settings)
    error('visualisationProfileValidate:invalidSettings', ...
        'settings must be a scalar struct.');
end
if ~isstruct(visualisationProfile.species) || ~isscalar(visualisationProfile.species)
    error('visualisationProfileValidate:invalidSpecies', ...
        'species must be a scalar struct.');
end

species = visualisationProfile.species;
if ~isfield(species, 'name')
    error('visualisationProfileValidate:missingSpeciesNames', ...
        'species.name is required.');
end

names = string(species.name(:));
n = numel(names);

checkVectorLength(species, 'visible', n, true);
checkVectorLength(species, 'fraction', n, true);
checkVectorLength(species, 'markerSize', n, true);
checkVectorLength(species, 'markerSizeLinked', n, true);

if ~isfield(species, 'color') || isempty(species.color)
    error('visualisationProfileValidate:missingSpeciesColor', ...
        'species.color is required and must be Nx3.');
end
if size(species.color, 1) ~= n || size(species.color, 2) ~= 3
    error('visualisationProfileValidate:invalidSpeciesColor', ...
        'species.color must be Nx3 with N = numel(species.name).');
end

end

function checkVectorLength(species, fieldName, n, required)
if nargin < 4
    required = true;
end
if ~isfield(species, fieldName)
    if required
        error('visualisationProfileValidate:missingSpeciesField', ...
            'species.%s is required.', fieldName);
    end
    return;
end
value = species.(fieldName);
if numel(value) ~= n
    error('visualisationProfileValidate:invalidSpeciesFieldLength', ...
        'species.%s must have %d elements.', fieldName, n);
end
end
