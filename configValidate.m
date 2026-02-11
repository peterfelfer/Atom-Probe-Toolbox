function configValidate(config)
% CONFIGVALIDATE Basic structural validation for toolbox configuration structs.
%
% configValidate(config)

if ~isstruct(config) || ~isscalar(config)
    error('configValidate:invalidType', ...
        'Configuration must be a scalar struct.');
end

fieldNames = fieldnames(config);
for i = 1:numel(fieldNames)
    fieldName = fieldNames{i};
    if startsWith(fieldName, '_')
        error('configValidate:invalidFieldName', ...
            'Top-level field ''%s'' is invalid. Leading underscore is reserved.', fieldName);
    end
end

end
