function hdf5FileCreateFromMetaDataList(fileName,metaData)
% hdf5FileCreateFromMetaDataList creates a hdf5 file with all the metadata
% that is provided in the metaData variable. This variable is a cell array
% with variable name, including path in the first column, variable value in
% the second and variable unit in the third column. Experimental data and
% data from the analysis can then be added to this file. 
%
%
%   hdf5FileCreateFromMetaDataList(fileName,metaData)
%
%   INPUTS:
%   fileName        full file name including path as string or char array
%   
%   metaData        cell array of metadata in the form {path/varName, value, unit}


%% creating hdf5 file 
fid = H5F.create(fileName);


%% creating dataset groups
% extracting groups form metadata list
% group list needs to be created such that a group is created before its
% subgroup. This demands extracting the individual components of each path

% deleting the individual variable names
groupsTmp = string(metaData(:,1)); 
delPos = regexp(groupsTmp,"/");
delPos = cellfun(@(x) x(end),delPos);
groupsTmp = extractBefore(groupsTmp,delPos+1);
groupsTmp = unique(groupsTmp);

% extracting the group hierachy
delPos = regexp(groupsTmp,"/");
groups = [];
for gr = 1:length(groupsTmp)
    groupDepth = length(delPos{gr});
    tmp = repmat(groupsTmp(gr),groupDepth,1);
    tmp = extractBefore(tmp,(delPos{gr})');
    groups = [groups; tmp];
end

% format group list
groups = unique(groups);
groups = sort(groups);
groups(groups == "") = [];
groups = groups + "/";

% cretaing individual groups
for gr = 1:length(groups)
    groupID(gr) = H5G.create(fid,groups(gr),...
        'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(groupID(gr));
end



%% writing metadata
for m = 1:size(metaData,1)
    % formatting of metadata for writing
    if ~isempty(metaData{m,2})
        
        % get dataset path in hdf5 file from metadata text
        idxLastSlash = find(metaData{m,1} == '/');
        idxLastSlash = max(idxLastSlash);
        
        path = metaData{m,1}(1:idxLastSlash);
        attribute = metaData{m,1}(idxLastSlash+1:end);
        
        if isdatetime(metaData{m,2}) % for dates
            data = datestr(metaData{m,2});
            
        elseif iscategorical(metaData{m,2}) % for enums (categoricals in Matlab)
            data = string(metaData{m,2});
            
        elseif islogical(metaData{m,2})
            if metaData{m,2}
                data = 'true';
            else 
                data = 'false';
            end
            
        else
            data = metaData{m,2};
        end
        
        h5writeatt(fileName,path,attribute,data,'TextEncoding','UTF-8');
    end      
end

H5F.close(fid);