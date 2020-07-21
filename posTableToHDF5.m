function posTableToHDF5(fileName,pos,metaData,ionTable,rangeTable)
% posToHDF5 exports a pos variable to an HDF 5 file, optionally including
% metadata about ranges and ions identified in the data.
%
% fileName as text
%
% pos as table with entries ionIdx, x, y, z, mc and more if needed
%
% metadata as cell aray with {variableName, value, unit}
%
% TODO: implement export of *.epos
%       implementation of export of categorical arrays


%% creating hdf5 file 
fid = H5F.create(fileName);


%% creating dataset groups

% extracting groups form metadata list
groupsTmp = string(meta{:,1}); 
delPos = regexp(groupsTmp,"/");
delPos = cellfun(@(x) x(end),delPos);
groups = extractBefore(groupsTmp,delPos+1);

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
        
        h5writeatt(fileName,path,attribute,data);
    end
    
    
end


%% writing pos data
posColumnNames = pos.Properties.VariableNames;
for col = 1:width(pos)
    dataSetName = [dataPath posColumnNames{col}]; % a dataset in HDF5 is one consistent piece of data. i.e. each column of the pos variable is one 'dataset'
    
    % formatting data for hdf5 wrtie if necessary
    if ismember(posColumnNames{col},{'x','y','z','mc'}) %float types
        data = table2array(pos(:,col));
        dataType = 'single';
        dataPath = '/atomProbeTomography/reconstruction/';
        type_id = H5T.copy('H5T_NATIVE_SINGLE');
        
        H5D.write(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL',plist,data);
        H5D.close(dset_id);
        
    elseif ismember(posColumnNames{col},{'ionIdx','chargeState','isotope','ionComplexity'}) %integer types
        data = table2array(pos(:,col));
        dataType = 'uint32';
        
    elseif ismember(posColumnNames{col},{'ion','atom'}) % categorical arrays converted to string arrays
        %data = string(table2array(pos(:,col)));
        %dataType = 'string';
        data = [1];
        dataType = 'single';
        
    end

    
    % write attributes of the column
    % h5writeatt(fileName,dataSetName,'unit',pos.Properties.VariableUnits{col});
    % h5writeatt(fileName,dataSetName,'description',pos.Properties.VariableDescriptions{col});
end


%% writing ions 




%% writing ranges


%% close file
H5F.close(fid);