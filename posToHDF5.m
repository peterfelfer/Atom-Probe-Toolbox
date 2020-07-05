function posToHDF5(fileName,pos,metaData,ionTable,rangeTable)
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


%% writing pos data
posColumnNames = pos.Properties.VariableNames;
hitData = '/hitData/';
for col = 1:width(pos)
    dataSetName = [hitData posColumnNames{col}]; % a dataset in HDF5 is one consistent piece of data. i.e. each column of the pos variable is one 'dataset'
    
    % formatting data for hdf5 wrtie if necessary
    if ismember(posColumnNames{col},{'x','y','z','mc'}) %float types
        data = table2array(pos(:,col));
        dataType = 'single';
        
    elseif ismember(posColumnNames{col},{'ionIdx','chargeState','isotope','ionComplexity'}) %integer types
        data = table2array(pos(:,col));
        dataType = 'uint32';
        
    elseif ismember(posColumnNames{col},{'ion','atom'}) % categorical arrays converted to string arrays
        %data = string(table2array(pos(:,col)));
        %dataType = 'string';
        data = [1];
        dataType = 'single';
        
    end
    % write the individual columns
    h5create(fileName,dataSetName,size(data),'DataType',dataType);
    h5write(fileName,dataSetName,data);
    
    % write attributes of the column
    h5writeatt(fileName,dataSetName,'unit',pos.Properties.VariableUnits{col});
    %h5writeatt(fileName,dataSetName,'description',pos.Properties.VariableDescriptions{col});
end


%% writing metadata
for m = 1:size(metaData,1)
    
    % formatting of metadata for writing
    if isempty(metaData{m,2})
        data = '';
    else
        if isdatetime(metaData{m,2})
            data = datestr(metaData{m,2});
        else
            data = metaData{m,2};
        end
    end
    
    h5writeatt(fileName,'/',metaData{m,1},data);
end