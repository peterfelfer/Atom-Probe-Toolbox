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


%% writing pos data ==> APT specific!
numEntries = height(pos);
isDecomposed = numEntries > max(pos.ionIdx); % check if pos file is decomposed

if isDecomposed
    dataPath = '/atomProbeTomography/reconstruction/atom/';
else
    dataPath = '/atomProbeTomography/reconstruction/ion/';
end

posColumnNames = pos.Properties.VariableNames;
for col = 1:width(pos)    
    % formatting data for hdf5 wrtie if necessary
    if ismember(posColumnNames{col},{'x','y','z'}) %float type coords
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath posColumnNames{col}],[numEntries 1]);
        h5write(fileName,[dataPath posColumnNames{col}],data);
        h5writeatt(fileName,[dataPath posColumnNames{col}],'unit','nm','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'chargeState','isotope','ionComplexity'}) %float type coords
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath posColumnNames{col}],[numEntries 1]);
        h5write(fileName,[dataPath posColumnNames{col}],data);
        
    elseif ismember(posColumnNames{col},{'mc'})
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath 'massToChargeState'],[numEntries 1]);
        h5write(fileName,[dataPath 'massToChargeState'],data);
        h5writeatt(fileName,[dataPath 'massToChargeState'],'unit','Da','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'ionIdx'})
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath 'fieldEvaporationSequenceIndex'],[numEntries 1]);
        h5write(fileName,[dataPath 'fieldEvaporationSequenceIndex'],data);
        
    elseif ismember(posColumnNames{col},{'ion','atom'}) % categorical arrays converted to string arrays
        %data = string(table2array(pos(:,col)));
        %dataType = 'string';
        data = [1];
        dataType = 'single';
        
        
        % epos variables
    elseif ismember(posColumnNames{col},{'tof'})
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/','unit','ns','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'VDC'})
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/standingVoltage',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/standingVoltage',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/standingVoltage','unit','V','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'VP'})
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/','unit','V','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'detx'})
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/x/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/x/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/x/','unit','mm','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'dety'})
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/y/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/y/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/y/','unit','mm','TextEncoding','UTF-8');
        
    end

    
    % write attributes of the column
    % h5writeatt(fileName,dataSetName,'unit',pos.Properties.VariableUnits{col});
    % h5writeatt(fileName,dataSetName,'description',pos.Properties.VariableDescriptions{col});
end


%% writing ions 




%% writing ranges


%% close file
%H5F.close(fid);