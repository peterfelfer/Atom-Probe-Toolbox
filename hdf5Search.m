function result = hdf5Search(topFolder,attributes)
%HDF5SEARCH looks in a set of HDF5 files for 
%  

%  topFolder = folder that contains all of the HDF5 files
%  attributes = struct array with attributes, containing location and name
%               of attribute, 
%                 You can extract this attribute array by loading a metaDataFile array in Matlab and extracting the desired variables
%  result = table with the result 
%  Detailed explanation goes here


%% look for HDF5 files in the folder

files = dir([topFolder '/**/*.h5']);



%% look for the variables and their location in the HDF5 file

% get location and attribute name
for i = 1:length(attributes)
        attrNamePos = find(attributes{i} == '/', 1, 'last');
        attrNameTBD = extractAfter(attributes{i},attrNamePos);
        attrLoc = extractBefore(attributes{i},(attrNamePos + 1));
        locList{i,1} = attrLoc;
        nameList{i,1} = attrNameTBD;  
end


%% search for the attributes 

result = table();
fileListCell = struct2cell(files);
fileListCell2 = fileListCell';
result = addvars(result, fileListCell2(:,1), 'NewVariableNames', {'fileName'});
for j = 1:length(attributes) %every single attribute
    loc = locList{j,1};
    attrName = nameList{j,1};
    
    for k = 1:length(files)
        fileName = files(k).name; % find file name
        
        %check if parameter and location are there 
        attrInfo = h5info(fileName, loc); 
        if ~isempty(attrInfo.Attributes)
            allAttributes = {attrInfo.Attributes.Name}';
            attrExist = strcmp(attrName, allAttributes);
        else
            attrExist = 0;
        end

        if sum(attrExist) == 1 
            attr = h5readatt(fileName, loc, attrName);
            attrList{k} = attr; 
        else
            %display('no attribute in .h5 file');
            attrList{k} = 0;
        end

    end
    
    result = addvars(result, attrList', 'NewVariableNames', {attrName});
    
end


% clearvars -except result;
end

