%[text] # Workflow HDF5 Data Analysis
%[text] This workflow shows how to analyse the database of HDF5 files. 
%[text] The following points are explained and can be calculated with this workflow:
%[text] - Search for attributes form the backbone structure
%[text] - Find the range McBegin und McEnd of a certain element
%[text] - IonList of every ion that is stored in the HDF5 database
%[text] - Calculate the concentration of each data set and store it in a table
%[text] - Extract the concentration of a specific Element \
%%
%[text] Storage folder for HDF5 files
database = 'L:\existingHDF5files\'; %[control:editfield:27fd]{"position":[12,35]}
cd(database); % make the database to your current folder
%%
%[text] ### Search for attributes from the backbone structure
%[text] Load the metadata file, that hast the structure like the HDF5file then add the attributes you want to look for.
metaOri = metaDataReadTextFile('L:\experiment metadata list v3.txt');
attr1 = 'pulseType'; %[control:editfield:2652]{"position":[9,20]}
attr2 = 'identifier'; %[control:editfield:54f6]{"position":[9,21]}
%attr3 = 'ambientTemperature'; %[control:editfield:16fc]{"position":[10,30]}

attrList = {attr1 attr2}; % add if you want more attributes

for i=1:length(attrList)
    indexC = strfind(metaOri(:,1),attrList{i});
    index = find(not(cellfun('isempty',indexC)));
    attributes{1,i} = metaOri{index};
end
%[text] Search for the attributes
result = hdf5Search(database,attributes);
clearvars -except result;
  %[control:button:6a25]{"position":[1,2]}
%[text] ### 
%%
%[text] ### Find the range McBegin und McEnd of a certain element
%[text] 
fileList = dir('*.h5'); % create a List with all HDF5 files
loc1 = '/atomProbeTomography/massToChargeRange/'; % find the place in the backbone structure where the range information is stored
attrElement = 'Al'; % type in the element %[control:editfield:99f6]{"position":[15,19]}
attr1Name = 'begin';
attr2Name = 'end';

for s = 1:length(fileList)
        
        fileName = fileList(s).name;
        fileInfo = h5info(fileName, loc1);

    for i = 1:length(fileInfo.Groups) %looks inside an H5 file for information
        loci = [loc1, num2str(i)];
        nameElement = h5readatt(fileName, loci, 'name');
        if strcmp(nameElement, attrElement) >= 1
            attr1 = h5readatt(fileName, loci, attr1Name);
            attr2 = h5readatt(fileName, loci, attr2Name);
            attr1List(s) = attr1;
            attr2List(s) = attr2;
            break
        end
        
    end
end

result = table();  
fileListCell = struct2cell(fileList);
fileListCell2 = fileListCell';
result = addvars(result, fileListCell2(:,1), attr1List', attr2List');
result.Properties.VariableNames = {'name', attr1Name, attr2Name};
clearvars -except result;
  %[control:button:8670]{"position":[1,2]}
%%
%[text] ### IonList of every ion that is stored in the HDF5
%[text] The code extracts from each HDF5 file the ionTable and deletes duplicate entries
fileList = dir('*.h5');
ionTableAllH5 = table();

for i = 1:length(fileList)
    fileName = fileList(i).name; 
    loc = '/atomProbeTomography/identifiedIon/';
    fileInfo = h5info(fileName, loc);
    % check if the HDF5 list has a stored ionList
    if ~isempty(fileInfo.Groups)
        ionTable = ionTableFromHDF5(fileName);
        ionTableAllH5 = [ionTableAllH5; ionTable]; 
    end 
end
%[text] Delete doubled entries and the color
ionTableAllH5 = removevars(ionTableAllH5, 'ion');
ionTableH5 = removevars(ionTableAllH5, 'color');
ionTableH5 = unique(ionTableH5);
clearvars -except ionTableH5

% if you want to have just a categorical array with the unique ions 
% ionTableH5UniqueIon = unique(ionTableH5.ionName);



  %[control:button:7ab0]{"position":[1,2]}
%%
%[text] ### Calculate the concentration of each data set and store it in a table
fileList = dir('*.h5');
% create result table
result = table();       
fileListCell = struct2cell(fileList)';
result = addvars(result, fileListCell(:,1), 'NewVariableNames','fileName');
result(4:end,:) = []; %%%%%%%%%%%%%%%%%%%%%%%

for i= 1:3  %length(fileList) %%%%%%%%
    % pos erstellen
    fileName = fileList(i).name; 
    posIn = posTableFromHDF5(fileName);
    rangeTable = rangeTableFromHDF5(fileName);
    pos = posAllocateRange(posIn, rangeTable, 'decompose');
    %calculate concentration
    conc = posCalculateConcentrationSimple(pos,37,{'Ga', 'unranged'},fileName);
    concList{i} = conc;
    clearvars -except concList result fileList 
end
result = addvars(result, concList', 'NewVariableNames','concentration');

  %[control:button:8a8b]{"position":[1,2]}
%[text] #### Extract the concentration of a specific Element
%[text] After the result table with the entire concentration is calculated - the concentration of a specific element is extracted.
resultElement = table();
fileListCell = struct2cell(fileList)';
resultElement = addvars(resultElement, fileListCell(:,1), 'NewVariableNames','fileName');
resultElement(4:end,:) = [];

for j = 1:3 %length(fileList) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    concHspecific = result.concentration{j,1}.H'; %%%%%%%%%%%%%%%%% change here the H for another element e.g. result.concentration{j,1}.Al
    concH(j, 1) = concHspecific(1); 
    concH(j, 2) = concHspecific(2);
    concH(j, 3) = concHspecific(3);
    
end
resultElement = addvars(resultElement, concH(:,1), concH(:,2), concH(:,3),  'NewVariableNames', {'count' 'concentration' 'variance'});
clearvars -except result resultElement
  %[control:button:2dee]{"position":[1,2]}
%%
%[text] ### Calculate the concentration of certain element in entire data sets
fileList = dir('*.h5');

for i= 1:3 %length(fileList9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pos erstellen
    fileName = fileList(i).name; 
    posIn = posTableFromHDF5(fileName);
    rangeTable = rangeTableFromHDF5(fileName);
    pos = posAllocateRange(posIn, rangeTable, 'decompose');
    conc = posCalculateConcentrationSimple(pos,37,{'Ga', 'unranged'},fileName);
    ion = conc.H;
    ionConv = ion';
    ionAllConc(i, 1) = ionConv(1); 
    ionAllConc(i, 2) = ionConv(2);
    ionAllConc(i, 3) = ionConv(3);
end

result = table();       
fileListCell = struct2cell(fileList)';


result = addvars(result,ionAllConc(:, 1), ionAllConc(:, 2), ionAllConc(:, 2));
result.Properties.VariableNames{1} = 'count';
result.Properties.VariableNames{2} = 'concentration';
result.Properties.VariableNames{3} = 'variance';
  %[control:button:2fc1]{"position":[1,2]}
%[text] 
%[text] 
%[text] 
%[text] 
%[text]  
%[text] 
%[text] 
%[text] 
%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":37.2}
%---
%[control:editfield:27fd]
%   data: {"defaultValue":"''","label":"database","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:2652]
%   data: {"defaultValue":"'pulseType'","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:54f6]
%   data: {"defaultValue":"'contactPerson'","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:16fc]
%   data: {"defaultValue":"'ambientTemperature'","label":"attribute1","run":"Section","valueType":"Char"}
%---
%[control:button:6a25]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:99f6]
%   data: {"defaultValue":"''","label":"attrElement","run":"Nothing","valueType":"Char"}
%---
%[control:button:8670]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:7ab0]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:8a8b]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:2dee]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:2fc1]
%   data: {"label":"Run","run":"Section"}
%---
