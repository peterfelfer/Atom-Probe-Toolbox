%[text] # HDF5 data I/O and file searching
%[text] **reading of a metadata text file, creating a hdf5 file, adding information, performing analysis**
%[text] A convenient way to make APT data searchable is the use of the HDF5 file format to store experimental data, analysis data, and metadata in one file. Large collections of these files can then be searched for certain properties and analysed using scripts in the toolbox. At current, the toolbox supports the following data to be imported/exported to HDF5:
%[text] - pos data (ionIdx, x, y, z, m/c, tof, detx, dety, VDC, VP, pulse delta / pulse idx, hit multiplicity)
%[text] - range data (name, ion, mcbegin, mcend, color)
%[text] - ion data (ion type, chargestate)
%[text] - images (e.g. SEM image of tip)
%[text] - % mesh objects (e.g. isosurfaces)
%[text] - % ROIs \
%[text] The data names and attributes are defined in the metadata files or can be supplied in Matlab. The metadata reference file for the toolbox can be found in meta.
%[text] **Reminder: HDF5 files are meant for storage and permanence. Therefore deletion utilities do not exist.**
%%
%[text] ## Creating an HDF5 file from metadata
%[text] As a first step, the user has to specify the path and the file name of the metadata's text file. The user can specify both manually or by executing the code \[path fileName\] = uigetfile;. In case of the latter, the user can select a metadata text file within the popped up selection window. 
%[text] NOTE: The user may have to change the displayed file type in the buttom right of the selection window, so that also the wanted \*.txt files are displayed.
[fileName path] = uigetfile;
meta = metaDataReadTextFile([path fileName]); % input of the metadata text file and its file path %[output:2862c64b]
  %[control:button:1cc0]{"position":[1,2]}
%%
%[text] Subsequently to the conversion into a cell array, the metadata list can be saved as an HDF5 file with the function *hdf5FileCreateFromMetaDataList*. By doing so, the metadata is readable for an HDF5 reader and can be used for data mining and other purposes.
%[text] As first input, the wanted path and file name must be defined as a string (e.g., '\\Desktop\\test.h5'). Secondly, the previously generated metadata cell array must be given as input. 
%[text] NOTE (**for all upcoming sections**): If the file should be saved in the current folder, the specification of the folder path is dispensable. Furthermore, the suffixes \*.h5, \*.hdf5 are both possible.
fileName = "4322_RC15_F1_1.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:8d82]{"position":[12,31]}
metaData = meta; % definition of the previously generated metadata cell array
hdf5FileCreateFromMetaDataList(fileName,metaData); % generates and saves metadata as an HDF5 file
  %[control:button:922a]{"position":[1,2]}
%[text] 
%%
%[text] ## Creating / writing to and reading from an HDF5
%[text] With an HDF5 file already present, it is possible to add data (e.g., pos/epos data, range/ion table, images) to it and also to read it and extract data from the file.
%[text] ### Working with a pos/epos variable
%[text] #### Adding a pos/epos table to an HDF5 file
%[text] The function *posTableAddToHDF5* enables the user to add pos or epos data to an existing HDF5 file with the appropriate groups structure. As inputs the name of the file, in which the pos data should be added and the pos variable itself must be specified.
%[text] The file name must be the full file name including the path as a string or char array.
%[text] NOTE: If the file name does not yet exist in the specified path, a new HDF5 file will be created.
fileName = "4322_RC15_F1_1.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:8d19]{"position":[12,31]}
posTableAddToHDF5(fileName,pos); % adds a pos or epos table to the HDF5 file
  %[control:button:9ed8]{"position":[1,2]}
%[text] NOTE: Again, if the file should be saved in the current folder, the specification of the folder path is dispensable. Furthermore, the suffixes \*.h5, \*.hdf5 are both possible.
%%
%[text] ### Working with a range table
%[text] #### Adding a range table to an HDF5 file
%[text] For the addition of a range table to an HDF5 file, the range table needs to be created by executing the function *rangesExtractFromMassSpec* (functionality explained in the live script **FirstSteps**). As input the current mass spectrum (often called *spec* in the workspace) is needed. Subsequently, the function *rangeTableAddToHDF5* adds the data located in the range table to the HDF5 file. The file name (*fileName*, including the path as a string or char array) and the range table (*rangeTable*) need to be specified. 
%rangeTable = rangesExtractFromMassSpec(spec); % extracts ranges from a mass spec and creates a range table
%fileName = "test.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:95d4]{"position":[13,22]}
rangeTableAddToHDF5(fileName,rangeTable); % adds a range table to an existing HDF5 file
  %[control:button:8466]{"position":[1,2]}
%%
%[text] #### Reading a range table from an HDF5 file
%[text] With the function *rangeTableFromHDF5*, the user can extract the ranges and create a range table in the workspace. As input - as usual - only the file name including the path as a string or char array is needed.
fileName = "test.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:358f]{"position":[12,21]}
rangeTable = rangeTableFromHDF5(fileName); % creates a range table in the workspace from an HDF5 file
  %[control:button:3a2c]{"position":[1,2]}
%[text] 
%%
%[text] ### Working with an ion table
%[text] #### Adding an ion table to an HDF5 file
%[text] Analogous to the feature presented in the previous sections, the user can also add an ion table to an existing HDF5 file by executing the function *ionTableAddToHDF5*.
%[text] In order to add an ion table, the user has to extract the ions from a current mass spectrum. This procedure is explained in the live script **FirstSteps**. 
ionTable = ionsExtractFromMassSpec(spec); % creates an ion table of a current mass spectrum
%fileName = "test.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:8f73]{"position":[13,22]}
ionTableAddToHDF5(fileName,ionTable); % adds an ion table to an existing HDF5 file
  %[control:button:07ee]{"position":[1,2]}
%%
%[text] #### Reading an ion table from an HDF5 file
%[text] With the function *ionTableFromHDF5*, the user can create an ion table in the workspace on the basis of an HDF5 file. As input the file name, including the path as a string or char array, is needed.
fileName = "test.h5"; % specification of the file name, including its path, no inverted commas needed %[control:editfield:31a7]{"position":[12,21]}
ionTable = ionTableFromHDF5(fileName); % creates an ion table in the workspace from an HDF5 file
  %[control:button:775c]{"position":[1,2]}
%[text] ### 
%%
%[text] ### Adding/reading an image to the HDF5 file 
%[text] #### Adding an image to the HDF5 file
%[text] For adding an image to a HDF5 file, the user needs to load the image into the Matlab workspace. Afterwards, by using *imageAddToHDF5,* the user can add the image to the HDF5 file. The user can add an image type or contrast mechanism to to give more details about the image. 
image = imread('FDM.jpg'); % or .jpg .png ... %[control:editfield:832f]{"position":[16,25]}
fileName = "test.h5"; % specification of the HDF5 file name, including its path, no inverted commas needed %[control:editfield:7935]{"position":[12,21]}
imageType = 'Matlab'; % image Type can be the kind of microscope/camera with whom the image was taken %[control:editfield:5dd0]{"position":[13,21]}
contrastMechanism = 'No'; % e.g. 'Secondary Electron' %[control:editfield:5fa4]{"position":[21,25]}
horizontalPixelScale = 1;
verticalPixelScale = 1;

imageAddToHDF5(fileName,image,imageType,contrastMechanism,...
    horizontalPixelScale,verticalPixelScale);
  %[control:button:46be]{"position":[1,2]}
%%
%[text] #### Reading an image of a HDF5 file
%[text] To read an image out of a HDF5 file, the number of images stored inside the file is important.
fileName = 'test.h5'; % specification of the HDF5 file name, including its path, no inverted commas needed %[control:editfield:9cca]{"position":[12,21]}
imageInfo = h5info(fileName, '/sample');
imageInfo = imageInfo.Groups.Datasets
numberImg = size(imageInfo,1)
  %[control:button:28d9]{"position":[1,2]}
%[text] The variable *imageInfo* contains the properties of the images including the attributes provided by the *imageAddToHDF5* function. 
%%
%[text] Now, the user can choose which image should be extracted form the HDF5 file.
numImg = 3; % the number of the desired image %[control:slider:5cfa]{"position":[10,11]}
image = imageFromHDF5(fileName, numImg);
% imshow(image); % opens the image
  %[control:button:82d9]{"position":[1,2]}
%%
%[text] ## Searching folder structures with HDF5 files
%[text] One of the benefits when using HDF5 files is the possibility to perform data mining operations. Since all important (meta-)data is stored within an HDF5 file, the user can search for specific parameters (e.g., temperature, at which the experiment was conducted) within a vast number of files. A list of such files, which fulfils these properties, is generated as an output.
%[text] ### Finding \*.h5 files that fulfil certain properties
%[text] The function *hdf5FileFindByAttribute* allows the scanning of vast number of \*.h5 files with respect to certain properties. As input, the folder (*topFolder*), in which the files are stored (including subfolders), must be specified. The attributes are parsed by combinations of keyword and value. For numeric values, comparisons can be used. If comparisons are used, the value input needs to be char type examples:
%[text] **String type:** '/sample/materialType', 'aluminium'
%[text] **Numeric with comparison:** '/atomProbeTomography/experiment/specimenTemperature', '\<= 40'
%[text] **Numeric with exact value:** '/atomProbeTomography/experiment/specimenTemperature', 40
topFolder = uigetdir(); % opens a selection window
%[text] In the following, the user can scan the folder either on the basis of one, two, or three attributes. In order to do so, it is necessary to uncomment the desired section (by marking the section and pressing the '%x'-sign or hitting Ctrl+T on the keyboard). 
%[text] The output *files* is a structure, which contains various fields (e.g., name, folder, date, and others).
%[text] For one attribute-value-pair:
attribute1 = '/atomProbeTomography/experiment/specimenTemperature'; % specification of first attribute %[control:editfield:3ea6]{"position":[14,67]}
value1 = '40';     % specification of first value corresponding to first attribute %[control:editfield:89e2]{"position":[10,14]}
files = hdf5FileFindByAttribute(topFolder,...   % generates output on the basis of attribute-value input arguments                               
    attribute1, value1);
  %[control:button:4c78]{"position":[1,2]}
%[text] For two attribute-value-pairs:
% attribute1 = ''; % specification of first attribute %[control:editfield:7d70]{"position":[16,18]}
% value1 = '';     % specification of first value corresponding to first attribute %[control:editfield:2cb2]{"position":[12,14]}
% attribute2 = ""; % specification of second attribute  %[control:editfield:9879]{"position":[16,18]}
% value2 = "";     % specification of second value corresponding to second attribute %[control:editfield:359e]{"position":[12,14]}
% files = hdf5FileFindByAttribute(topFolder,...   % generates output on the basis of attribute-value input arguments
%     attribute1, value1,...
%     attribute2, value2);
%   %[control:button:045c]{"position":[3,4]}
%[text] For three attribute-value-pairs:

% attribute1 = ''; % specification of first attribute %[control:editfield:1c45]{"position":[16,18]}
% value1 = '';     % specification of first value corresponding to first attribute %[control:editfield:384a]{"position":[12,14]}
% attribute2 = ""; % specification of second attribute  %[control:editfield:4b1e]{"position":[16,18]}
% value2 = "";     % specification of second value corresponding to second attribute %[control:editfield:7a7d]{"position":[12,14]}
% attribute3 = ""; % specification of third attribute %[control:editfield:9225]{"position":[16,18]}
% value3 = "";     % specification of third value corresponding to third attribute %[control:editfield:1982]{"position":[12,14]}
% files = hdf5FileFindByAttribute(topFolder,...   % generates output on the basis of attribute-value input arguments
%     attribute1, value1,...
%     attribute2, value2,...
%     attribute3, value3);
%   %[control:button:3085]{"position":[3,4]}

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":18.9}
%---
%[control:button:1cc0]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:8d82]
%   data: {"defaultValue":"\"4322_RC15_F1_1.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:922a]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:8d19]
%   data: {"defaultValue":"\"4322_RC15_F1_1.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:9ed8]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:95d4]
%   data: {"defaultValue":"\"test.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:8466]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:358f]
%   data: {"defaultValue":"\"test.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:3a2c]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:8f73]
%   data: {"defaultValue":"\"test.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:07ee]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:31a7]
%   data: {"defaultValue":"\"test.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:button:775c]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:832f]
%   data: {"defaultValue":"'FDM.jpg'","label":"Edit field","run":"Section","valueType":"Char"}
%---
%[control:editfield:7935]
%   data: {"defaultValue":"\"test.h5\"","label":"fileName","run":"Nothing","valueType":"String"}
%---
%[control:editfield:5dd0]
%   data: {"defaultValue":"'Matlab'","label":"Edit field","run":"Section","valueType":"Char"}
%---
%[control:editfield:5fa4]
%   data: {"defaultValue":"'No'","label":"Edit field","run":"Section","valueType":"Char"}
%---
%[control:button:46be]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:9cca]
%   data: {"defaultValue":"'test.h5'","label":"Edit field","run":"Section","valueType":"Char"}
%---
%[control:button:28d9]
%   data: {"label":"Run","run":"Section"}
%---
%[control:slider:5cfa]
%   data: {"defaultValue":3,"label":"numImg","max":100,"min":1,"run":"Nothing","runOn":"ValueChanging","step":1}
%---
%[control:button:82d9]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:3ea6]
%   data: {"defaultValue":"'\/atomProbeTomography\/experiment\/specimenTemperature'","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:89e2]
%   data: {"defaultValue":"'40'","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:button:4c78]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:7d70]
%   data: {"defaultValue":"''","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:2cb2]
%   data: {"defaultValue":"''","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:9879]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:editfield:359e]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:button:045c]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:1c45]
%   data: {"defaultValue":"''","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:384a]
%   data: {"defaultValue":"''","label":"attribute1","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:4b1e]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:editfield:7a7d]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:editfield:9225]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:editfield:1982]
%   data: {"defaultValue":"\"\"","label":"attribute1","run":"Nothing","valueType":"String"}
%---
%[control:button:3085]
%   data: {"label":"Run","run":"Section"}
%---
%[output:2862c64b]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Index exceeds the number of array elements. Index must not exceed 0.\n\nError in <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('metaDataReadTextFile', 'Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\metaDataReadTextFile.m', 63)\" style=\"font-weight:bold\">metaDataReadTextFile<\/a> (<a href=\"matlab: opentoline('Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\metaDataReadTextFile.m',63,0)\">line 63<\/a>)\n    varFormat = fileByLine{li}(isBrackets(1)+1:isBrackets(2)-1);"}}
%---
