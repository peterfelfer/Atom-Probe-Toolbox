function ionTable = ionsExtractFromHDF5(fileName)
% ionsExtractFromHDF5 extracts a list of ions from an HDF5 file and puts
% them in a table , akin to ionsExtractFromMassSpec
%
% USAGE: ionList = ionsExtractFromHDF5(fileName);
%
% INPUT:    fileName = Name of hdf5 file incl. path
%
% OUTPUT:   ionTable = table with the following columns: ionName [string],
%           chargeState [int], ion [categorical], color [float x3]
%           
%           Alternatively false as an output if the file does not contain
%           range information, i.e. is not an atom probe data file
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

% get content structure in HDF5 file
ionInfo = h5info(fileName,'/atomProbeTomography/identifiedIon');
ionInfo = ionInfo.Groups;

numIon = size(ionInfo,1);

% extract each ion from HDF5 file
for i = 1:numIon
    attributes = struct2cell(ionInfo(i).Attributes)';
    
    isIon = strcmp(attributes(:,1),'ion');
    isColor = strcmp(attributes(:,1),'color');
    
    ionFullName = attributes{isIon,end};
    ionColor = attributes{isColor,end};
    
    [it, cs] = ionConvertName(ionFullName);
    
    % populate table
    ionName(i,:) = string(ionConvertName(it.element));
    chargeState(i,:) = cs;
    ion{i,:} = it.element;
    color(i,:) = ionColor;
    
    
end

ionName = categorical(ionName);
ionTable = table(ionName,chargeState,ion,color);