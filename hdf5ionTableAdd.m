function hdf5ionTableAdd(fileName,ionTable)
% hdf5ionTableAdd can be used to write an ion table into an HDF5 file.
%
% hdf5ionTableAdd(fileName,ionTable)
%
% INPUT
% fileName:    full file name including path as string or char array of existing HDF5 file
%
% ionTable:    table of ions usually as created by extraction from a
%              mass spectrum plot via the function 
%              ionsExtractFromMassSpec(spec)
%
% Written are only the ion types, not the individual isotopic combinations.
% The ions are written in string format in the identifiedIon group:
% /atomProbeTomography/identifiedIon/1/ion [string] = NULL [] % e.g. Fe2 O3 ++
% /atomProbeTomography/identifiedIon/1/color [float32] = NULL [] % RGB 0-1 e.g. 0.23 , 0.4, 1.0
%
% all ion information is written as attributes in consecutive groups
% /atomProbeTomography/identifiedIon/1/ion
% /atomProbeTomography/identifiedIon/2/ion
% etc...

numIon = height(ionTable);

basePath = "/atomProbeTomography/identifiedIon/";

% open the hdf5 file
fid = H5F.open(fileName,'H5F_ACC_RDWR','H5P_DEFAULT');

for i=1:numIon
    path = basePath + num2str(i);    
    
    % create group
    groupID(i) = H5G.create(fid,path,...
        'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(groupID(i));
    
    % write ion name and color
    h5writeatt(fileName,path,"ion",ionConvertName(ionTable.ion{i}, ionTable.chargeState(i)));
    h5writeatt(fileName,path,"color",ionTable.color(i,:));
end


H5F.close(fid);