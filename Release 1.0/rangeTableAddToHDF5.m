function rangeTableAddToHDF5(fileName,rangeTable)
% hdf5rangeTableAdd can be used to add range information to an existing HDF5 file
%
% rangeTableAddToHDF5(fileName,rangeTable)
%
% INPUT    
% fileName:         full file name including path as string or char array
%
% rangeTable:       table of ranges usually as created by extraction from a
%                   mass spectrum plot via rangesExtractFromMassSpec(spec)
%
% For each range, written are:
%
% /atomProbeTomography/massToChargeRange/1/name [string] % e.g. Fe2 O3
% /atomProbeTomography/massToChargeRange/1/ion [string] % e.g. 56Fe2 16O3 ++
% /atomProbeTomography/massToChargeRange/1/begin [float32] % e.g. 32.34
% /atomProbeTomography/massToChargeRange/1/end [float32] % e.g. 33.34
% /atomProbeTomography/massToChargeRange/1/color [float32] % RGB 0-1 e.g. 0.23 , 0.4, 1.0
%
% All ranges information is written as attributes in consecutive groups
% /atomProbeTomography/massToChargeRange/1/name
% /atomProbeTomography/massToChargeRange/2/name
% etc...
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

numRng = height(rangeTable);

basePath = "/atomProbeTomography/massToChargeRange/";

% open the hdf5 file
fid = H5F.open(fileName,'H5F_ACC_RDWR','H5P_DEFAULT');

for r = 1:numRng
    path = basePath + num2str(r);
    
    % create group
    groupID(r) = H5G.create(fid,path,...
        'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(groupID(r));
    
    % write range name
    h5writeatt(fileName,path, "name", string(rangeTable.rangeName(r)));
    
    % write range ion
        % check: if rangeName is an individual Name and does not have a specific
        % ion 
    if ~isempty(rangeTable.ion{r})
        rangeIon = string(ionConvertName(rangeTable.ion{r},rangeTable.chargeState(r)));
        h5writeatt(fileName,path, "ion", rangeIon);
    else
        rangeIon = string(rangeTable.rangeName(r));
        h5writeatt(fileName,path, "ion", rangeIon);
    end
    
    % write range bounds
    h5writeatt(fileName,path, "begin", rangeTable.mcbegin(r));
    h5writeatt(fileName,path, "end", rangeTable.mcend(r));
    
    % write visualisation color
    h5writeatt(fileName,path, "color", rangeTable.color(r,:));
    
end


H5F.close(fid);