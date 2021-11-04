function newMeta = metaDataReadLEAPAttributes(meta, exNum, exList)
%UNTITLED Summary of this function goes here
%   This function adds the output parameter from the LEAP to the meta Data
%   file list 
% INPUT
% meta = meta data List 
% exNum = specificNumber of Experiment
% exList = experiment List with the output Data from the LEAP 


% get important row out of exList
intRow = exList(exList.ExperimentID == exNum, :);

% sort Variables 


meta(contains(meta(:,1),'experiment/uniqueIdentifier'),2) = intRow.RunID;
meta(contains(meta(:,1),'specimen/uniqueIdentifier'),2) = intRow.Specimen;






app.MetaWS(contains(app.MetaWS(:,1),'/project/identifier'),2)







end

