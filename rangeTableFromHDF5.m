function rangeTable = rangeTableFromHDF5(fileName)
% rangeTableFromHDF5 extracts a list of ranges from an HDF5 file and puts
% them in a table , akin to rangesExtractFromMassSpec
%
% rangeTable = rangeTableFromHDF5(fileName);
%
% INPUT    
% fileName:     full file name including path as string or char array
%
% OUTPUT   
% rangeTable:   table with the following columns: ionName [string],
%               chargeState [int], ion [categorical], color [float x3]
%           
% NOTE: Alternatively false as an output if the file does not contain
% range information, i.e. is not an atom probe data file
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

% get content structure in HDF5 file
rangeInfo = h5info(fileName,'/atomProbeTomography/massToChargeRange');
rangeInfo = rangeInfo.Groups;

numRng = size(rangeInfo,1);


% extract each range from HDF5 file
for r = 1:numRng
    attributes = struct2cell(rangeInfo(r).Attributes)';
    
    isName = strcmp(attributes(:,1),'name');
    isIon = strcmp(attributes(:,1),'ion');
    isBegin = strcmp(attributes(:,1),'begin');
    isEnd = strcmp(attributes(:,1),'end');
    isColor = strcmp(attributes(:,1),'color');
    
    % convert ion name to table and charge state value
    [it, cs] = ionConvertName(attributes{isIon,end});
    
    rangeName(r,:) = string(attributes{isName,end});
    chargeState(r,:) = cs;
    mcbegin(r,:) = attributes{isBegin,end};
    mcend(r,:) = attributes{isEnd,end};
    volume(r,:) = 0;
    ion{r,:} = it;
    color(r,:) = attributes{isColor,end};
    
end

rangeName = categorical(rangeName);

% assemble range table
rangeTable = table(rangeName,chargeState,mcbegin,mcend,volume,ion,color);
rangeTable = sortrows(rangeTable,'mcbegin','ascend');