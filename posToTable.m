function pos = posToTable(fileName)
% posToTable returns data from a .pos or .epos file and returns a table with the
% entries and units from an atom probe measurement
%
% pos = posToTable(fileName)
% pos = posToTable()
%
% INPUT
% fileName: is the name of the .pos or .epos file, it is optional. 
%           A dialog box will pop up if no file name is given.
% 
% OUTPUT
% pos: is the variable that contains the entire data from the atom probe for further analysis,
%       table
%
% hint: in case of an epos file as input, be aware to select (*.epos) as
%       displayed data type in the dialog box

if ~exist('fileName','var')
    [file, path, idx] = uigetfile({'*.pos';'*.epos'},'Select a pos file');
    fileName = [path file];
    disp(['file ' file ' loaded']);
end

[~ , ~, ext] = fileparts(fileName);

switch ext
    
    case '.pos'
        idx = 1;
        
    case '.POS'
        idx = 1;
        
    case '.epos'
        idx = 2;
        
    case '.EPOS'
        idx = 2;
        
    otherwise
        idx = 0;
end

if ~idx
    pos = 0;
    return
end

fid = fopen(fileName, 'r');

if idx == 1
    pos = fread(fid, inf, '4*float32', 'b');
    numAtoms = length(pos)/4;
    pos=reshape(pos, [4 numAtoms]);
    pos=double(pos');
    ionIdx = (1:numAtoms)';
    
    pos = table(ionIdx,pos(:,1),pos(:,2),pos(:,3),pos(:,4));
    pos.Properties.VariableNames = {'ionIdx','x','y','z','mc'};
    pos.Properties.VariableUnits = {'1','nm','nm','nm','Da'};
 
    
    
        %% reads .epos
elseif idx == 2
    
    %% Reads through the file made of 9 floats, with 8 byte stride (the two
    %% integers at the end)
    pos = fread(fid, inf, '9*float32', 8 ,'ieee-be');
    
    
    %% reads the pulse info from epos
    fseek(fid,36,'bof');
    
    pul = fread(fid, inf, '2*uint32', 36 ,'ieee-be');
    
    %% Makes an array with the list of floats
    numAtoms = length(pos)/9;
    
    pos = reshape(pos, [9 numAtoms]);
    pos = double(pos');
    ionIdx = (1:numAtoms)';
    
    
    pulse = reshape(pul, [2,numAtoms]);
    pulse = pulse';
    
    pos = table(ionIdx,pos(:,1),pos(:,2),pos(:,3),pos(:,4),pos(:,5),pos(:,6),pos(:,7),pos(:,8),pos(:,9),pulse(:,1),pulse(:,2),(1:numAtoms)');
    pos.Properties.VariableNames = {'ionIdx','x','y','z','mc','tof','VDC','VP','detx','dety','deltaP','multi','#atom'};
    pos.Properties.VariableUnits = {'1','nm','nm','nm','Da','ns','V','V','mm','mm','1','1','1'};    
    
end

% close file
fclose(fid);
