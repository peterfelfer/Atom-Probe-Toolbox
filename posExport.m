function posExport(posIn, fileName)
% posExport writes data from a pos table to a .pos or .epos file
%
% tableToPos(posIn, fileName)
% tableToPos(posIn)
%
% INPUT
% posIn: is the pos table variable that contains the data from the atom probe for writing to a file
% fileName: is the name of the output .pos or .epos file, it is optional. 
%           A dialog box will pop up if no file name is given.
%
% OUTPUT
% None
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
%%

if nargin < 2
    [file, path, ~] = uiputfile({'*.pos';'*.epos'},'Save as pos file');
    fileName = [path file];
    disp(['Saving to file ' file]);
end

[~, ~, ext] = fileparts(fileName);

switch ext
    case '.pos'
        idx = 1;
    case '.epos'
        idx = 2;
    otherwise
        error('Invalid file extension. Only .pos and .epos files are supported.');
end

fid = fopen(fileName, 'w');

if idx == 1
    posData = posIn{:, 2:5};
    posData = posData';
    fwrite(fid, posData(:), 'float32', 'b');
    
elseif idx == 2

    posData = posIn{:, 2:10};
    posData = posData';
    frewind(fid);
    fwrite(fid, posData, '9*float32',8, 'ieee-be');

    frewind(fid);
    pulseData = posIn{:, 11:12};
    pulseData = pulseData';
    fwrite(fid, pulseData, '2*uint32', 36, 'ieee-be');

end

fclose(fid);
end
