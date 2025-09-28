function rangesToRRNG(rangeTable,fileName)
% rangesToRRNG save the current rangeTable as .rrng file in the location
% specified by fileName
%   
% INPUT:
%   rangeTable = rangeTable extracted from a mass spectra
%   fileName = fileName as a string with storage location if needed e.g.
%              "Z:\Test\Test.rrng"
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


%% first Part Ions
expIons = unique(rangeTable.rangeName);
numIons = length(expIons);
rrngText = "[Ions]" + newline + "Number=" + numIons;
% enumeration of ions
for i = 1:numIons
    rrngText = rrngText + newline + "Ion" + i + "=" + string(expIons(i));
end
%% second part ranges

numRange = length(rangeTable.rangeName);
hexColor = strings(1,numRange)';

% create color in hex
for k = 1 : numRange
    hexStr = "";
    for c = 1:3
        hexSingle = dec2hex(round(255*rangeTable.color(k,c)));
        if hexSingle == '0'
            hexSingle = string(hexSingle) + string(hexSingle);
        end
        hexStr = hexStr + string(hexSingle);
        hexColor(k,1) = hexStr;
    end
end

% create ranges part
rrngText = rrngText + newline + "[Ranges]" + newline + "Number=" + numRange;
for j = 1:numRange
    rrngText = rrngText + newline + "Range" + j + "=" + string(rangeTable.mcbegin(j)) + " " + string(rangeTable.mcend(j))...
        + " " + "Vol:" + string(rangeTable.volume(j)) ...
        + " " + string(rangeTable.rangeName(j))+ ":" + string(rangeTable.chargeState(j)) ...
        + " Color:" + hexColor(j);
end

%% third part save as .rrng file

if ~exist('fileName','var')             
    % Create filepath
    [file path] = uiputfile('*.rrng','Save *.rrng file to');
    fileName = [path file];
end

% creates the file
fid = fopen(fileName,'wt');
fprintf(fid, rrngText);



end

