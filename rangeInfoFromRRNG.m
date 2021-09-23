function [mcBegin, mcEnd, ionName, ionVolume, ionColor] = rangeInfoFromRRNG(rrngLine)
% rangeInfoFromRRNG interpetes a line from an RRNG file and extracts the
% information needed to create an ion and a range in a range figure

% constants and data preparation
HEXRANGE = 256;
rrngLine = char(rrngLine);
rrngLine = [rrngLine ' '];

if strcmp(rrngLine(1:5),'Range')
    % range limits
    rangeLimits = str2num(extractBefore(extractAfter(rrngLine,'='),'Vol'));
    mcBegin = min(rangeLimits);
    mcEnd = max(rangeLimits);
    
    % ion type - find match in ion list
    ionName = strtrim(extractAfter(extractAfter(extractBefore(rrngLine,'Color'),'Vol'),' '));
    ionName = ionName(ionName ~= ':');
    ionName = ionConvertName(ionConvertName(ionName)); %back and forth conversion for sanity checking anf formatting
    
    %ion volume
    ionVolume = str2num(extractBefore(rrngLine(strfind(rrngLine,'Vol')+4:end),' '));
    
    % ion color
    colHex = extractBefore(rrngLine(strfind(rrngLine,'Color')+6:end),' ');
    ionColor = [hex2dec(colHex(1:2)) hex2dec(colHex(3:4)) hex2dec(colHex(5:6))]/HEXRANGE;
else
    warning('parsed line is not an rrng range');
    mcBegin = [];
    mcEnd = [];
    ionName = [];
    ionVolume = [];
    ionColor = [];
end