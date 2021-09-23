function [h, txt, colorScheme] = rangeAddFromRRNG(spec, colorScheme, isotopeTable, rrng, ionList, useRrngColor)
% takes a range line from an RRNG file and interprets it. It creates a
% range and an ion if the ion does not yet exist in the plot. Ions are take
% from the ion list provided and compared to the ions in the range file. 
%
% 
% the output color scheme has the colors from the range file. Pre-existing
% colors are either overwritten or missing colors are added. the initial
% colorScheme can also be empty.

rrng = char(rrng);

% default is the use of the exisiting color scheme
if not(exist('useRrngColor','var'))
    useRrngColor = false;
end

% constants and data preparation
HEXRANGE = 256;
rrng = [rrng ' '];

% constants for ion fitting
IONFITTYPE = 'most abundant';



% interpret line
%ion volume TODO - GET INTO RANGE
vol = str2num(extractBefore(rrng(strfind(rrng,'Vol')+4:end),' '));

% ion color
colHex = extractBefore(rrng(strfind(rrng,'Color')+6:end),' ');
ionColor = [hex2dec(colHex(1:2)) hex2dec(colHex(3:4)) hex2dec(colHex(5:6))]/HEXRANGE;

% range limits
rangeLimits = str2num(extractBefore(extractAfter(rrng,'='),'Vol'));

% ion type - find match in ion list
ionStr = strtrim(extractAfter(extractAfter(extractBefore(rrng,'Color'),'Vol'),' '));
ionStr = ionStr(ionStr ~= ':');
ionStr = ionConvertName(ionConvertName(ionStr)); %back and forth conversion for sanity checking

% checking if an ion from the ion list can be identified in the range
potentialIons = ionList(ionList.ion == ionStr,:);



if isempty(potentialIons)
    disp(['ion ' ionStr ' is not on ion list. pls expand list to include ion']);
    isIn = 0;
else
    isIn = potentialIons.mc > rangeLimits(1) & potentialIons.mc < rangeLimits(2);
    potentialIons = potentialIons(isIn,:);
    chargeState = nnz(char(potentialIons.ionIsotopic(1)) == '+');
end


% adding color in the color scheme
% if the ion does not exist in the color scheme, it is always added
if not(isempty(potentialIons))
    ionName = char(potentialIons.ion(1));
    if not(any(colorScheme.ion == ionName))
        colorScheme.ion(end+1) = ionName;
        colorScheme.color(end,:) = ionColor;
        
        % if it exists, override needs to be activated
    elseif any(colorScheme.ion == potentialIons.ion(1)) & useRrngColor
        colorScheme.color(colorScheme.ion == potentialIons.ion(1),:) = ionColor;
    end
end


% adding range
if not(any(isIn))
    warning('no ion found for range - nothing added');
    h = [];
    txt = [];
    addIon = false;
elseif nnz(isIn) > 1
    warning('multiple ions found in range - only ion added');
    h = [];
    txt = [];
    addIon = true;
else
    ion = char(potentialIons.ionIsotopic);
    
    [h, txt] = rangeAdd(spec,colorScheme,ion,[rangeLimits(1) rangeLimits(2)]);
    h.UserData.ionVolume = vol;
    h.UserData.ionVolumeUnit = 'nm3';
    addIon = true;
end



% adding ion
if addIon
    % check if ion already exists
    il = ionsExtractFromMassSpec(spec);
    if isempty(il)
        ionAdd(spec,char(potentialIons.ion(1)),chargeState,isotopeTable,colorScheme,0,0,IONFITTYPE,0.1);
        
    elseif not(any(il.chargeState == chargeState & il.ionName == char(potentialIons.ion(1))))
        ionAdd(spec,char(potentialIons.ion(1)),chargeState,isotopeTable,colorScheme,0,0,IONFITTYPE,0.1);
    end
end




