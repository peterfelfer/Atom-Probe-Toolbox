%[text] # Workflow Load RRNG file into your mass Spectra
%[text] 
%[text] With the following workflow you can load a .rrng file into your mass Spectra. 
%[text] Please run the code in your command window.
%[text] 
% Load in the pos table and create the mass spectra
%pos = posToTable();
%spec = massSpecPlot(posIn,.01);

% Load the .rrng file, colorScheme and the isotopeTable
[file path] = uigetfile('*.rrng');
% load('colorScheme.mat')
% load('isotopeTable_naturalAbundances.mat')

% extract the Elements from the .rrng file and creat an ionList
elements = elementsExtractFromText(fileread([path file]));
ionList = ionsCreateComplex(elements,[1 2],isotopeTable,[1 2 3]); %[output:4ccf19a1]

%split the .rrng file into strings
rrng = string(splitlines(fileread([path file])));
rrng(rrng == "") = [];
useRrngColor = true;
%%
%[text] Create for a range information from the .rrng file and at it to the spectrum 
for i=1:length(rrng)
    [mcBegin, mcEnd, ionName, ionVolume, ionColor] = rangeInfoFromRRNG(rrng(i));
    if not(isempty(mcBegin))
    [h, txt, colorScheme] = rangeAddFromRangeInfo(spec, colorScheme, isotopeTable, ionList, mcBegin, mcEnd, ionName, ionVolume, ionColor, useRrngColor);
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":23.6}
%---
%[output:4ccf19a1]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'isotopeTable'."}}
%---
