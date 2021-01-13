function rangeTable = rangesExtractFromFile(fileName)
% rangesExtractFromFile loads all ranges and additional information 
% from a MATLAB figure containing a mass spectrum 
%
% rangeTable = rangesExtractFromMassSpec(spec)
%
% INPUT
% fileName:     MATLAB figure containing an area plot that displays the mass spectrum (histogram of m/c
%               frequencies)either in raw counts or normalised to bin width
%               and total ion count
%
% OUTPUT 
% rangeTable:   table with allocated ranges of the ions and additional 
%               information(charge state, corresponding color code)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

if ~exist('fileName','var')
    [file, path] = uigetfile('*.fig','select figure file with ranges');
    fileName = [path '/' file];
end

% loading of ranges from figure
rangeFigureVisibility = 'invisible'; 
rngFig = openfig(fileName,rangeFigureVisibility);

spec = findobj( get(rngFig,'Children'), '-depth', 2, 'DisplayName', 'mass spectrum');

rangeTable = rangesExtractFromMassSpec(spec);

close(rngFig);