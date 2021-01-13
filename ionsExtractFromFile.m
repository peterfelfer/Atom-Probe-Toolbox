function ionTable = ionsExtractFromFile(fileName)
% ionsExtractFromMassSpec pulls all ions and corresponding information 
% from a mass spectrum plot and gets all plots connected to the mass 
% spectrum
%
% ionTable = ionsExtractFromMassSpec(spec)
%
% INPUT
% fileName: MATLAB figure file containing area plot that displays the mass 
%           spectrum (histogram of m/c frequencies)either in raw counts or 
%           normalised to bin width and total ion count
%       
%
% OUTPUT 
% ionTable: table with allocated ions and additional information
%           (charge state, corresponding color code)
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

ionTable = ionsExtractFromMassSpec(spec);

close(rngFig);