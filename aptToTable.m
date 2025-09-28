function [posIn, header] = aptToTable(fileName)
% posToTable returns data from a .pos or .epos file and returns a table with the
% entries and units from an atom probe measurement
%
% [posIn, header] = aptToTable(fileName)
% [posIn, header] = aptToTable()
%
% INPUT
% fileName: is the name of the .pos or .epos file, it is optional. 
%           A dialog box will pop up if no file name is given.
% 
% OUTPUT
% posIn: is the variable that contains the entire data from the atom probe for further analysis,
%       table
%
% header: struct containing the additional information
%
% hint: in case of an epos file as input, be aware to select (*.epos) as
%       displayed data type in the dialog box
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg
% but using the Paraprobe APT loader by Markus Kühbach, MPIE
%%

if ~exist('fileName','var')
    [file, path] = uigetfile({'*.apt'},'Select a pos file');
    fileName = [path file];
    disp(['file ' file ' loading']);
end


%using paraprobe transcoder to load file
apt = PARAPROBE_Transcoder2(fileName); 
numAtoms = length(apt.Position(1,:));
ionIdx = (1:numAtoms)';

if isempty(apt.tof)
    apt.tof = zeros(size(apt.Mass));
end
if isempty(apt.Vref)
    apt.Vref = zeros(size(apt.Mass));
end
if isempty(apt.Vap)
    apt.Vap = zeros(size(apt.Mass));
end
if isempty(apt.XDet_mm)
    apt.XDet_mm = zeros(size(apt.Mass));
end
if isempty(apt.YDet_mm)
    apt.YDet_mm = zeros(size(apt.Mass));
end

%TODO: incorporate full epos functionality
posIn = table(ionIdx,apt.Position(1,:)',apt.Position(2,:)',apt.Position(3,:)',apt.Mass',apt.tof',apt.Vref',apt.Vap',apt.XDet_mm',apt.YDet_mm');%,pulse(:,1),pulse(:,2),(1:numAtoms)');
posIn.Properties.VariableNames = {'ionIdx','x','y','z','mc','tof','VDC','VP','detx','dety'};%,'deltaP','multi','atomNum'};
posIn.Properties.VariableUnits = {'1','nm','nm','nm','Da','ns','V','V','mm','mm'};%,'1','1','1'};    
end