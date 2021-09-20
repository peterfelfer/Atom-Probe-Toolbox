function pos = posTableFromHDF5(fileName)
% posTableFromHDF5 extract the pos/epos table from an h5 file
%   Detailed explanation goes here


%% check for pos / epos file 

attval = h5readatt(fileName,'/atomProbeTomography/reconstruction','isDecomposed');

if attval == 'true'
    dataPath = '/atomProbeTomography/reconstruction/atom/';
else
    dataPath = '/atomProbeTomography/reconstruction/ion/';
end

dataPathDetectorEvent = '/atomProbeTomography/detectorEvent/';

%% get the data out of the h5 file 

% stored in /atomProbeTomography/reconstruction/atom/
ionIdx = h5read(fileName, [dataPath 'fieldEvaporationSequenceIndex']);
x = h5read(fileName, [dataPath 'x']);
y = h5read(fileName, [dataPath 'y']);
z = h5read(fileName, [dataPath 'z']);
mc = h5read(fileName, [dataPath 'massToChargeState']);

%for bug fixing
bug = length(x);

% stored in the group '/atomProbeTomography/detectorEvent/'
tof = h5read(fileName, [dataPathDetectorEvent 'timeOfFlight']);
VDC = h5read(fileName, [dataPathDetectorEvent 'standingVoltage']);
VP = h5read(fileName, [dataPathDetectorEvent 'voltagePulseAmplitude']);
detx = h5read(fileName, [dataPathDetectorEvent 'x']);
dety = h5read(fileName, [dataPathDetectorEvent 'y']);
deltaP = h5read(fileName, [dataPathDetectorEvent 'pulseInterval']);

multi = (1 : bug)'; % not yet stored in the h5 files
atomNum = (1 : bug)';

% h5write does not support writing strings therefore just the pos file is
% stored in an h5 data set
% ion = (1 : bug)';
% chargeState = h5read(fileName, [dataPath 'chargeState']);
% 
% if attval == 'true' 
%     atom = (1 : bug)';
%     isotope = h5read(fileName, [dataPath 'isotope']);
%     ionComplexity = h5read(fileName, [dataPath 'ionComplexity']);
% end


if attval == true
    pos = table(ionIdx, x, y, z, mc, tof, VDC, VP, detx, dety, deltaP, multi, atomNum, ion, chargeState, atom, isotope, ionComplexity); % complete decomposed epos file
else 
    pos = table(ionIdx, x, y, z, mc, tof, VDC, VP, detx, dety, deltaP, multi, atomNum, ion, chargeState); % complete raw epos file
end


end
