function pos = posTableFromHDF5(fileName)
% posTableFromHDF5 extract the pos/epos table from an h5 file
%   Detailed explanation goes here


%% check for pos / epos file 


attval = h5readatt(fileName,'/atomProbeTomography/reconstruction','isDecomposed');


if attval == "false" 
    dataPath = '/atomProbeTomography/reconstruction/atom/';
else
    dataPath = '/atomProbeTomography/reconstruction/ion/';
end

dataPathDetectorEvent = '/atomProbeTomography/detectorEvent/';
% check for posType
    attrInfo = h5info(fileName,'/atomProbeTomography/detectorEvent'); 
    allDatasets = {attrInfo.Name}.';
    tof = strcmp('timeOfFlight', allDatasets);
    if sum(tof) == 1 
        posType = "epos";
    else 
        posType = "pos";
    end

%% get the data out of the h5 file 

% stored in /atomProbeTomography/reconstruction/atom/ or /ion
ionIdx = h5read(fileName, [dataPath 'fieldEvaporationSequenceIndex']);
x = h5read(fileName, [dataPath 'x']);
y = h5read(fileName, [dataPath 'y']);
z = h5read(fileName, [dataPath 'z']);
mc = h5read(fileName, [dataPath 'massToChargeState']); % end for pos file

if posType == "epos"
    % stored in the group '/atomProbeTomography/detectorEvent/'
    tof = h5read(fileName, [dataPathDetectorEvent 'timeOfFlight']);
    VDC = h5read(fileName, [dataPathDetectorEvent 'standingVoltage']);
    VP = h5read(fileName, [dataPathDetectorEvent 'voltagePulseAmplitude']);
    detx = h5read(fileName, [dataPathDetectorEvent 'x']);
    dety = h5read(fileName, [dataPathDetectorEvent 'y']);
    deltaP = h5read(fileName, [dataPathDetectorEvent 'pulseInterval']);
    %check if multi is stored in the HDF5 file
    multiCheck = strcmp('multi', allDatasets);
    if sum(multiCheck) == 1
        multi = h5read(fileName, [dataPathDetectorEvent 'multi']);
    end
    
    atomNumCheck = strcmp('atomNum', allDatasets);
    if sum(atomNumCheck) == 1
        atomNum = h5read(fileName, [dataPathDetectorEvent 'atomNum']);
    end
    % end for epos file
end

%for bug fixing
bug = length(x);

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


if posType == "pos"
    pos = table(ionIdx, x, y, z, mc); % complete initial pos file
else
    if sum(multiCheck) + sum(atomNumCheck)== 2
        pos = table(ionIdx, x, y, z, mc, tof, VDC, VP, detx, dety, deltaP, multi, atomNum); % complete initial epos file
    else
        pos = table(ionIdx, x, y, z, mc, tof, VDC, VP, detx, dety, deltaP);
    end
end


end
