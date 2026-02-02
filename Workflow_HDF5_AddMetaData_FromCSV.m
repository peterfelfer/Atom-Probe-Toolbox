%[text] # Workflow add metaData information from atom probe to HDF5 file
%[text] This workflow shows how to add meta Data information in a semi-automated way. Semi-automated because the user has to extract a .csv file from the LEAP 4000XHR software system with all important metaData for each dataset. The workflow takes this .csv file and adds the information to the already existing HDF5 files. To get to the right metaData information the code uses the unique identifier number for each atom probe measurement. 
%[text] 1. Get the metaData information from the .csv file
%[text] 2. Add additional information to the HDF5 files \
%%
%[text] 1\. Get the meta Data information from the .csv file
currentFolder = 'L:\existingHDF5files'; % Go to your database folder %[control:editfield:93aa]{"position":[17,39]}
searchAttributeType = '/atomProbeTomography/experiment/pulseType';  %[control:editfield:2a6b]{"position":[23,66]}
% add the file name of the .csv file
fileName = 'AllExperiments_202304.csv'; %[control:editfield:862f]{"position":[12,39]}
% import it to Matlab
opts = detectImportOptions(fileName); %[output:22b5684f]
opts = setvartype(opts,{'PulseFraction___', 'LaserPolarization', 'LaserPulseEnergy_nJ_'}, 'double');
opts = setvartype(opts,{'Result', 'Operator', 'LocalElectrode'}, 'categorical');
disp([opts.VariableNames' opts.VariableTypes'])
attrLEAP = readtable(fileName, opts);
clearvars opts

  %[control:button:2a91]{"position":[1,2]}
%%
%[text] 2\. Add additional information to the HDF5 files
% get the HDF5 file list
fileListH5 = dir('*.h5');

for k = 1:height(fileListH5)

    fileName = fileListH5(k).name;
    fileID = str2double(fileName(6:9));
    dataInfo = attrLEAP(attrLEAP.ExperimentID == fileID, :);

% writes attributes in existing HDF5 file
    h5writeatt(fileName,'/experiment','uniqueIdentifier', dataInfo.RunID);
    h5writeatt(fileName,'/specimen', 'name', dataInfo.Specimen);
    h5writeatt(fileName,'/project', 'identifier', double(string(dataInfo.Project)));
    h5writeatt(fileName,'/atomProbeTomography/experiment', 'specimenTemperatureActual', dataInfo.Temperature_K_);
    h5writeatt(fileName,'/experiment', 'instrumentOperator', double(string(dataInfo.Operator)));
    h5writeatt(fileName,'/atomProbeTomography/instrument','apertureUniqueIdentifier', double(string(dataInfo.LocalElectrode)));
    h5writeatt(fileName, '/atomProbeTomography/experiment', 'pulseFraction', double(string(dataInfo.PulseFraction___)));
    h5writeatt(fileName, '/atomProbeTomography/instrument', 'flightPathLength', dataInfo.FlightPathLength_mm_);
    h5writeatt(fileName, '/atomProbeTomography/experiment', 'laserPulseEnergy', dataInfo.LaserPulseEnergy_nJ_);
    h5writeatt(fileName, '/atomProbeTomography/experiment', 'pulsePeriod', 1/dataInfo.PulseFrequency_Hz_*10^9);
    h5writeatt(fileName, '/atomProbeTomography/experiment', 'pulseFraction', double(string(dataInfo.PulseFraction___)));
    
    
% write if it is laser or voltage pulsed experiment based on if LaserPulseEnergynJ is NaN
    t = isnan(dataInfo.LaserPulseEnergy_nJ_);
    if t == 1
            h5writeatt(fileName, '/atomProbeTomography/experiment','pulseType', 'voltagePulsed');
        else
            h5writeatt(fileName, '/atomProbeTomography/experiment','pulseType', 'laserPulsed');
    end

end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":10.4}
%---
%[control:editfield:93aa]
%   data: {"defaultValue":"'L:\\existingHDF5files\\LaserExperiments'","label":"desti","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:2a6b]
%   data: {"defaultValue":"'\/atomProbeTomography\/experiment\/pulseType'","label":"desti","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:862f]
%   data: {"defaultValue":"'AllExperiments_202304.csv'","label":"fileName","run":"Section","valueType":"Char"}
%---
%[control:button:2a91]
%   data: {"label":"Run","run":"Section"}
%---
%[output:22b5684f]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('detectImportOptions')\" style=\"font-weight:bold\">detectImportOptions<\/a>\nUnable to find or open 'AllExperiments_202304.csv'. Check the path and filename or file permissions."}}
%---
