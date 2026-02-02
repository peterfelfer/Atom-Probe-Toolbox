%[text] # Workflow find measurement conditions for each measurement
%[text] This workflow shows how to extract data from the variable VacPresLEAP made by the Workflow\_Import\_VacPresData\_LEAP and add them for each individual measurement to a result table, for Laser measurements the variable is called peakLaser and for voltage measurements peakVol.
%[text] The Variable VacPresLEAP contains the vacuum and temperature levels depending on time. With the begin and end time of a measurement, the vacuum level and the temperature are extracted from the VacPresLEAP variable and a mean value for the measurement is calculated.
%[text] In Addition the Worfkflow shows how to add additional information for each measurement searched in the HDF5 files for example pulseFraction, pulseFrequency or pulseEnergy for Laser Measurements.
%[text] If it is the same code - for both measurement types - you can choose at the beginning, which variable you want to use.
%[text] 
%[text] 1. Load Data
%[text] 2. Find analysis vacuum (VacPresLEAP & HDF5 database)
%[text] 3. Find temperature (VacPresLEAP & HDF5 database)
%[text] 4. Find pulseFraction, pulsePeriod - voltage measurement (HDF5 database)
%[text] 5. Find pulsePeriod, pulseEnergy - Laser measurement (HDF5 database)
%[text] 6. Save variables in workspace \
%%
%[text] ## 1. Load data
%[text] peakDataAll is calculated by the Workflow\_Calculate\_H.mlx and calculates for each measurement the concentration of the H peak - for Laser experiments the variable is called peakLaser, for Voltage experiments peakVoltage
%[text] VacPresLEAP is calculated by the Workflow\_Import\_VacPresData\_LEAP and has the vacuum level over time stored. 
load('L:\existingHDF5files\Laserexperiments\peakDataAll_H_New.mat');
peakLaser = peakDataAll;
load('L:\existingHDF5files\VoltageExperiments\peakDataAll_H_Voltage_All.mat');
peakVoltage = peakDataAll;
load('L:\existingHDF5files\VacPressures_CryoTemps\2016_2023_VacPres_LEAP.mat');
clearvars peakDataAll
%%
%[text] ## 2. Find analysis vacuum for each measurement and add it to the HDF5 file
%%%%%%%%%%%%%%%%%%%% voltage or laser?
peakDataAll = peakLaser;
% peakDataAll = peakVoltage;


% create Vector to store the vacuum level
vacAnalysis = zeros(height(peakDataAll),1);

% Go through every measurement and find the vacuum
for i = 1:height(peakDataAll)

    % Find start and end point for each measurement
    attrInfo = h5info(peakDataAll.name{i});
    beginMeas = h5readatt(peakDataAll.name{i}, '/experiment', 'beginDateTime');
    endMeas = h5readatt(peakDataAll.name{i}, '/experiment', 'endDateTime');
    beginDate = datetime(beginMeas, "Format","uuuu-MM-dd HH:mm:ss.SSS");
    endDate = datetime(endMeas, "Format","uuuu-MM-dd HH:mm:ss.SSS");

    %Find the vacuum during the measurement
    idxBeg = VacPresLEAP.dateTime>=beginDate;
    idxEnd = VacPresLEAP.dateTime<=endDate;
    idx = logical(idxBeg + idxEnd -1);
    VacPresMeas = VacPresLEAP(idx,:);

    % Calculate Mean value in the analysis chamber during measurement
    vacLev = mean(VacPresMeas.AnalysisUHV);

    % Write Vacuumlevel into the HDF5 file
    h5writeatt(peakDataAll.name{i}, '/atomProbeTomography/experiment','measurementChamberPressure', vacLev);
    vacAnalysis(i,1) = vacLev;
end
% Add Vacuum Level to peakDataAll matrix
peakDataAll = addvars(peakDataAll, vacAnalysis, 'NewVariableNames','AnalysisUHV' );

%%%%%%%%%%%%%%%%%%%% voltage or laser?
peakLaser = peakDataAll;
% peakVoltage = peakDataAll;
%%
%[text] ## 3. Find actual specimen temperature
%[text] The excel sheet gives as single value the specimen temperature actual for each measurement. With the import of the different vacuum level - the specimen temperature during the entire time is also accessible
%%%%%%%%%%%%%%%%%%%% voltage or laser?
% peakDataAll = peakLaser;
peakDataAll = peakVoltage; 

tempAnalysis = zeros(height(peakDataAll),1);

 for k = 1:height(peakDataAll)

    attrInfo = h5info(peakDataAll.name{k});
    beginMeas = h5readatt(peakDataAll.name{k}, '/experiment', 'beginDateTime');
    endMeas = h5readatt(peakDataAll.name{k}, '/experiment', 'endDateTime');

    beginDate = datetime(beginMeas, "Format","uuuu-MM-dd HH:mm:ss.SSS");
    endDate = datetime(endMeas, "Format","uuuu-MM-dd HH:mm:ss.SSS");


    idxBeg = VacPresLEAP.dateTime>=beginDate;
    idxEnd = VacPresLEAP.dateTime<=endDate;
    idx = logical(idxBeg + idxEnd -1);
    VacPresMeas = VacPresLEAP(idx,:);

    tempAnalysisMeas = mean(VacPresMeas.CryoTempChB);
    tempAnalysis(k,1) = tempAnalysisMeas;
 end

 peakDataAll = addvars(peakDataAll, tempAnalysis, 'NewVariableNames','SpecTemp' );

%%%%%%%%%%%%%%%%%%%% voltage or laser?
% peakLaser = peakDataAll;
peakVoltage = peakDataAll;
 
%%
%[text] ## 4. Find pulseFraction and pulsePeriod - Voltage measurement
pulsePeriod = zeros(height(peakVoltage),1);
pulseFraction = zeros(height(peakVoltage),1);

for l= 1:height(peakVoltage)
pulsePeriod(l, 1) = 1/h5readatt(peakVoltage.name{l}, '/atomProbeTomography/experiment', 'pulsePeriod')*10^6;
pulseFraction(l, 1) = h5readatt(peakVoltage.name{l}, '/atomProbeTomography/experiment', 'pulseFraction');
end

peakVoltage = addvars(peakVoltage, pulseFraction, pulsePeriod, 'NewVariableNames',{'pulseFraction', 'pulsePeriod'} );
%%
%[text] ## 5. Find pulsePeriod and pulseEnergy - Laser measurement
pulsePeriod = zeros(height(peakLaser),1);
pulseEnergy = zeros(height(peakLaser),1);

for l= 1:height(peakLaser)
pulsePeriod(l, 1) = 1/h5readatt(peakLaser.name{l}, '/atomProbeTomography/experiment', 'pulsePeriod')*10^6;
pulseEnergy(l, 1) = h5readatt(peakLaser.name{l}, '/atomProbeTomography/experiment', 'laserPulseEnergy');
end

peakLaser = addvars(peakLaser,pulsePeriod, pulseEnergy, 'NewVariableNames',{'pulsePeriod', 'pulseEnergy'});
%%
%[text] ## 6. Save variables in workspace
save("peakVoltage", "peakVoltage", '-mat')
save("peakLaser", "peakLaser", '-mat')

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":26.6}
%---
