function meta = metaDataReadLEAPAttributes(meta, exNum, exList)
%metaDataReadLEAPAttributes adds the extracted attributes of a LEAP
% measurement to the corresponding meta Data file
%
% INPUT
% meta = meta data List 
% exNum = specificNumber of Experiment
% exList = experiment List with the output Data from the LEAP 
%
% OUTPUT
% meta = meta Data File with the Attributes


% get important row out of exList
intRow = exList(exList.ExperimentID == exNum, :);

% sort Variables 
meta(contains(meta(:,1),'specimen/uniqueIdentifier'),2) = intRow.Specimen;
meta(contains(meta(:,1),'project/identifier'),2) = cellstr(char(intRow.Project));
meta(contains(meta(:,1),'experiment/specimenTemperatureActual'),2) = num2cell(intRow.Temperature_K_);
% start time
t = datetime(intRow.Date, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS');
date = datestr(t, 'yyyy-mm-dd');
time = datestr(t, 'HH:MM:SS.FFF');
meta(contains(meta(:,1),'experiment/beginDateTime'),2) = cellstr([date ' ' time]);
% end time
endDateTime = intRow.Date + minutes(intRow.RunTime_Hrs_*60);
endDate = datestr(endDateTime, 'yyyy-mm-dd');
endTime = datestr(endDateTime, 'HH:MM:SS.FFF');
meta(contains(meta(:,1),'experiment/endDateTime'),2) = cellstr([endDate ' ' endTime]);

meta(contains(meta(:,1),'experiment/result'),2) = {intRow.Result};
meta(contains(meta(:,1),'experiment/status'),2) = {'finished'};
meta(contains(meta(:,1),'experiment/uniqueIdentifier'),2) = intRow.RunID;
meta(contains(meta(:,1),'experiment/instrumentOperator'),2) = cellstr(char(intRow.Operator));
if ~isnan(intRow.PulseFraction___)
    meta(contains(meta(:,1),'experiment/pulseFraction'),2) = num2cell(intRow.PulseFraction___);
    meta(contains(meta(:,1),'experiment/pulseType'),2) = {'voltagePulsed'};
end
meta(contains(meta(:,1),'experiment/pulsePeriod'),2) = num2cell((1/intRow.PulseFrequency_Hz_)*10^9);
meta(contains(meta(:,1),'instrument/flightPathLength'),2) = num2cell(intRow.FlightPathLength_mm_);
if ~isnan(intRow.LaserPulseEnergy_nJ_)
    meta(contains(meta(:,1),'experiment/laserPulseEnergy'),2) = num2cell(intRow.LaserPulseEnergy_nJ_);
    meta(contains(meta(:,1),'experiment/pulseType'),2) = {'laserPulsed'};
end
meta(contains(meta(:,1),'experiment/notes'),2) = intRow.Comments;

end

