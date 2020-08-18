function hdf5posTableAdd(fileName,pos)
% hdf5posTableAdd enables the user to add pos data 
% to an existing HDF5 file with the appropriate groups structure.
%
% hdf5posTableAdd(fileName,pos)
%
% INPUT 
% fileName:         full file name including path as string or char array
%
% pos:              pos table variable. Should contain the following entries:
%                   ionIdx,x,y,z,mc,tof,VDC,VP,detx,dety,deltaP,multi,ion,atom,isotope,chargeState,ionComplexity
%                   decomposition state is automatically determined and the
%                   corresponding data is written either into the group:
%                   '/atomProbeTomography/reconstruction/ion/'
%                   (undecomposed)
%                   '/atomProbeTomography/reconstruction/atom/'
%                   (decomposed)


%% writing pos data ==> APT specific!
numEntries = height(pos);
isDecomposed = numEntries > max(pos.ionIdx); % check if pos file is decomposed

if isDecomposed
    dataPath = '/atomProbeTomography/reconstruction/atom/';
else
    dataPath = '/atomProbeTomography/reconstruction/ion/';
end

% hdf5 write for individual variables
posColumnNames = pos.Properties.VariableNames;
for col = 1:width(pos)    
    if ismember(posColumnNames{col},{'x','y','z'}) %float type coords
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath posColumnNames{col}],[numEntries 1]);
        h5write(fileName,[dataPath posColumnNames{col}],data);
        h5writeatt(fileName,[dataPath posColumnNames{col}],'unit','nm','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'chargeState','isotope','ionComplexity'}) % int type coords
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath posColumnNames{col}],[numEntries 1]);
        h5write(fileName,[dataPath posColumnNames{col}],data);
        
    elseif ismember(posColumnNames{col},{'mc'}) % mass to charge state in Da
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath 'massToChargeState'],[numEntries 1]);
        h5write(fileName,[dataPath 'massToChargeState'],data);
        h5writeatt(fileName,[dataPath 'massToChargeState'],'unit','Da','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'ionIdx'}) % field evaporation index in the sequence
        data = table2array(pos(:,col));
        h5create(fileName,[dataPath 'fieldEvaporationSequenceIndex'],[numEntries 1]);
        h5write(fileName,[dataPath 'fieldEvaporationSequenceIndex'],data);
        
    elseif ismember(posColumnNames{col},{'ion','atom'}) % categorical arrays converted to string arrays
        %data = string(table2array(pos(:,col)));
        %dataType = 'string';
        data = [1];
        dataType = 'single';
        
        
        % epos variables
    elseif ismember(posColumnNames{col},{'tof'}) % export of ion flight times
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/timeOfFlight/','unit','ns','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'VDC'}) % export of experiment standing voltage
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/standingVoltage',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/standingVoltage',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/standingVoltage','unit','V','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'VP'}) % export of experiment pulse voltage
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/voltagePulseAmplitude/','unit','V','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'detx'}) % export of detector x hit coordinates
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/x/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/x/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/x/','unit','mm','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'dety'})% export of detector y hit coordinates
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/y/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/y/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/y/','unit','mm','TextEncoding','UTF-8');
        
    elseif ismember(posColumnNames{col},{'deltaP'})% export of detector pulse trigger intervals
        data = table2array(pos(:,col));
        h5create(fileName,'/atomProbeTomography/detectorEvent/pulseInterval/',[numEntries 1]);
        h5write(fileName,'/atomProbeTomography/detectorEvent/pulseInterval/',data);
        h5writeatt(fileName,'/atomProbeTomography/detectorEvent/pulseInterval/','unit','1','TextEncoding','UTF-8');
        
    end

end