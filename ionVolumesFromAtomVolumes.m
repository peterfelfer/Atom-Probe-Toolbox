function [ionVolumeList, ionsOut] = ionVolumesFromAtomVolumes(ions,isotopeTable)
% ionVolumesFromAtomVolumes outputs a list of ion volumes based on the ions
% that are contained in the variable ions. 
% As input parameters table, categorical or an array of elements is
% possible.
%
% [ionVolumeList, ionsOut] = ionVolumesFromAtomVolumes(ions,isotopeTable);
%
% INPUT
%
% ions:             table, categorical or any other array that contains
%                   strings/charakters
%
% isotopeTable:     the parsed isotope table is the basis of the atom 
%                   density
%
% OUTPUT
% 
% ionVolumeList:    table where ionName and the depending ionVolume is
%                   stored
% ionsOut:          pos/ionList table with the correspondig ion Volume for
%                   each ion
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% input validation -> transfer everything in table format
    
if istable(ions)
    cat = categorical(categories(ions.ion));
elseif isstring(ions)
    ions = split(ions);
    cat = categorical(ions);
else
    cat = categorical(cellstr(ions));
end

ionsCat = table(cat, 'VariableNames', {'ionName'});

%% find AtomVolumes in isotopeTable and create ionVolumeList

ionVolume = zeros(height(ionsCat),1);
for i = 1:height(ionsCat)
    ionTableTest = ionConvertName(string(ionsCat.ionName(i)));
    k = 0;
    for j = 1:height(ionTableTest)
        volumeElement = mean(isotopeTable.atomDensity(isotopeTable.element == ionTableTest.element(j)));
        k = k+volumeElement;
    end
    ionVolume(i,1) = k;
end

ionVolumeList = addvars(ionsCat, ionVolume); 

%% Join both tables and maintain the order

[ionsOut, order] = outerjoin(ions,ionVolumeList,"Type","left","LeftKeys","ion",...
    "RightKeys","ionName",'RightVariables','ionVolume');
ionsOut(order,:) = ionsOut;

    
end
