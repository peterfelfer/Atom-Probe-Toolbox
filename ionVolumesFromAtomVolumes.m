function ionVolumeList = ionVolumesFromAtomVolumes(ions,isotopeTable)
% ionVolumesFromAtomVolumes outputs a list of ion volumes based on the ions
% that are contained in the variable ions. 
% As input parameters table, categorical or an array of elements is
% possible.
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
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
%% input validation
    
if istable(ions)
    ions = categorical(ions.ion(ions.ionComplexity >= 1));
else
    ions = categorical(ions);
end
cat = categorical(categories(ions));
ions = table(cat, 'VariableNames', {'ionName'});


%% find AtomVolumes in isotopeTable
for  i = 1:height(ions)
    ionTableTest = ionConvertName(string(ions.ionName(i)));
    k = 0;
    for j = 1:height(ionTableTest)
        volumeElement = mean(isotopeTable.atomDensity(isotopeTable.element == ionTableTest.element(j)));
        k = k+volumeElement;
    end
    ionVolume(i,1) = k; 
end

ionVolumeList = addvars(ions, ionVolume);
end
