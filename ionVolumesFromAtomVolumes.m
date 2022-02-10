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
    ions = categorical(ions.ion);
elseif isstring(ions)
    ions = split(ions);
    ions = categorical(ions);
else
    ions = categorical(ions);
end
cat = categorical(categories(ions));
ions = table(cat, 'VariableNames', {'ionName'});


%% find AtomVolumes in isotopeTable
% include waitbar since this could take a while
wb = waitbar(0,'building ion volumes');

for  i = 1:height(ions)
    ionTableTest = ionConvertName(string(ions.ionName(i)));
    k = 0;
    for j = 1:height(ionTableTest)
        volumeElement = mean(isotopeTable.atomDensity(isotopeTable.element == ionTableTest.element(j)));
        k = k+volumeElement;
    end
    ionVolume(i,1) = k;
    waitbar(i/length(ions),wb)
end

ionVolumeList = addvars(ions, ionVolume); % Ionen Liste mit gewichten 

%% attach ionVolumes to pos
% for  i = 1:height(ions) % Laenge von der Liste - gehe jeden Eitnrag durch 
%     ion = ions.ionName (i); % suche dir den ersten eintrag
%     % k = 0; %
%     for j = 1:height(ionVolumeList) % Laenge ion Volume liste
%     ionVolumeForPos = ionVolumeList.ionVolume(ion == ionVolumeList.ionName (j)); % suche entsprechendes volumen 
%     ;
%     end
%     ionVolumeForPosList(i,1) = k;
% end
% 
% ions = addvars(ionVolumeForPosList);

%% Search inside the pos/ion table for the ion and add the Volume of the ion
h = height(ions);
volume = zeros(h,1); 
ions = addvars(ions, volume); % Add a zero Volume column to the pos file
% include waitbar since this could take a while
wb = waitbar(0,'attach ionvolumes to input file');

for i = 1:height(ionVolumeList)
    for j = 1:height(ions)
        if ionVolumeList.ionName(i) == ions.ion(j)
            ions.testVolumen(j) = ionVolumeList.ionVolume(i);
        end 
    end 
    waitbar(i/length(ionVolumeList),wb)
end


    
end
