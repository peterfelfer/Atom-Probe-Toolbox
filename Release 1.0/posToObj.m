function posToObj(pos,species,fileName)
% posToObj exports pos file's coordinates to an obj file, so that it can be
% read into various computer graphics programs such as Blender.
%
% posToObj(pos)                     exports obj file
% posToObj(pos,species)             exports obj file of selected species
% posToObj(pos,{})                  exports obj file of species selected dl
% posToObj(pos,species,fileName)    export to specific file name       
% posToObj(pos,~,fileName)          export whole variable to file name
% 
% INPUT
% pos =         pos is the pos variable to be exported
% 
% fileName =    name of the file
% 
% species =     optional with the format cellstr {'Fe','Cr',...} or 
%               string ["Fe", "Cr",..]
%               if species is an empty cell variable {}, an dialogue will
%               pop up with the possible selections
% 
% OUTPUT
% .obj file stored with the fileName
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

% get file name if it is not provided
if ~exist('fileName','var')
   getFileName = true; 
elseif isempty(fileName)
   getFileName = true;
end
if getFileName
    [file path] = uiputfile('*.obj','writing obj files');
    fileName = [path file(1:end-4)];
end


% determine which coordinates are to be written out
if ~exist('species','var')
    % all coordinates are used
    coords = [pos.x, pos.y, pos.z];
    
elseif iscell(species)
    % find if its ionic or atomic decomposed
    isAtomic = any(pos.Properties.VariableNames == "atom");
    
    if ~isempty(species)
        % select only the relevant atoms / ions
        if isAtomic
            isIn = ismember(pos.atom, species);
        else
            isIn = ismember(pos.ion, species);
        end
        
    else
        % select atoms / ions in a dialogue
        if isAtomic
            types = categories(pos.atom);
            [idx,tf] = listdlg('ListString',types);
            if tf
                isIn = ismember(pos.atom, types(idx));
            else 
                return
            end
        else % for ions
            types = categories(pos.ion);
            [idx,tf] = listdlg('ListString',types);
            if tf
                isIn = ismember(pos.ion, types(idx));
            else 
                return
            end
        end 
    end
    
    % create coordinate variable
    coords = [pos.x(isIn), pos.y(isIn), pos.z(isIn)];

else
    error('species format not recognized')
end


% write obj as text file
fid = fopen([fileName '.obj'],'wt');
if( fid == -1 )
    error('Cant open file.');
    return;
end
fprintf(fid, '# written by APT toolbox obj exporter, (c) Peter Felfer \n');
fprintf(fid, 'v %f %f %f\n',coords(:,1:3)');
fclose(fid);
clear fid;