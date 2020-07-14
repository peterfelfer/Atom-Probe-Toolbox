function patchToObj(patch,objName,fileName)
% exports obj file into Wavefront obj
% if a list (vector) of patches is parsed it will be saved in one obj file 
% object names can be parsed in 'objNames'
%
% patchToObj();
%           opens a "Save *.obj file to" window if selected object is a
%           patch
% patchToObj(patch);
%           opens a "Save *.obj file to" window
% patchToObj(patch,[],fileName);
%           saves file in current folder
% patchToObj(patch,objName,fileName);
%           saves file in current folder
%
% INPUTS:
% patch: structure with faces (f) and vertices (v)
%
% objName: object names can be parsed in 'objNames'
%
% fileName: desired name of saved file with .obj as suffix

if ~exist('patch','var')         
    object = gco;
    patch.vertices = get(object,'vertices');
    patch.faces = get(object,'faces');
end
numPatch = length(patch); 

if ~exist('fileName','var')             % Create filepath
    [file path] = uiputfile('*.obj','Save *.obj file to');
    
    fileName = [path file];
end

if ~exist('objName','var')
    for p = 1:numPatch
        objName{p} = ['mesh' num2str(p)];
    end
elseif numPatch == 1
    objNameCell{1} = objName; % Create a cell array from a string
    objName = objNameCell;
end

fid = fopen(fileName,'wt');

if( fid == -1 )
    error('Cant open file.');
    return;
end

fprintf(fid, '# openVA obj exporter, (c) Peter Felfer, The University of Sydney 2012 \n');


%% writing each patch
offset = 0;

for p = 1:numPatch
    fprintf(fid, 'v %f %f %f\n', patch(p).vertices');
    patch(p).faces =  patch(p).faces + offset;
    offset = offset + length(patch(p).vertices(:,1));
    
end

%% writing faces for individual objects
for p = 1:numPatch
    oName = ['o ' objName{p} '\n'];
    fprintf(fid, oName);
    fprintf(fid, 'f %u %u %u\n', patch(p).faces');
end


fclose(fid);
clear fid;
end