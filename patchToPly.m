function patchToPly(fv,vertColors,fileName,comment)
% patchToPly saves a patch to a ply file (polygon file format). This format
% contains information about the coloring of individual vertices and can be
% read by many popular 3d programs. A file definition can be found on: http://paulbourke.net/dataformats/ply/
%
% patchToPly(fv);                
%       opens a "Save *.ply file to" window
%
% patchToPly(fv, vertColors);    
%       opens a "Save *.ply file to" window
%
% patchToPly(fv, vertColors, 'filename.ply'); 
%       saves filename.ply in current folder
% 
% patchToPly(fv, vertColors, 'filename.ply', 'comment');
%       saves filename.ply in current folder, including a comment
%
% INPUTS:
%   fv:          structure with faces (f) and vertices (v)
%   
%   vertColors:  represents color code for all vertices;
%                must be a n-by-3 array with n = number of vertices 
%                color code values must be all integers between 0 and 255
%                or all values <1
%   
%   fileName:    desired filename, input as character array with .ply suffix
%
%   comment:     optional, inserts comment about values used as limits into ply
%                input as character array
%
%
% 
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
%
%
if ~exist('fileName','var')
    [file path] = uiputfile('*.ply','Save *.ply file to');
    
    fileName = [path file];
    
    if ~file
        disp('no file selected');
        return
    end
end



if ~exist('vertColors','var') & ~exist('fv','var')
    object = gco;
    
    fv.vertices = get(object,'vertices');
    fv.faces = get(object,'faces');
    
    
    cmap = get(gcf,'Colormap');
    lim = caxis(gca);
    
    vals = get(object,'FaceVertexCdata');


    vals(vals > lim(2)) = lim(2);
    vals(vals < lim(1)) = lim(1);
    
    vertColors = value2VertColor(vals,cmap,lim);
    
    comment = ['value min: ' num2str(min(vals),2) ' at/nm2 '...
    'max: ' num2str(max(vals)) ' at/nm2']; % inserts comment about values used as limits into ply
    
    
end


numVerts = length(fv.vertices(:,1));
numFaces = length(fv.faces(:,1));

% indexing is 0 based in *.ply
fv.faces = fv.faces - 1;


%% writing header
header = ['ply \nformat ascii 1.0 \n'];
header = [header 'comment openVA ply exporter, (c) Peter Felfer, The University of Sydney 2013 \n'];

if exist('comment','var')
    header = [header comment '\n'];
end

% vertex properties

header = [header 'element vertex ' num2str(numVerts) ' \n'];
header = [header 'property float x \nproperty float y \nproperty float z \n'];

if exist('vertColors','var')
    header = [header 'property uchar red \nproperty uchar green \nproperty uchar blue \n'];
    
    % vertex colors are from 0 - 255
    
    if max(vertColors)<=1
        vertColors = round(vertColors* 255);
    end
    if max(vertColors)>255
        error('vertColors values exceed 255')
    end
    
    
end

header = [header 'element face ' num2str(numFaces) ' \n'];
header = [header 'property list uchar uint vertex_indices\n'];
header = [header 'end_header \n'];



%% file writing
if fileName
    fid = fopen(fileName,'wt');


    if( fid == -1 )
        error('Cant open file.');
        return;
    end
    
else
    return
    
end

fprintf(fid, header);


if exist('vertColors','var')
    for v = 1:numVerts
        fprintf(fid, '%f %f %f ', fv.vertices(v,:)');
        fprintf(fid, '%u %u %u\n', vertColors(v,:)');
    end
else
    
    fprintf(fid, '%f %f %f\n', fv.vertices');
end

fprintf(fid, '3 %u %u %u\n', fv.faces');


fclose(fid);
clear fid;
end


















