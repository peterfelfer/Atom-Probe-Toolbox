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

if ~exist('fv','var') || isempty(fv)
    fv = gco;
end

vertColorsAuto = [];
commentAuto = '';
ax = [];

if isgraphics(fv)
    [fv, vertColorsAuto, commentAuto, ax] = fvFromHandleForPly(fv);
end

if (~exist('vertColors','var') || isempty(vertColors)) && ~isempty(vertColorsAuto)
    vertColors = vertColorsAuto;
end

if (~exist('comment','var') || isempty(comment)) && ~isempty(commentAuto)
    comment = commentAuto;
end

if isstruct(fv) && numel(fv) > 1
    [fv, vertColors] = combineFvAndColors(fv, vertColors);
end

if (~exist('vertColors','var') || isempty(vertColors)) && isstruct(fv) && isfield(fv, 'facevertexcdata')
    vertColors = fv.facevertexcdata;
end

if exist('vertColors','var') && ~isempty(vertColors)
    [fv, vertColors] = normalizeVertColors(fv, vertColors, ax);
end

if ~exist('fileName','var')
    [file path] = uiputfile('*.ply','Save *.ply file to');
    
    fileName = [path file];
    
    if ~file
        disp('no file selected');
        return
    end
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

function [fv, vertColors, comment, ax] = fvFromHandleForPly(h)
    vertColors = [];
    comment = '';
    ax = [];

    meshes = struct('vertices', {}, 'faces', {}, 'colors', {});
    comments = {};

    if isgraphics(h, 'axes')
        ax = h;
        patches = findobj(h, 'Type', 'patch');
        surfaces = findobj(h, 'Type', 'surface');
        patches = flipud(patches);
        surfaces = flipud(surfaces);
        for i = 1:numel(patches)
            [m, c, com] = meshFromPatch(patches(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
        for i = 1:numel(surfaces)
            [m, c, com] = meshFromSurface(surfaces(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
    elseif isgraphics(h, 'figure')
        patches = findobj(h, 'Type', 'patch');
        surfaces = findobj(h, 'Type', 'surface');
        patches = flipud(patches);
        surfaces = flipud(surfaces);
        for i = 1:numel(patches)
            [m, c, com] = meshFromPatch(patches(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
        for i = 1:numel(surfaces)
            [m, c, com] = meshFromSurface(surfaces(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
    elseif isgraphics(h, 'surface')
        ax = ancestor(h, 'axes');
        for i = 1:numel(h)
            [m, c, com] = meshFromSurface(h(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
    elseif isgraphics(h, 'patch')
        ax = ancestor(h, 'axes');
        for i = 1:numel(h)
            [m, c, com] = meshFromPatch(h(i));
            meshes(end+1) = struct('vertices', m.vertices, 'faces', m.faces, 'colors', c); %#ok<AGROW>
            if ~isempty(com)
                comments{end+1} = com; %#ok<AGROW>
            end
        end
    end

    if isempty(meshes)
        fv = struct('vertices', [], 'faces', []);
        return;
    end

    [fv, vertColors] = combineMeshes(meshes);

    if ~isempty(comments)
        comment = comments{1};
    end
end

function [mesh, vertColors, comment] = meshFromPatch(h)
    mesh.vertices = get(h, 'vertices');
    mesh.faces = get(h, 'faces');
    cdata = [];
    if isprop(h, 'FaceVertexCData')
        cdata = get(h, 'FaceVertexCData');
    end
    comment = '';
    if ~isempty(cdata)
        ax = ancestor(h, 'axes');
        [mesh, vertColors] = normalizeVertColors(mesh, cdata, ax);
        if isvector(cdata) || size(cdata, 2) == 1
            comment = buildComment(cdata);
        end
    else
        vertColors = [];
    end
end

function [mesh, vertColors, comment] = meshFromSurface(h)
    x = get(h, 'XData');
    y = get(h, 'YData');
    z = get(h, 'ZData');
    c = [];
    if isprop(h, 'CData')
        c = get(h, 'CData');
    end
    if ~isempty(c)
        [f, v, cdata] = surf2patch(x, y, z, c, 'triangles');
    else
        [f, v] = surf2patch(x, y, z, 'triangles');
        cdata = [];
    end
    mesh.vertices = v;
    mesh.faces = f;
    comment = '';
    if ~isempty(cdata)
        ax = ancestor(h, 'axes');
        [mesh, vertColors] = normalizeVertColors(mesh, cdata, ax);
        if isvector(cdata) || size(cdata, 2) == 1
            comment = buildComment(cdata);
        end
    else
        vertColors = [];
    end
end

function [meshOut, vertColorsOut] = normalizeVertColors(meshIn, vertColorsIn, ax)
    meshOut = meshIn;
    vertColorsOut = [];

    if isempty(vertColorsIn)
        return;
    end

    nVerts = size(meshIn.vertices, 1);
    nFaces = size(meshIn.faces, 1);

    if size(vertColorsIn, 2) == 3
        if size(vertColorsIn, 1) == nVerts
            vertColorsOut = vertColorsIn;
        elseif size(vertColorsIn, 1) == nFaces
            [meshOut, vertColorsOut] = expandFaceColors(meshIn, vertColorsIn);
        elseif size(vertColorsIn, 1) == 1
            vertColorsOut = repmat(vertColorsIn, nVerts, 1);
        end
    else
        vals = double(vertColorsIn(:));
        if numel(vals) == nVerts
            vertColorsOut = valuesToRgb(vals, ax);
        elseif numel(vals) == nFaces
            faceColors = valuesToRgb(vals, ax);
            [meshOut, vertColorsOut] = expandFaceColors(meshIn, faceColors);
        elseif numel(vals) == 1
            vertColorsOut = repmat(valuesToRgb(vals, ax), nVerts, 1);
        end
    end

    if ~isempty(vertColorsOut)
        if max(vertColorsOut(:)) > 1
            vertColorsOut = vertColorsOut / 255;
        end
        vertColorsOut = min(max(vertColorsOut, 0), 1);
    end
end

function rgb = valuesToRgb(values, ax)
    cmap = parula(256);
    lim = [];
    if ~isempty(ax) && isgraphics(ax, 'axes')
        try
            cmap = colormap(ax);
            lim = caxis(ax);
        catch
        end
    end
    values = double(values(:));
    if isempty(lim) || any(~isfinite(lim))
        lim = [min(values) max(values)];
    end
    if lim(1) == lim(2)
        lim(2) = lim(1) + 1;
    end
    values(values > lim(2)) = lim(2);
    values(values < lim(1)) = lim(1);
    rgb = value2VertColor(values, cmap, lim);
end

function comment = buildComment(vals)
    vals = double(vals(:));
    if isempty(vals)
        comment = '';
        return;
    end
    comment = ['value min: ' num2str(min(vals),2) ' at/nm2 '...
        'max: ' num2str(max(vals)) ' at/nm2'];
end

function [meshOut, vertColorsOut] = expandFaceColors(meshIn, faceColors)
    faces = meshIn.faces;
    numFaces = size(faces, 1);
    vOut = meshIn.vertices(faces(:), :);
    fOut = reshape(1:(numFaces * 3), 3, [])';
    cOut = repelem(faceColors, 3, 1);
    meshOut = struct('vertices', vOut, 'faces', fOut);
    vertColorsOut = cOut;
end

function [fvOut, vertColorsOut] = combineMeshes(meshes)
    vertices = [];
    faces = [];
    colors = [];
    offset = 0;
    for i = 1:numel(meshes)
        v = meshes(i).vertices;
        f = meshes(i).faces;
        faces = [faces; f + offset]; %#ok<AGROW>
        vertices = [vertices; v]; %#ok<AGROW>
        if ~isempty(meshes(i).colors)
            colors = [colors; meshes(i).colors]; %#ok<AGROW>
        elseif ~isempty(colors)
            colors = [colors; zeros(size(v, 1), 3)]; %#ok<AGROW>
        end
        offset = offset + size(v, 1);
    end
    fvOut = struct('vertices', vertices, 'faces', faces);
    vertColorsOut = colors;
end

function [fvOut, vertColorsOut] = combineFvAndColors(fvArray, vertColorsIn)
    vertices = [];
    faces = [];
    colors = [];
    offset = 0;
    totalVerts = 0;
    for i = 1:numel(fvArray)
        totalVerts = totalVerts + size(fvArray(i).vertices, 1);
    end
    hasInputColors = exist('vertColorsIn', 'var') && ~isempty(vertColorsIn);

    for i = 1:numel(fvArray)
        v = fvArray(i).vertices;
        f = fvArray(i).faces;
        faces = [faces; f + offset]; %#ok<AGROW>
        vertices = [vertices; v]; %#ok<AGROW>
        if ~hasInputColors
            if isfield(fvArray(i), 'facevertexcdata') && ~isempty(fvArray(i).facevertexcdata)
                colors = [colors; fvArray(i).facevertexcdata]; %#ok<AGROW>
            elseif ~isempty(colors)
                colors = [colors; zeros(size(v, 1), 3)]; %#ok<AGROW>
            end
        end
        offset = offset + size(v, 1);
    end

    fvOut = struct('vertices', vertices, 'faces', faces);

    if hasInputColors
        if size(vertColorsIn, 1) ~= totalVerts
            error('patchToPly:colorSizeMismatch', 'vertColors size does not match combined vertices.');
        end
        vertColorsOut = vertColorsIn;
    else
        vertColorsOut = colors;
    end
end
