function patchTo3mf(inputArg, fileName, varargin)
% patchTo3mf exports a patch (FV) structure or graphics object to a 3MF file.
%
% patchTo3mf();
%           opens a "Save *.3mf file to" window if selected object is a patch
% patchTo3mf(fvOrHandle);
%           opens a "Save *.3mf file to" window
% patchTo3mf(fvOrHandle, fileName);
%           saves file in current folder or specified path
%
% INPUTS:
% inputArg: patch handle, axes handle, or FV structure
% fileName: desired name of saved file with .3mf suffix
%
% name-value options:
%   'unit'   3MF unit: 'micrometer' (default), 'millimeter', 'meter'
%   'scale'  scale factor applied to vertices (default: 1)
%
% Behavior:
% - axes handle: exports all patch objects as separate 3MF objects
% - single patch handle: exports that mesh only
% - array of patch handles: combines into one object
% - FV array: combines into one object
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if ~exist('inputArg', 'var') || isempty(inputArg)
    inputArg = gco;
end

if exist('fileName', 'var') && (ischar(fileName) || isstring(fileName))
    if strcmpi(string(fileName), "unit") || strcmpi(string(fileName), "scale")
        varargin = [{fileName}, varargin];
        fileName = [];
    end
end

opts = struct('unit', 'micrometer', 'scale', 1);
opts = parseOptions(opts, varargin{:});

if ~exist('fileName', 'var') || isempty(fileName)
    [file, path] = uiputfile('*.3mf', 'Save *.3mf file to');
    if isequal(file, 0)
        return;
    end
    fileName = [path file];
end

if ~contains(fileName, '.')
    fileName = [fileName '.3mf'];
end

[meshes, exportMode] = resolveMeshes(inputArg);

if isempty(meshes)
    error('patchTo3mf:noMeshes', 'No patch objects found to export.');
end

if strcmp(exportMode, 'combine')
    meshes = combineMeshes(meshes);
end

% Apply scale
for i = 1:numel(meshes)
    meshes(i).vertices = meshes(i).vertices .* opts.scale;
end

write3mfFile(fileName, meshes, opts.unit);
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('patchTo3mf:invalidOptions', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k + 1};
        switch key
            case "unit"
                opts.unit = char(val);
            case "scale"
                opts.scale = double(val);
            otherwise
                error('patchTo3mf:invalidOption', 'Unknown option "%s".', key);
        end
    end
end

function [meshes, exportMode] = resolveMeshes(inputArg)
    meshes = struct('vertices', {}, 'faces', {}, 'colors', {});
    exportMode = 'single';

    if isstruct(inputArg)
        if numel(inputArg) > 1
            exportMode = 'combine';
        end
        for i = 1:numel(inputArg)
            meshes(end+1) = meshFromStruct(inputArg(i), []); %#ok<AGROW>
        end
        return;
    end

    if isgraphics(inputArg, 'axes')
        exportMode = 'separate';
        patches = findobj(inputArg, 'Type', 'patch');
        patches = flipud(patches);
        for i = 1:numel(patches)
            meshes(end+1) = meshFromPatch(patches(i)); %#ok<AGROW>
        end
        surfaces = findobj(inputArg, 'Type', 'surface');
        surfaces = flipud(surfaces);
        for i = 1:numel(surfaces)
            meshes(end+1) = meshFromSurface(surfaces(i)); %#ok<AGROW>
        end
        return;
    end

    if isgraphics(inputArg, 'patch')
        if numel(inputArg) > 1
            exportMode = 'combine';
        else
            exportMode = 'single';
        end
        for i = 1:numel(inputArg)
            meshes(end+1) = meshFromPatch(inputArg(i)); %#ok<AGROW>
        end
        return;
    end

    if isgraphics(inputArg, 'surface')
        if numel(inputArg) > 1
            exportMode = 'combine';
        else
            exportMode = 'single';
        end
        for i = 1:numel(inputArg)
            meshes(end+1) = meshFromSurface(inputArg(i)); %#ok<AGROW>
        end
        return;
    end

    if isgraphics(inputArg, 'figure')
        exportMode = 'separate';
        patches = findobj(inputArg, 'Type', 'patch');
        patches = flipud(patches);
        for i = 1:numel(patches)
            meshes(end+1) = meshFromPatch(patches(i)); %#ok<AGROW>
        end
        surfaces = findobj(inputArg, 'Type', 'surface');
        surfaces = flipud(surfaces);
        for i = 1:numel(surfaces)
            meshes(end+1) = meshFromSurface(surfaces(i)); %#ok<AGROW>
        end
        return;
    end
end

function mesh = meshFromPatch(h)
    v = get(h, 'Vertices');
    f = get(h, 'Faces');
    cdata = [];
    if isprop(h, 'FaceVertexCData')
        cdata = get(h, 'FaceVertexCData');
    end
    faceColor = 'none';
    if isprop(h, 'FaceColor')
        faceColor = get(h, 'FaceColor');
    end
    ax = ancestor(h, 'axes');

    mesh = meshFromArrays(v, f, cdata, faceColor, ax);
end

function mesh = meshFromSurface(h)
    x = get(h, 'XData');
    y = get(h, 'YData');
    z = get(h, 'ZData');
    c = [];
    if isprop(h, 'CData')
        c = get(h, 'CData');
    end
    faceColor = 'none';
    if isprop(h, 'FaceColor')
        faceColor = get(h, 'FaceColor');
    end
    ax = ancestor(h, 'axes');

    if ~isempty(c)
        [f, v, cdata] = surf2patch(x, y, z, c, 'triangles');
    else
        [f, v] = surf2patch(x, y, z, 'triangles');
        cdata = [];
    end

    mesh = meshFromArrays(v, f, cdata, faceColor, ax);
end

function mesh = meshFromStruct(fv, ax)
    v = fv.vertices;
    f = fv.faces;
    cdata = [];
    faceColor = 'none';
    if isfield(fv, 'facevertexcdata')
        cdata = fv.facevertexcdata;
    elseif isfield(fv, 'FaceVertexCData')
        cdata = fv.FaceVertexCData;
    end
    if isfield(fv, 'facecolor')
        faceColor = fv.facecolor;
    end
    mesh = meshFromArrays(v, f, cdata, faceColor, ax);
end

function mesh = meshFromArrays(vertices, faces, cdata, faceColor, ax)
    [facesTri, faceMap] = triangulateFaces(faces);
    verticesOut = vertices;
    colors = [];

    if ~isempty(cdata)
        [colors, verticesOut, facesTri] = resolveColors(vertices, facesTri, faceMap, cdata, faceColor, ax);
    elseif isnumeric(faceColor) && numel(faceColor) == 3
        colors = repmat(reshape(faceColor, 1, 3), size(vertices, 1), 1);
    end

    colors = normalizeColors(colors);

    mesh = struct('vertices', verticesOut, 'faces', facesTri, 'colors', colors);
end

function [colors, verticesOut, facesOut] = resolveColors(vertices, facesTri, faceMap, cdata, faceColor, ax)
    nVerts = size(vertices, 1);
    nFaces = size(facesTri, 1);
    if isempty(faceMap)
        origFaceCount = 0;
    else
        origFaceCount = max(faceMap);
    end
    verticesOut = vertices;
    facesOut = facesTri;
    colors = [];

    if size(cdata, 2) == 3
        if size(cdata, 1) == nVerts
            colors = cdata;
            return;
        elseif size(cdata, 1) == origFaceCount
            faceColors = cdata(faceMap, :);
            [verticesOut, facesOut, colors] = expandFaceColors(vertices, facesTri, faceColors);
            return;
        elseif size(cdata, 1) == nFaces
            [verticesOut, facesOut, colors] = expandFaceColors(vertices, facesTri, cdata);
            return;
        end
    end

    if isvector(cdata) || size(cdata, 2) == 1
        vals = double(cdata(:));
        if numel(vals) == nVerts
            colors = scalarToRgb(vals, ax);
            return;
        elseif numel(vals) == origFaceCount
            faceVals = vals(faceMap);
            faceColors = scalarToRgb(faceVals, ax);
            [verticesOut, facesOut, colors] = expandFaceColors(vertices, facesTri, faceColors);
            return;
        elseif numel(vals) == nFaces
            faceColors = scalarToRgb(vals, ax);
            [verticesOut, facesOut, colors] = expandFaceColors(vertices, facesTri, faceColors);
            return;
        elseif isscalar(vals)
            colors = repmat(scalarToRgb(vals, ax), nVerts, 1);
            return;
        end
    end

    if isnumeric(faceColor) && numel(faceColor) == 3
        colors = repmat(reshape(faceColor, 1, 3), nVerts, 1);
    end
end

function rgb = scalarToRgb(values, ax)
    cmap = parula(256);
    clim = [];
    if ~isempty(ax) && isgraphics(ax, 'axes')
        try
            cmap = colormap(ax);
            clim = caxis(ax);
        catch
        end
    end
    values = double(values(:));
    if isempty(clim) || any(~isfinite(clim))
        clim = [min(values) max(values)];
    end
    if clim(1) == clim(2)
        clim(2) = clim(1) + 1;
    end
    values = min(max(values, clim(1)), clim(2));
    t = (values - clim(1)) ./ (clim(2) - clim(1));
    idx = 1 + t * (size(cmap, 1) - 1);
    rgb = interp1(1:size(cmap, 1), cmap, idx, 'linear');
end

function colors = normalizeColors(colors)
    if isempty(colors)
        return;
    end
    colors = double(colors);
    if max(colors(:)) > 1
        colors = colors / 255;
    end
    colors = min(max(colors, 0), 1);
end

function [facesTri, faceMap] = triangulateFaces(faces)
    facesTri = zeros(0, 3);
    faceMap = zeros(0, 1);
    if isempty(faces)
        return;
    end
    faces = faces(all(isfinite(faces), 2), :);
    faces = faces(all(faces > 0, 2), :);
    for i = 1:size(faces, 1)
        face = faces(i, :);
        face = face(face > 0);
        if numel(face) < 3
            continue;
        end
        if numel(face) == 3
            facesTri(end+1, :) = face; %#ok<AGROW>
            faceMap(end+1, 1) = i; %#ok<AGROW>
        else
            v1 = face(1);
            for k = 2:(numel(face) - 1)
                facesTri(end+1, :) = [v1, face(k), face(k+1)]; %#ok<AGROW>
                faceMap(end+1, 1) = i; %#ok<AGROW>
            end
        end
    end
end

function [verticesOut, facesOut, colorsOut] = expandFaceColors(vertices, facesTri, faceColors)
    numFaces = size(facesTri, 1);
    verticesOut = vertices(facesTri(:), :);
    facesOut = reshape(1:(numFaces * 3), 3, [])';
    colorsOut = repelem(faceColors, 3, 1);
end

function meshesOut = combineMeshes(meshes)
    if numel(meshes) <= 1
        meshesOut = meshes;
        return;
    end
    vertices = [];
    faces = [];
    colors = [];
    offset = 0;
    for i = 1:numel(meshes)
        v = meshes(i).vertices;
        f = meshes(i).faces;
        c = meshes(i).colors;
        faces = [faces; f + offset]; %#ok<AGROW>
        vertices = [vertices; v]; %#ok<AGROW>
        if ~isempty(c)
            colors = [colors; c]; %#ok<AGROW>
        elseif ~isempty(colors)
            colors = [colors; zeros(size(v, 1), 3)]; %#ok<AGROW>
        end
        offset = offset + size(v, 1);
    end
    meshesOut = struct('vertices', vertices, 'faces', faces, 'colors', colors);
end

function write3mfFile(fileName, meshes, unitName)
    tempDir = tempname;
    mkdir(tempDir);
    mkdir(fullfile(tempDir, '_rels'));
    mkdir(fullfile(tempDir, '3D'));

    contentTypesPath = fullfile(tempDir, '[Content_Types].xml');
    relsPath = fullfile(tempDir, '_rels', '.rels');
    modelPath = fullfile(tempDir, '3D', '3dmodel.model');

    writeContentTypes(contentTypesPath);
    writeRels(relsPath);
    writeModel(modelPath, meshes, unitName);

    zipBase = tempname;
    zip(zipBase, {'[Content_Types].xml', '_rels/.rels', '3D/3dmodel.model'}, tempDir);
    zipFile = [zipBase '.zip'];
    if exist(fileName, 'file')
        delete(fileName);
    end
    movefile(zipFile, fileName);
    rmdir(tempDir, 's');
end

function writeContentTypes(path)
    fid = fopen(path, 'wt');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">\n');
    fprintf(fid, '  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>\n');
    fprintf(fid, '  <Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>\n');
    fprintf(fid, '</Types>\n');
    fclose(fid);
end

function writeRels(path)
    fid = fopen(path, 'wt');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    fprintf(fid, '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n');
    fprintf(fid, '  <Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/>\n');
    fprintf(fid, '</Relationships>\n');
    fclose(fid);
end

function writeModel(path, meshes, unitName)
    hasColors = false;
    for i = 1:numel(meshes)
        if ~isempty(meshes(i).colors)
            hasColors = true;
            break;
        end
    end

    fid = fopen(path, 'wt');
    fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
    if hasColors
        fprintf(fid, '<model unit="%s" xml:lang="en-US" ', unitName);
        fprintf(fid, 'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02" ');
        fprintf(fid, 'xmlns:m="http://schemas.microsoft.com/3dmanufacturing/material/2015/02">\n');
    else
        fprintf(fid, '<model unit="%s" xml:lang="en-US" xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">\n', unitName);
    end

    fprintf(fid, '  <resources>\n');

    objIds = zeros(numel(meshes), 1);
    colorGroupId = 1;

    for i = 1:numel(meshes)
        objId = i;
        objIds(i) = objId;
        colors = meshes(i).colors;

        if ~isempty(colors)
            [colorGroupId, colorIndex] = writeColorGroup(fid, colors, colorGroupId);
            meshes(i).colorIndex = colorIndex;
            meshes(i).pid = colorGroupId - 1;
        else
            meshes(i).colorIndex = [];
            meshes(i).pid = [];
        end

        writeObject(fid, meshes(i), objId);
    end

    fprintf(fid, '  </resources>\n');
    fprintf(fid, '  <build>\n');
    for i = 1:numel(objIds)
        fprintf(fid, '    <item objectid="%d"/>\n', objIds(i));
    end
    fprintf(fid, '  </build>\n');
    fprintf(fid, '</model>\n');
    fclose(fid);
end

function [nextId, colorIndex] = writeColorGroup(fid, colors, groupId)
    colors8 = uint8(round(colors * 255));
    [uniqueColors, ~, idx] = unique(colors8, 'rows', 'stable');
    colorIndex = idx - 1;

    fprintf(fid, '    <m:colorgroup id="%d">\n', groupId);
    for i = 1:size(uniqueColors, 1)
        c = uniqueColors(i, :);
        fprintf(fid, '      <m:color color="#%02X%02X%02X"/>\n', c(1), c(2), c(3));
    end
    fprintf(fid, '    </m:colorgroup>\n');
    nextId = groupId + 1;
end

function writeObject(fid, mesh, objId)
    v = mesh.vertices;
    f = mesh.faces;
    hasColor = ~isempty(mesh.colors);

    fprintf(fid, '    <object id="%d" type="model">\n', objId);
    fprintf(fid, '      <mesh>\n');
    fprintf(fid, '        <vertices>\n');
    for i = 1:size(v, 1)
        fprintf(fid, '          <vertex x="%g" y="%g" z="%g"/>\n', v(i, 1), v(i, 2), v(i, 3));
    end
    fprintf(fid, '        </vertices>\n');
    fprintf(fid, '        <triangles>\n');

    if hasColor
        pid = mesh.pid;
        colorIndex = mesh.colorIndex;
        for i = 1:size(f, 1)
            v1 = f(i, 1) - 1;
            v2 = f(i, 2) - 1;
            v3 = f(i, 3) - 1;
            p1 = colorIndex(f(i, 1));
            p2 = colorIndex(f(i, 2));
            p3 = colorIndex(f(i, 3));
            fprintf(fid, '          <triangle v1="%d" v2="%d" v3="%d" pid="%d" p1="%d" p2="%d" p3="%d"/>\n', ...
                v1, v2, v3, pid, p1, p2, p3);
        end
    else
        for i = 1:size(f, 1)
            fprintf(fid, '          <triangle v1="%d" v2="%d" v3="%d"/>\n', ...
                f(i, 1) - 1, f(i, 2) - 1, f(i, 3) - 1);
        end
    end

    fprintf(fid, '        </triangles>\n');
    fprintf(fid, '      </mesh>\n');
    fprintf(fid, '    </object>\n');
end
