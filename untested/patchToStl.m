function patchToStl(fv, fileName, varargin)
% patchToStl exports a patch (FV) structure to an STL file (ASCII or binary).
%
% patchToStl()
%           opens a "Save *.stl file to" window if selected object is a patch
% patchToStl(fv)
%           opens a "Save *.stl file to" window
% patchToStl(fv, fileName)
%           saves file in current folder or specified path
% patchToStl(fv, fileName, 'format', 'binary')
%           exports binary STL
%
% INPUTS:
% fv:        structure with fields 'faces' and 'vertices'
% fileName:  desired name of saved file with .stl suffix
%
% name-value options:
%   'format'  'ascii' | 'binary' (default: 'ascii')
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if ~exist('fv', 'var') || isempty(fv)
    fv = gco;
end

if isgraphics(fv)
    fv = fvFromHandle(fv);
end

% Allow patchToStl(fv, 'format', 'binary')
if exist('fileName', 'var') && (ischar(fileName) || isstring(fileName))
    if strcmpi(string(fileName), "format")
        varargin = [{fileName}, varargin];
        fileName = [];
    end
end

format = 'ascii';
if ~isempty(varargin)
    if mod(numel(varargin), 2) ~= 0
        error('patchToStl:invalidOptions', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k + 1};
        switch key
            case "format"
                format = lower(char(val));
            otherwise
                error('patchToStl:invalidOption', 'Unknown option "%s".', key);
        end
    end
end

if ~exist('fileName', 'var') || isempty(fileName)
    [file, path] = uiputfile('*.stl', 'Save *.stl file to');
    if isequal(file, 0)
        return;
    end
    fileName = [path file];
end

if ~isfield(fv, 'vertices') || ~isfield(fv, 'faces')
    error('patchToStl:invalidInput', 'Input must have fields ''vertices'' and ''faces''.');
end

if isempty(fv.faces) || isempty(fv.vertices)
    error('patchToStl:emptyPatch', 'Faces or vertices are empty.');
end

if ~contains(fileName, '.')
    fileName = [fileName '.stl'];
end

faces = fv.faces;
vertices = fv.vertices;

faces = triangulateFaces(faces);
if isempty(faces)
    error('patchToStl:invalidFaces', 'No valid triangular faces found.');
end

switch lower(format)
    case 'ascii'
        fid = fopen(fileName, 'wt');
        if fid == -1
            error('patchToStl:fileOpenError', 'Cannot open file for writing.');
        end

        [~, solidName] = fileparts(fileName);
        fprintf(fid, 'solid %s\n', solidName);

        for i = 1:size(faces, 1)
            idx = faces(i, :);
            v1 = vertices(idx(1), :);
            v2 = vertices(idx(2), :);
            v3 = vertices(idx(3), :);

            n = cross(v2 - v1, v3 - v1);
            nNorm = norm(n);
            if nNorm > 0
                n = n ./ nNorm;
            else
                n = [0 0 0];
            end

            fprintf(fid, '  facet normal %g %g %g\n', n(1), n(2), n(3));
            fprintf(fid, '    outer loop\n');
            fprintf(fid, '      vertex %g %g %g\n', v1(1), v1(2), v1(3));
            fprintf(fid, '      vertex %g %g %g\n', v2(1), v2(2), v2(3));
            fprintf(fid, '      vertex %g %g %g\n', v3(1), v3(2), v3(3));
            fprintf(fid, '    endloop\n');
            fprintf(fid, '  endfacet\n');
        end

        fprintf(fid, 'endsolid %s\n', solidName);
        fclose(fid);

    case 'binary'
        fid = fopen(fileName, 'w', 'ieee-le');
        if fid == -1
            error('patchToStl:fileOpenError', 'Cannot open file for writing.');
        end

        header = uint8(zeros(1, 80));
        tag = uint8('Created by Atom Probe Toolbox');
        header(1:numel(tag)) = tag;
        fwrite(fid, header, 'uint8');
        fwrite(fid, uint32(size(faces, 1)), 'uint32');

        for i = 1:size(faces, 1)
            idx = faces(i, :);
            v1 = vertices(idx(1), :);
            v2 = vertices(idx(2), :);
            v3 = vertices(idx(3), :);

            n = cross(v2 - v1, v3 - v1);
            nNorm = norm(n);
            if nNorm > 0
                n = n ./ nNorm;
            else
                n = [0 0 0];
            end

            fwrite(fid, single(n), 'float32');
            fwrite(fid, single(v1), 'float32');
            fwrite(fid, single(v2), 'float32');
            fwrite(fid, single(v3), 'float32');
            fwrite(fid, uint16(0), 'uint16'); % attribute byte count
        end

        fclose(fid);
    otherwise
        error('patchToStl:invalidFormat', 'Format must be ''ascii'' or ''binary''.');
end
end

function facesTri = triangulateFaces(faces)
% Ensure faces are triangles (fan triangulation for polygons)

facesTri = zeros(0, 3);
if isempty(faces)
    return;
end

function fvOut = fvFromHandle(h)
% Convert patch/surface/axes/figure handles to a single FV struct

    meshes = struct('vertices', {}, 'faces', {});

    if isgraphics(h, 'axes')
        patches = findobj(h, 'Type', 'patch');
        surfaces = findobj(h, 'Type', 'surface');
        patches = flipud(patches);
        surfaces = flipud(surfaces);
        for i = 1:numel(patches)
            meshes(end+1) = fvFromPatch(patches(i)); %#ok<AGROW>
        end
        for i = 1:numel(surfaces)
            meshes(end+1) = fvFromSurface(surfaces(i)); %#ok<AGROW>
        end
    elseif isgraphics(h, 'figure')
        patches = findobj(h, 'Type', 'patch');
        surfaces = findobj(h, 'Type', 'surface');
        patches = flipud(patches);
        surfaces = flipud(surfaces);
        for i = 1:numel(patches)
            meshes(end+1) = fvFromPatch(patches(i)); %#ok<AGROW>
        end
        for i = 1:numel(surfaces)
            meshes(end+1) = fvFromSurface(surfaces(i)); %#ok<AGROW>
        end
    elseif isgraphics(h, 'surface')
        for i = 1:numel(h)
            meshes(end+1) = fvFromSurface(h(i)); %#ok<AGROW>
        end
    elseif isgraphics(h, 'patch')
        for i = 1:numel(h)
            meshes(end+1) = fvFromPatch(h(i)); %#ok<AGROW>
        end
    end

    if isempty(meshes)
        fvOut = struct('vertices', [], 'faces', []);
        return;
    end

    if numel(meshes) == 1
        fvOut = meshes;
        return;
    end

    vertices = [];
    faces = [];
    offset = 0;
    for i = 1:numel(meshes)
        v = meshes(i).vertices;
        f = meshes(i).faces;
        faces = [faces; f + offset]; %#ok<AGROW>
        vertices = [vertices; v]; %#ok<AGROW>
        offset = offset + size(v, 1);
    end
    fvOut = struct('vertices', vertices, 'faces', faces);
end

function fv = fvFromPatch(h)
    fv.vertices = get(h, 'vertices');
    fv.faces = get(h, 'faces');
end

function fv = fvFromSurface(h)
    x = get(h, 'XData');
    y = get(h, 'YData');
    z = get(h, 'ZData');
    [f, v] = surf2patch(x, y, z, 'triangles');
    fv.vertices = v;
    fv.faces = f;
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
    else
        v1 = face(1);
        for k = 2:(numel(face) - 1)
            facesTri(end+1, :) = [v1, face(k), face(k+1)]; %#ok<AGROW>
        end
    end
end
end
