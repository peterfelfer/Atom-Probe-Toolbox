function fv = stlToPatch(fileName, varargin)
% stlToPatch imports an STL file into a patch (FV) structure.
%
% fv = stlToPatch()
%           opens a "Select *.stl file" window
% fv = stlToPatch(fileName)
%           reads the specified STL file
% fv = stlToPatch(fileName, 'removeUnusedVertices', true)
%           removes duplicate/unused vertices and reindexes faces
%
% INPUTS:
% fileName:  STL file name
%
% name-value options:
%   'removeUnusedVertices'  true/false (default: false)
%
% OUTPUT:
% fv: structure with fields 'faces' and 'vertices'
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

removeUnusedVertices = false;
if ~isempty(varargin)
    removeUnusedVertices = parseOptions(varargin{:});
end

if ~exist('fileName', 'var') || isempty(fileName)
    [file, path] = uigetfile('*.stl', 'Select *.stl file');
    if isequal(file, 0)
        fv = struct('faces', [], 'vertices', []);
        return;
    end
    fileName = [path file];
end

if ~exist(fileName, 'file')
    error('stlToPatch:fileNotFound', 'File not found: %s', fileName);
end

fileInfo = dir(fileName);
fileSize = fileInfo.bytes;

fid = fopen(fileName, 'r', 'ieee-le');
if fid == -1
    error('stlToPatch:fileOpenError', 'Cannot open file.');
end

% Determine if binary STL based on file size
fseek(fid, 80, 'bof');
numTri = fread(fid, 1, 'uint32');
expectedSize = 84 + 50 * numTri;

fseek(fid, 0, 'bof');
if fileSize == expectedSize && numTri > 0
    [faces, vertices] = readBinaryStl(fid, numTri);
else
    fclose(fid);
    [faces, vertices] = readAsciiStl(fileName);
end

if removeUnusedVertices
    [vertices, faces] = mergeDuplicateVertices(vertices, faces);
end

fv = struct('faces', faces, 'vertices', vertices);
end

function [faces, vertices] = readBinaryStl(fid, numTri)
% Read binary STL file (already opened in ieee-le)

fseek(fid, 80, 'bof');
numTri = fread(fid, 1, 'uint32');
vertices = zeros(numTri * 3, 3);
faces = reshape(1:(numTri * 3), 3, [])';

for i = 1:numTri
    fread(fid, 3, 'float32'); % normal (ignored)
    v1 = fread(fid, 3, 'float32')';
    v2 = fread(fid, 3, 'float32')';
    v3 = fread(fid, 3, 'float32')';
    fread(fid, 1, 'uint16');  % attribute byte count

    base = (i - 1) * 3;
    vertices(base + 1, :) = v1;
    vertices(base + 2, :) = v2;
    vertices(base + 3, :) = v3;
end

fclose(fid);
end

function [faces, vertices] = readAsciiStl(fileName)
% Read ASCII STL file

txt = fileread(fileName);
tokens = regexp(txt, 'vertex\s+([-\d\.eE+]+)\s+([-\d\.eE+]+)\s+([-\d\.eE+]+)', 'tokens');

if isempty(tokens)
    error('stlToPatch:invalidStl', 'No vertices found in STL file.');
end

numVerts = numel(tokens);
vertices = zeros(numVerts, 3);
for i = 1:numVerts
    vertices(i, :) = str2double(tokens{i});
end

numTri = floor(numVerts / 3);
vertices = vertices(1:(numTri * 3), :);
faces = reshape(1:(numTri * 3), 3, [])';
end

function [verticesOut, facesOut] = mergeDuplicateVertices(vertices, faces)
% Merge duplicate vertices and reindex faces

[verticesOut, ~, idx] = unique(vertices, 'rows');
facesOut = idx(faces);
end

function removeUnusedVertices = parseOptions(varargin)
removeUnusedVertices = false;
if mod(numel(varargin), 2) ~= 0
    error('stlToPatch:invalidOptions', 'Options must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(string(varargin{k}));
    val = varargin{k + 1};
    if key == "removeunusedvertices"
        removeUnusedVertices = logical(val);
    else
        error('stlToPatch:invalidOption', 'Unknown option "%s".', key);
    end
end
end
