function fv = threeMfToPatch(fileName, varargin)
% threeMfToPatch imports a 3MF file into a patch (FV) structure.
%
% fv = threeMfToPatch()
%           opens a "Select *.3mf file" window
% fv = threeMfToPatch(fileName)
%           reads the specified 3MF file
% fv = threeMfToPatch(fileName, 'removeUnusedVertices', true)
%           removes duplicate/unused vertices and reindexes faces
%
% name-value options:
%   'removeUnusedVertices'  true/false (default: false)
%   'scale'                 scale factor applied to vertices (default: 1)
%   'combineObjects'        combine multiple objects into one FV (default: true)
%
% OUTPUT:
% fv: structure with fields 'faces' and 'vertices' (and 'facevertexcdata' if available)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

opts = struct('removeUnusedVertices', false, 'scale', 1, 'combineObjects', true);
opts = parseOptions(opts, varargin{:});

if ~exist('fileName', 'var') || isempty(fileName)
    [file, path] = uigetfile('*.3mf', 'Select *.3mf file');
    if isequal(file, 0)
        fv = struct('faces', [], 'vertices', []);
        return;
    end
    fileName = [path file];
end

if ~exist(fileName, 'file')
    error('threeMfToPatch:fileNotFound', 'File not found: %s', fileName);
end

tempDir = tempname;
mkdir(tempDir);
cleanupDir = onCleanup(@() cleanupTempDir(tempDir));
unzip(fileName, tempDir);

modelPath = fullfile(tempDir, '3D', '3dmodel.model');
if ~exist(modelPath, 'file')
    modelFiles = dir(fullfile(tempDir, '3D', '*.model'));
    if ~isempty(modelFiles)
        modelPath = fullfile(modelFiles(1).folder, modelFiles(1).name);
    else
        rmdir(tempDir, 's');
        error('threeMfToPatch:missingModel', 'No 3D model file found in 3MF.');
    end
end

doc = xmlread(modelPath);

materialGroups = parseMaterialGroupsXml(doc);
objects = parseObjectsXml(doc, materialGroups);
buildIds = parseBuildIdsXml(doc);

if isempty(objects)
    fv = struct('faces', [], 'vertices', []);
    return;
end

[objectsOut] = resolveBuildObjects(objects, buildIds);

if opts.combineObjects && numel(objectsOut) > 1
    objectsOut = combineObjects(objectsOut);
end

if numel(objectsOut) == 1
    fv = objectsOut(1);
else
    fv = objectsOut;
end

if opts.scale ~= 1
    for i = 1:numel(fv)
        fv(i).vertices = fv(i).vertices .* opts.scale;
    end
end

if opts.removeUnusedVertices
    for i = 1:numel(fv)
        [fv(i).vertices, fv(i).faces] = mergeDuplicateVertices(fv(i).vertices, fv(i).faces);
    end
end
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('threeMfToPatch:invalidOptions', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k + 1};
        switch key
            case "removeunusedvertices"
                opts.removeUnusedVertices = logical(val);
            case "scale"
                opts.scale = double(val);
            case "combineobjects"
                opts.combineObjects = logical(val);
            otherwise
                error('threeMfToPatch:invalidOption', 'Unknown option "%s".', key);
        end
    end
end

function materialGroups = parseMaterialGroupsXml(doc)
    materialGroups = containers.Map('KeyType', 'double', 'ValueType', 'any');
    nodes = doc.getElementsByTagName('*');
    for i = 0:(nodes.getLength - 1)
        node = nodes.item(i);
        if isElementNamed(node, 'colorgroup')
            idVal = getAttrDouble(node, 'id');
            if ~isfinite(idVal)
                continue;
            end
            colors = parseColorEntries(node, 'color');
            materialGroups(idVal) = colors;
        elseif isElementNamed(node, 'basematerials')
            idVal = getAttrDouble(node, 'id');
            if ~isfinite(idVal)
                continue;
            end
            colors = parseColorEntries(node, 'base');
            materialGroups(idVal) = colors;
        else
            continue;
        end
    end
end

function objects = parseObjectsXml(doc, materialGroups)
    objects = struct('id', {}, 'faces', {}, 'vertices', {}, 'facevertexcdata', {}, 'components', {});
    nodes = doc.getElementsByTagName('*');
    for i = 0:(nodes.getLength - 1)
        node = nodes.item(i);
        if ~isElementNamed(node, 'object')
            continue;
        end
        objId = getAttrDouble(node, 'id');
        if ~isfinite(objId)
            continue;
        end
        meshNode = findChildElement(node, 'mesh');
        compNode = findChildElement(node, 'components');

        obj = struct('id', objId, 'faces', [], 'vertices', [], 'facevertexcdata', [], 'components', []);
        if ~isempty(meshNode)
            defaultPid = getAttrDouble(node, 'pid');
            defaultP = getDefaultP(node);
            [vertices, faces, colors] = parseMeshXml(meshNode, materialGroups, defaultPid, defaultP);
            obj.vertices = vertices;
            obj.faces = faces;
            if ~isempty(colors)
                obj.facevertexcdata = colors;
            end
        elseif ~isempty(compNode)
            obj.components = parseComponentsXml(compNode);
        else
            continue;
        end
        objects(end+1) = obj; %#ok<AGROW>
    end
end

function comps = parseComponentsXml(compNode)
    comps = struct('objectId', {}, 'transform', {});
    nodes = compNode.getChildNodes;
    for i = 0:(nodes.getLength - 1)
        node = nodes.item(i);
        if ~isElementNamed(node, 'component')
            continue;
        end
        objId = getAttrDouble(node, 'objectid');
        if ~isfinite(objId)
            continue;
        end
        t = getAttrString(node, 'transform');
        if strlength(t) > 0
            vals = str2double(strsplit(t));
        else
            vals = [];
        end
        comps(end+1) = struct('objectId', objId, 'transform', vals); %#ok<AGROW>
    end
end

function [vertices, faces, colors] = parseMeshXml(meshNode, materialGroups, defaultPid, defaultP)
    vertices = [];
    faces = [];
    colors = [];

    vertContainer = findChildElement(meshNode, 'vertices');
    triContainer = findChildElement(meshNode, 'triangles');
    if isempty(vertContainer) || isempty(triContainer)
        return;
    end

    vertNodes = vertContainer.getChildNodes;
    verts = zeros(vertNodes.getLength, 3);
    vCount = 0;
    for i = 0:(vertNodes.getLength - 1)
        vNode = vertNodes.item(i);
        if ~isElementNamed(vNode, 'vertex')
            continue;
        end
        x = getAttrDouble(vNode, 'x');
        y = getAttrDouble(vNode, 'y');
        z = getAttrDouble(vNode, 'z');
        if any(~isfinite([x y z]))
            continue;
        end
        vCount = vCount + 1;
        verts(vCount, :) = [x y z];
    end
    vertices = verts(1:vCount, :);

    triNodes = triContainer.getChildNodes;
    faces = zeros(triNodes.getLength, 3);
    triPid = nan(triNodes.getLength, 1);
    triP = nan(triNodes.getLength, 3);
    hasColor = false;
    tCount = 0;

    for i = 0:(triNodes.getLength - 1)
        tNode = triNodes.item(i);
        if ~isElementNamed(tNode, 'triangle')
            continue;
        end
        v1 = getAttrDouble(tNode, 'v1');
        v2 = getAttrDouble(tNode, 'v2');
        v3 = getAttrDouble(tNode, 'v3');
        if any(~isfinite([v1 v2 v3]))
            continue;
        end
        tCount = tCount + 1;
        faces(tCount, :) = [v1 v2 v3] + 1;

        pid = getAttrDouble(tNode, 'pid');
        p1 = getAttrDouble(tNode, 'p1');
        p2 = getAttrDouble(tNode, 'p2');
        p3 = getAttrDouble(tNode, 'p3');
        p = getAttrDouble(tNode, 'p');
        if ~isfinite(p)
            p = getAttrDouble(tNode, 'pindex');
        end

        if ~isfinite(pid) && isfinite(defaultPid)
            pid = defaultPid;
        end
        if ~isfinite(p) && isfinite(defaultP)
            p = defaultP;
        end

        if isfinite(pid) && isfinite(p1) && isfinite(p2) && isfinite(p3)
            hasColor = true;
            triPid(tCount) = pid;
            triP(tCount, :) = [p1 p2 p3];
        elseif isfinite(pid) && isfinite(p)
            hasColor = true;
            triPid(tCount) = pid;
            triP(tCount, :) = [p p p];
        end
    end

    faces = faces(1:tCount, :);
    triPid = triPid(1:tCount);
    triP = triP(1:tCount, :);

    if hasColor
        [vertices, faces, colors] = expandVertexColors(vertices, faces, triPid, triP, materialGroups);
    end
end

function val = getAttrDouble(node, name)
    val = NaN;
    if isempty(node) || ~node.hasAttributes
        return;
    end
    attr = node.getAttributes.getNamedItem(char(name));
    if isempty(attr)
        return;
    end
    val = str2double(char(attr.getValue));
end

function val = getAttrString(node, name)
    val = "";
    if isempty(node) || ~node.hasAttributes
        return;
    end
    attr = node.getAttributes.getNamedItem(char(name));
    if isempty(attr)
        return;
    end
    val = string(char(attr.getValue));
end

function tf = isElementNamed(node, name)
    tf = false;
    if isempty(node) || node.getNodeType ~= node.ELEMENT_NODE
        return;
    end
    nodeName = string(char(node.getNodeName));
    localName = "";
    try
        localName = string(char(node.getLocalName));
    catch
    end
    if strlength(localName) > 0
        tf = strcmpi(localName, name);
    else
        if contains(nodeName, ":")
            nodeName = extractAfter(nodeName, ":");
        end
        tf = strcmpi(nodeName, name);
    end
end

function child = findChildElement(parent, name)
    child = [];
    if isempty(parent)
        return;
    end
    nodes = parent.getChildNodes;
    for i = 0:(nodes.getLength - 1)
        node = nodes.item(i);
        if isElementNamed(node, name)
            child = node;
            return;
        end
    end
end

function [vOut, fOut, cOut] = expandVertexColors(vertices, faces, triPid, triP, materialGroups)
    numFaces = size(faces, 1);
    vOut = zeros(numFaces * 3, 3);
    cOut = zeros(numFaces * 3, 3);
    fOut = reshape(1:(numFaces * 3), 3, [])';

    for i = 1:numFaces
        idx = faces(i, :);
        base = (i - 1) * 3;
        vOut(base + 1, :) = vertices(idx(1), :);
        vOut(base + 2, :) = vertices(idx(2), :);
        vOut(base + 3, :) = vertices(idx(3), :);

        pid = triPid(i);
        if isKey(materialGroups, pid)
            colors = materialGroups(pid);
            p1 = triP(i, 1) + 1;
            p2 = triP(i, 2) + 1;
            p3 = triP(i, 3) + 1;
            if p1 <= size(colors, 1) && p2 <= size(colors, 1) && p3 <= size(colors, 1)
                cOut(base + 1, :) = colors(p1, :);
                cOut(base + 2, :) = colors(p2, :);
                cOut(base + 3, :) = colors(p3, :);
            end
        end
    end
end

function objectsOut = combineObjects(objects)
    vertices = [];
    faces = [];
    colors = [];
    offset = 0;
    for i = 1:numel(objects)
        v = objects(i).vertices;
        f = objects(i).faces;
        faces = [faces; f + offset]; %#ok<AGROW>
        vertices = [vertices; v]; %#ok<AGROW>
        if isfield(objects(i), 'facevertexcdata') && ~isempty(objects(i).facevertexcdata)
            colors = [colors; objects(i).facevertexcdata]; %#ok<AGROW>
        elseif ~isempty(colors)
            colors = [colors; zeros(size(v, 1), 3)]; %#ok<AGROW>
        end
        offset = offset + size(v, 1);
    end
    objectsOut = struct('faces', faces, 'vertices', vertices);
    if ~isempty(colors)
        objectsOut.facevertexcdata = colors;
    end
end

function [verticesOut, facesOut] = mergeDuplicateVertices(vertices, faces)
    [verticesOut, ~, idx] = unique(vertices, 'rows');
    facesOut = idx(faces);
end

function colors = parseColorEntries(parentNode, elementName)
    colors = [];
    childNodes = parentNode.getChildNodes;
    for j = 0:(childNodes.getLength - 1)
        child = childNodes.item(j);
        if ~isElementNamed(child, elementName)
            continue;
        end
        hex = getAttrString(child, 'color');
        if strlength(hex) == 0
            hex = getAttrString(child, 'displaycolor');
        end
        if strlength(hex) == 7 && startsWith(hex, "#")
            hex = extractAfter(hex, 1);
        end
        if strlength(hex) ~= 6
            continue;
        end
        colors(end+1, :) = [ ...
            hex2dec(char(extractBefore(extractAfter(hex, 0), 3))), ...
            hex2dec(char(extractBefore(extractAfter(hex, 2), 3))), ...
            hex2dec(char(extractBefore(extractAfter(hex, 4), 3))) ...
            ] / 255; %#ok<AGROW>
    end
end

function defaultP = getDefaultP(node)
    defaultP = getAttrDouble(node, 'p');
    if ~isfinite(defaultP)
        defaultP = getAttrDouble(node, 'pindex');
    end
end

function buildIds = parseBuildIdsXml(doc)
    buildIds = [];
    nodes = doc.getElementsByTagName('*');
    for i = 0:(nodes.getLength - 1)
        node = nodes.item(i);
        if ~isElementNamed(node, 'build')
            continue;
        end
        childNodes = node.getChildNodes;
        for j = 0:(childNodes.getLength - 1)
            child = childNodes.item(j);
            if ~isElementNamed(child, 'item')
                continue;
            end
            objId = getAttrDouble(child, 'objectid');
            if isfinite(objId)
                buildIds(end+1) = objId; %#ok<AGROW>
            end
        end
    end
end

function objectsOut = resolveBuildObjects(objects, buildIds)
    if isempty(objects)
        objectsOut = objects;
        return;
    end
    objMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1:numel(objects)
        objMap(objects(i).id) = objects(i);
    end
    if isempty(buildIds)
        buildIds = [objects.id];
    end
    objectsOut = struct('faces', {}, 'vertices', {}, 'facevertexcdata', {});
    for i = 1:numel(buildIds)
        objId = buildIds(i);
        if ~isKey(objMap, objId)
            continue;
        end
        mesh = resolveObject(objMap, objId, []);
        if ~isempty(mesh.vertices) && ~isempty(mesh.faces)
            objectsOut(end+1) = mesh; %#ok<AGROW>
        end
    end
end

function mesh = resolveObject(objMap, objId, stack)
    if any(stack == objId)
        mesh = struct('faces', [], 'vertices', [], 'facevertexcdata', []);
        return;
    end
    obj = objMap(objId);
    if ~isempty(obj.vertices) && ~isempty(obj.faces)
        mesh = struct('faces', obj.faces, 'vertices', obj.vertices, 'facevertexcdata', []);
        if isfield(obj, 'facevertexcdata')
            mesh.facevertexcdata = obj.facevertexcdata;
        end
        return;
    end
    stack = [stack, objId];
    if isempty(obj.components)
        mesh = struct('faces', [], 'vertices', [], 'facevertexcdata', []);
        return;
    end
    meshes = struct('faces', {}, 'vertices', {}, 'facevertexcdata', {});
    for i = 1:numel(obj.components)
        comp = obj.components(i);
        if ~isKey(objMap, comp.objectId)
            continue;
        end
        childMesh = resolveObject(objMap, comp.objectId, stack);
        if isempty(childMesh.vertices) || isempty(childMesh.faces)
            continue;
        end
        if ~isempty(comp.transform)
            childMesh.vertices = applyTransform(childMesh.vertices, comp.transform);
        end
        meshes(end+1) = childMesh; %#ok<AGROW>
    end
    if isempty(meshes)
        mesh = struct('faces', [], 'vertices', [], 'facevertexcdata', []);
        return;
    end
    mesh = combineObjects(meshes);
end

function vOut = applyTransform(vertices, vals)
    if numel(vals) ~= 12
        vOut = vertices;
        return;
    end
    T = reshape(vals, 4, 3)';
    R = T(:, 1:3);
    t = T(:, 4)';
    vOut = (vertices * R') + t;
end

function cleanupTempDir(tempDir)
    if exist(tempDir, 'dir')
        try
            rmdir(tempDir, 's');
        catch
        end
    end
end
