function [hull, sliceInfo] = scanHull(X, alpha, numSeg, nGon, capLoops, varargin)
% scanHull computes a detector-space hull for APT data and returns FV patches.
%
% [hull, sliceInfo] = scanHull(X, alpha, numSeg, nGon, capLoops)
% [hull, sliceInfo] = scanHull(X, ..., 'showPlot', true)
%
% INPUT
% X:            Nx2 detector coordinates or pos table with detx/dety
% alpha:        alpha value for alpha shape (default: 3)
% numSeg:       number of slices along ion index (default: 10)
% nGon:         number of angular segments per slice (default: 64)
% capLoops:     number of concentric loops for caps (default: 4)
%
% name-value options:
%   'showPlot'      true/false (default: false)
%   'useVoronoi'    true/false (default: false)
%   'noiseLimit'    Voronoi filter noise limit (default: 0.10)
%   'rMax'          max radius for ray intersection (default: max radius)
%   'debug'         true/false (default: false)
%
% OUTPUT
% hull:       struct array with fields vertices and faces (mantle, top, bottom)
% sliceInfo:  [idxAvg, area] per slice
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

if istable(X)
    if all(ismember({'detx','dety'}, X.Properties.VariableNames))
        X = [X.detx, X.dety];
    else
        error('scanHull:invalidInput', 'Table input must contain detx and dety.');
    end
end

if size(X, 2) > 2
    X = X(:, 1:2);
end

if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 3;
end
if ~exist('numSeg', 'var') || isempty(numSeg)
    numSeg = 10;
end
if ~exist('nGon', 'var') || isempty(nGon)
    nGon = 64;
end
if ~exist('capLoops', 'var') || isempty(capLoops)
    capLoops = 4;
end

opts = struct('showPlot', false, 'useVoronoi', false, 'noiseLimit', 0.10, ...
    'rMax', [], 'debug', false);
opts = parseOptions(opts, varargin{:});

addGeom2dToPath();

X = double(X);
numAtoms = size(X, 1);
atPerSlice = numAtoms / numSeg;

if isempty(opts.rMax)
    r = hypot(X(:,1), X(:,2));
    opts.rMax = max(r);
end

sliceInfo = zeros(numSeg, 2);
hullSlices = cell(numSeg, 1);

for slice = 1:numSeg
    if opts.debug
        fprintf('segment %d of %d\n', slice, numSeg);
    end

    idxBeg = max(1, floor(atPerSlice * (slice - 1)) + 1);
    idxEnd = max(idxBeg, floor(atPerSlice * slice) - 1);
    idxAvg = round((idxBeg + idxEnd) / 2);

    sliceX = X(idxBeg:idxEnd, 1:2);

    if opts.useVoronoi
        sliceX = VoronoiFilter2D(sliceX, opts.noiseLimit);
    end

    [area, alphaShape] = alphavol(sliceX, alpha);
    sliceInfo(slice, :) = [idxAvg, area];

    loop = [sliceX(alphaShape.bnd(:,1),1), sliceX(alphaShape.bnd(:,2),2)];
    loop = [loop; loop(1,:)];

    inc = 2*pi()/nGon;
    outline = zeros(nGon, 2);
    for n = 1:nGon
        [x, y] = pol2cart(inc * n, opts.rMax);
        edge = [0 0 x y];
        try
            tmp = intersectEdgePolygon(edge, loop);
            outline(n,:) = tmp(1,:);
        catch
            if opts.debug
                warning('scanHull:segmentFailed', 'error in slice %d, segment %d', slice, n);
            end
        end
    end

    hullSlices{slice} = [outline, double(ones(size(outline,1),1)) * double(idxAvg)];
end

upperEdge = hullSlices{1};
upperEdge(:,3) = ones(size(upperEdge,1),1);

lowerEdge = hullSlices{end};
lowerEdge(:,3) = ones(size(lowerEdge,1),1) * numAtoms;

hullSlices = [{upperEdge} hullSlices {lowerEdge}];

mantle.vertices = [];
mantle.faces = [];
for s = 1:length(hullSlices)
    mantle.vertices = [mantle.vertices; hullSlices{s}]; %#ok<AGROW>
end

for s = 0:length(hullSlices)-2
    for n = 1:nGon
        mantle.faces = [mantle.faces; ...
            n + s*nGon, n + 1 + s*nGon, n + (s+1)*nGon; ...
            n + s*nGon, n + (s+1)*nGon, n - 1 + (s+1)*nGon]; %#ok<AGROW>
    end
end

rim = hullSlices{1}(:,1:2);
topCap = buildCap(rim, capLoops, nGon);
topCap.vertices(:,3) = ones(size(topCap.vertices,1),1);

rim = hullSlices{end}(:,1:2);
bottomCap = buildCap(rim, capLoops, nGon);
bottomCap.vertices(:,3) = ones(size(bottomCap.vertices,1),1) * numAtoms;

if opts.showPlot
    figure('Name','detector hull');
    p = patch(mantle);
    set(p,'FaceColor',[0 1 1]);
    hold on;
    tp = patch(topCap);
    set(tp,'FaceColor',[1 0 0]);
    bp = patch(bottomCap);
    set(bp,'FaceColor',[0 0 1]);
    axis equal;
end

hull = [mantle, topCap, bottomCap];
end

function opts = parseOptions(opts, varargin)
    if mod(numel(varargin), 2) ~= 0
        error('scanHull:invalidOptions', 'Options must be name-value pairs.');
    end
    for k = 1:2:numel(varargin)
        key = lower(string(varargin{k}));
        val = varargin{k + 1};
        switch key
            case "showplot"
                opts.showPlot = logical(val);
            case "usevoronoi"
                opts.useVoronoi = logical(val);
            case "noiselimit"
                opts.noiseLimit = double(val);
            case "rmax"
                opts.rMax = double(val);
            case "debug"
                opts.debug = logical(val);
            otherwise
                error('scanHull:invalidOption', 'Unknown option "%s".', key);
        end
    end
end

function addGeom2dToPath()
    if exist('intersectEdgePolygon', 'file') == 2 && exist('alphavol', 'file') == 2
        return;
    end
    rootDir = fileparts(mfilename('fullpath'));
    toolboxRoot = fileparts(rootDir);
    addpath(fullfile(toolboxRoot, 'utilities_geom2d'));
    addpath(fullfile(toolboxRoot, 'utilities_geom2d', 'polygons2d'));
    addpath(fullfile(toolboxRoot, 'utilities_geom2d', 'geom2d'));
    addpath(fullfile(toolboxRoot, 'utilities'));
end

function cap = buildCap(rim, capLoops, nGon)
    cap.vertices = [];
    for pt = 1:size(rim, 1)
        [theta, rho] = cart2pol(rim(pt,1), rim(pt,2));
        for lp = 1:capLoops
            [x, y] = pol2cart(theta, rho / capLoops * lp);
            cap.vertices(end+1,:) = [x y]; %#ok<AGROW>
        end
    end

    cap.vertices(end+1,:) = [0, 0];
    numVerts = size(cap.vertices, 1);
    tris = [];
    for pt = 1:(nGon - 1)
        tris(end+1,:) = [numVerts, 1 + (pt - 1) * capLoops, 1 + pt * capLoops]; %#ok<AGROW>
        for lp = 1:(capLoops - 1)
            tris(end+1,:) = [1+(pt-1)*capLoops+(lp-1), 1+(pt-1)*capLoops+lp, 1+pt*capLoops+(lp-1)]; %#ok<AGROW>
            tris(end+1,:) = [1+(pt-1)*capLoops+lp, 1+pt*capLoops+lp, 1+pt*capLoops+(lp-1)]; %#ok<AGROW>
        end
    end
    tris(end+1,:) = [numVerts, 1 + (nGon - 1) * capLoops, 1 + 0 * capLoops]; %#ok<AGROW>
    for lp = 1:(capLoops - 1)
        tris(end+1,:) = [1+(nGon-1)*capLoops+(lp-1), 1+(nGon-1)*capLoops+lp, 1+0*capLoops+(lp-1)]; %#ok<AGROW>
        tris(end+1,:) = [1+(nGon-1)*capLoops+lp, 1+0*capLoops+lp, 1+0*capLoops+(lp-1)]; %#ok<AGROW>
    end
    cap.faces = tris;
end
