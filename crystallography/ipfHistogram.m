function [histData, ax] = ipfHistogram(fv, orientationInput, varargin)
% IPFHISTOGRAM Create orientation histogram in the IPF standard triangle
%
%   [HISTDATA, AX] = IPFHISTOGRAM(FV, CRYSTALORIENTATION) creates a histogram
%   of surface normal orientations in the standard triangle (IPF) for cubic
%   crystals. The histogram shows the frequency distribution of crystallographic
%   orientations of the mesh surface.
%
%   [HISTDATA, AX] = IPFHISTOGRAM(FV, CRYSTALORIENTATION, 'property', PROP, ...)
%   creates a property-weighted histogram where PROP is a per-vertex array of
%   values (e.g., concentration) to be aggregated by orientation.
%
%   Input:
%       fv - Mesh structure with fields:
%            fv.vertices - Nx3 array of vertex coordinates
%            fv.faces    - Mx3 array of face indices (triangles)
%
%       orientationInput - Either:
%                          - 3x3 rotation matrix (crystal to world)
%                          - crystalOrientation object
%
%   Optional Name-Value Pairs:
%       'property'     - Nx1 array of property values per vertex (default: [])
%       'binning'      - 'angular' (default) or 'adaptive'
%       'resolution'   - Angular resolution in degrees for angular binning (default: 2)
%       'nBins'        - Number of bins for adaptive binning (default: 500)
%       'aggregation'  - 'count', 'mean', 'sum' (default: 'count' if no property,
%                        'mean' if property provided)
%       'areaWeighted' - Logical, weight by face area (default: false)
%       'colormap'     - Colormap name or Nx3 array (default: 'parula')
%       'plot'         - Logical, create plot (default: true)
%       'showColorbar' - Logical, show colorbar (default: true)
%
%   Output:
%       histData - Structure with fields:
%                  .bins          - Bin center positions [h, k, l]
%                  .counts        - Count/weight per bin
%                  .values        - Aggregated property values per bin
%                  .crystalNormals - Nx3 crystal normals for all vertices
%                  .binIndices    - Bin assignment for each vertex
%
%       ax - Axes handle (if plot is true)
%
%   Example:
%       % Simple count histogram
%       [X,Y,Z] = sphere(50);
%       fv = surf2patch(X,Y,Z,'triangles');
%       ori = crystalOrientation(eye(3), 'cubic');
%       [hist, ax] = ipfHistogram(fv, ori);
%
%       % Property-weighted histogram
%       concentration = rand(size(fv.vertices, 1), 1);
%       [hist, ax] = ipfHistogram(fv, ori, 'property', concentration, ...
%                                 'resolution', 1, 'aggregation', 'mean');
%
%       % Area-weighted adaptive binning
%       [hist, ax] = ipfHistogram(fv, ori, 'binning', 'adaptive', ...
%                                 'nBins', 1000, 'areaWeighted', true);
%
%   See also: crystalOrientation, ipfColor, ipfMesh

    % Determine input type and extract rotation matrix
    if isa(orientationInput, 'crystalOrientation')
        R = orientationInput.getRotationMatrix();
    elseif isnumeric(orientationInput) && all(size(orientationInput) == [3, 3])
        R = orientationInput;
    else
        error('orientationInput must be a 3x3 rotation matrix or a crystalOrientation object');
    end

    % Parse inputs
    p = inputParser;
    addRequired(p, 'fv', @(x) isstruct(x) && isfield(x,'vertices') && isfield(x,'faces'));
    addRequired(p, 'orientationInput');
    addParameter(p, 'property', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
    addParameter(p, 'binning', 'angular', @(x) ismember(lower(x), {'angular', 'adaptive'}));
    addParameter(p, 'resolution', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'nBins', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'aggregation', '', @(x) isempty(x) || ismember(lower(x), {'count', 'mean', 'sum'}));
    addParameter(p, 'areaWeighted', false, @islogical);
    addParameter(p, 'colormap', 'parula', @(x) ischar(x) || (isnumeric(x) && size(x,2)==3));
    addParameter(p, 'plot', true, @islogical);
    addParameter(p, 'showColorbar', true, @islogical);
    parse(p, fv, orientationInput, varargin{:});

    property = p.Results.property(:);
    binningMode = lower(p.Results.binning);
    resolution = p.Results.resolution;
    nBins = round(p.Results.nBins);
    aggregation = lower(p.Results.aggregation);
    areaWeighted = p.Results.areaWeighted;
    cmapName = p.Results.colormap;
    doPlot = p.Results.plot;
    showColorbar = p.Results.showColorbar;

    % Set default aggregation
    if isempty(aggregation)
        if isempty(property)
            aggregation = 'count';
        else
            aggregation = 'mean';
        end
    end

    % Validate property size
    vertices = fv.vertices;
    faces = fv.faces;
    nVertices = size(vertices, 1);

    if ~isempty(property) && length(property) ~= nVertices
        error('Property array must have one value per vertex (%d vertices)', nVertices);
    end

    %% Step 1: Calculate vertex normals
    vertexNormals = calculateVertexNormals(vertices, faces);

    %% Step 2: Transform to crystal coordinates
    R_world2crystal = R';
    crystalNormals = (R_world2crystal * vertexNormals')';

    %% Step 3: Reduce to fundamental zone (standard triangle)
    reducedNormals = reduceToFundamentalZone(crystalNormals);

    %% Step 4: Calculate vertex weights if area-weighted
    if areaWeighted
        vertexWeights = calculateVertexWeights(vertices, faces);
    else
        vertexWeights = ones(nVertices, 1);
    end

    %% Step 5: Create bins and assign vertices to bins
    if strcmp(binningMode, 'angular')
        [binCenters, binEdges] = createAngularBins(resolution);
    else
        [binCenters, binEdges] = createAdaptiveBins(reducedNormals, nBins);
    end

    binIndices = assignToBins(reducedNormals, binCenters);

    %% Step 6: Aggregate values in each bin
    nBinsActual = size(binCenters, 1);
    counts = zeros(nBinsActual, 1);
    values = zeros(nBinsActual, 1);

    for b = 1:nBinsActual
        inBin = (binIndices == b);
        weights = vertexWeights(inBin);
        counts(b) = sum(weights);

        if ~isempty(property)
            propVals = property(inBin);
            switch aggregation
                case 'mean'
                    if sum(weights) > 0
                        values(b) = sum(propVals .* weights) / sum(weights);
                    else
                        values(b) = NaN;
                    end
                case 'sum'
                    values(b) = sum(propVals .* weights);
                case 'count'
                    values(b) = counts(b);
            end
        else
            values(b) = counts(b);
        end
    end

    %% Build output structure
    histData.bins = binCenters;
    histData.binEdges = binEdges;
    histData.counts = counts;
    histData.values = values;
    histData.crystalNormals = crystalNormals;
    histData.reducedNormals = reducedNormals;
    histData.binIndices = binIndices;
    histData.aggregation = aggregation;
    histData.areaWeighted = areaWeighted;

    %% Step 7: Visualization
    if doPlot
        ax = plotHistogramTriangle(histData, cmapName, showColorbar, aggregation, ~isempty(property));
    else
        ax = [];
    end
end

%% Helper Functions

function normals = calculateVertexNormals(vertices, faces)
% CALCULATEVERTEXNORMALS Compute vertex normals by averaging face normals

    nVertices = size(vertices, 1);
    nFaces = size(faces, 1);

    normals = zeros(nVertices, 3);

    for i = 1:nFaces
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);

        p1 = vertices(v1, :);
        p2 = vertices(v2, :);
        p3 = vertices(v3, :);

        edge1 = p2 - p1;
        edge2 = p3 - p1;
        faceNormal = cross(edge1, edge2);

        normals(v1, :) = normals(v1, :) + faceNormal;
        normals(v2, :) = normals(v2, :) + faceNormal;
        normals(v3, :) = normals(v3, :) + faceNormal;
    end

    % Normalize
    for i = 1:nVertices
        n = norm(normals(i, :));
        if n > 0
            normals(i, :) = normals(i, :) / n;
        end
    end
end

function reduced = reduceToFundamentalZone(directions)
% REDUCETOFUNDAMENTALZONE Reduce directions to standard triangle using cubic symmetry
%   Returns [h, k, l] with 0 <= h <= k <= l

    nDirs = size(directions, 1);
    reduced = zeros(nDirs, 3);

    for i = 1:nDirs
        dir = directions(i, :);
        n = norm(dir);
        if n > 0
            dir = dir / n;
        end

        % Take absolute values and sort ascending
        dirAbs = abs(dir);
        dirSorted = sort(dirAbs, 'ascend');
        reduced(i, :) = dirSorted;
    end
end

function weights = calculateVertexWeights(vertices, faces)
% CALCULATEVERTEXWEIGHTS Calculate area-based weights for each vertex

    nVertices = size(vertices, 1);
    nFaces = size(faces, 1);
    weights = zeros(nVertices, 1);

    for i = 1:nFaces
        v1 = faces(i, 1);
        v2 = faces(i, 2);
        v3 = faces(i, 3);

        p1 = vertices(v1, :);
        p2 = vertices(v2, :);
        p3 = vertices(v3, :);

        % Face area = 0.5 * |cross product|
        edge1 = p2 - p1;
        edge2 = p3 - p1;
        faceArea = 0.5 * norm(cross(edge1, edge2));

        % Distribute 1/3 of face area to each vertex
        weights(v1) = weights(v1) + faceArea / 3;
        weights(v2) = weights(v2) + faceArea / 3;
        weights(v3) = weights(v3) + faceArea / 3;
    end
end

function [binCenters, binEdges] = createAngularBins(resolution)
% CREATEANGULARBINS Create uniform angular bins in the standard triangle
%   Resolution is in degrees

    % Convert resolution to step in (h,k) parameter space
    % For the standard triangle: 0 <= h <= k <= 1 (with l=1 normalized)
    % Angular resolution corresponds to spacing in h,k space

    % Approximate: 1 degree ~ 0.0175 radians, and the triangle spans ~35 degrees
    % So step size ~ resolution / 35
    step = resolution / 45;  % Empirical scaling
    step = max(0.01, min(0.2, step));  % Clamp to reasonable range

    % Generate grid points in the triangular domain
    hVals = 0:step:1;
    kVals = 0:step:1;

    binCenters = [];
    for h = hVals
        for k = kVals
            if h <= k  % Only include points in the triangle
                binCenters = [binCenters; h, k, 1];
            end
        end
    end

    % Normalize bin centers
    for i = 1:size(binCenters, 1)
        binCenters(i, :) = binCenters(i, :) / norm(binCenters(i, :));
    end

    binEdges = step;  % Store step size as edge info
end

function [binCenters, binEdges] = createAdaptiveBins(reducedNormals, nBins)
% CREATEADAPTIVEBINS Create adaptive bins based on data density
%   Uses k-means clustering for adaptive bin placement

    % Normalize all directions
    norms = sqrt(sum(reducedNormals.^2, 2));
    norms(norms == 0) = 1;
    normalizedDirs = reducedNormals ./ norms;

    % Limit nBins to number of unique directions
    uniqueDirs = unique(round(normalizedDirs * 1000) / 1000, 'rows');
    nBins = min(nBins, size(uniqueDirs, 1));

    if nBins < 2
        binCenters = mean(normalizedDirs, 1);
        binEdges = [];
        return;
    end

    % Use k-means clustering
    try
        [~, binCenters] = kmeans(normalizedDirs, nBins, ...
            'MaxIter', 100, 'Replicates', 3, 'Display', 'off');
    catch
        % Fallback to uniform if kmeans fails
        [binCenters, binEdges] = createAngularBins(2);
        return;
    end

    % Ensure bin centers are in fundamental zone
    binCenters = reduceToFundamentalZone(binCenters);

    binEdges = [];  % Adaptive bins don't have regular edges
end

function binIndices = assignToBins(reducedNormals, binCenters)
% ASSIGNTOBINS Assign each direction to nearest bin center

    nDirs = size(reducedNormals, 1);
    nBins = size(binCenters, 1);
    binIndices = zeros(nDirs, 1);

    % Normalize directions
    for i = 1:nDirs
        n = norm(reducedNormals(i, :));
        if n > 0
            reducedNormals(i, :) = reducedNormals(i, :) / n;
        end
    end

    % Normalize bin centers
    for i = 1:nBins
        n = norm(binCenters(i, :));
        if n > 0
            binCenters(i, :) = binCenters(i, :) / n;
        end
    end

    % Find nearest bin for each direction using dot product (cosine similarity)
    for i = 1:nDirs
        dir = reducedNormals(i, :);
        similarities = binCenters * dir';
        [~, binIndices(i)] = max(similarities);
    end
end

function ax = plotHistogramTriangle(histData, cmapName, showColorbar, aggregation, hasProperty)
% PLOTHISTOGRAMTRIANGLE Render the histogram as filled contour on standard triangle

    figure('Color', 'w');
    ax = axes;
    hold(ax, 'on');

    bins = histData.bins;
    values = histData.values;

    % Remove bins with no data
    validBins = ~isnan(values) & histData.counts > 0;
    bins = bins(validBins, :);
    values = values(validBins);

    if isempty(bins)
        warning('No valid data to plot');
        return;
    end

    % Get 2D coordinates for bin centers
    [xBins, yBins] = direction2stereo(bins);

    % Create a fine mesh for interpolation
    nGrid = 100;
    [meshVertices, meshFaces, meshColors] = createInterpolatedMesh(xBins, yBins, values, nGrid);

    % Plot the mesh
    if ~isempty(meshVertices)
        patch(ax, 'Faces', meshFaces, 'Vertices', meshVertices, ...
              'FaceVertexCData', meshColors, ...
              'FaceColor', 'interp', ...
              'EdgeColor', 'none');
    end

    % Draw triangle boundary
    nBoundary = 100;
    t = linspace(0, 1, nBoundary);

    % Edge [001] to [011]: h=0, k varies 0 to 1
    bx1 = zeros(1, nBoundary); by1 = zeros(1, nBoundary);
    for i = 1:nBoundary
        [bx1(i), by1(i)] = direction2stereoSingle(0, t(i), 1);
    end
    plot(ax, bx1, by1, 'k-', 'LineWidth', 1.5);

    % Edge [011] to [111]: k=1, h varies 0 to 1
    bx2 = zeros(1, nBoundary); by2 = zeros(1, nBoundary);
    for i = 1:nBoundary
        [bx2(i), by2(i)] = direction2stereoSingle(t(i), 1, 1);
    end
    plot(ax, bx2, by2, 'k-', 'LineWidth', 1.5);

    % Edge [111] to [001]: h=k, both vary 1 to 0
    bx3 = zeros(1, nBoundary); by3 = zeros(1, nBoundary);
    for i = 1:nBoundary
        val = 1 - t(i);
        [bx3(i), by3(i)] = direction2stereoSingle(val, val, 1);
    end
    plot(ax, bx3, by3, 'k-', 'LineWidth', 1.5);

    % Add vertex labels
    [x001, y001] = direction2stereoSingle(0, 0, 1);
    [x011, y011] = direction2stereoSingle(0, 1, 1);
    [x111, y111] = direction2stereoSingle(1, 1, 1);

    text(ax, x001 - 0.03, y001 - 0.02, '[001]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
    text(ax, x011 + 0.03, y011 - 0.02, '[011]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    text(ax, x111, y111 + 0.03, '[111]', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Set colormap
    if ischar(cmapName)
        colormap(ax, cmapName);
    else
        colormap(ax, cmapName);
    end

    % Add colorbar
    if showColorbar
        cb = colorbar(ax);
        if hasProperty
            switch aggregation
                case 'mean'
                    cb.Label.String = 'Mean Property Value';
                case 'sum'
                    cb.Label.String = 'Sum of Property Values';
                case 'count'
                    cb.Label.String = 'Count';
            end
        else
            cb.Label.String = 'Count';
        end
    end

    % Format axes
    axis(ax, 'equal');
    axis(ax, 'off');
    title(ax, 'IPF Histogram', 'FontSize', 14);

    hold(ax, 'off');
end

function [x, y] = direction2stereo(dirs)
% DIRECTION2STEREO Convert multiple directions to stereographic coordinates

    n = size(dirs, 1);
    x = zeros(n, 1);
    y = zeros(n, 1);

    for i = 1:n
        [x(i), y(i)] = direction2stereoSingle(dirs(i,1), dirs(i,2), dirs(i,3));
    end
end

function [x, y] = direction2stereoSingle(h, k, l)
% DIRECTION2STEREOSINGLE Convert single direction to stereographic coordinates

    % Normalize
    norm_d = sqrt(h^2 + k^2 + l^2);
    if norm_d > 0
        h = h / norm_d;
        k = k / norm_d;
        l = l / norm_d;
    end

    % Stereographic projection
    theta = acos(l);
    if (h^2 + k^2) > 0
        phi = atan2(k, h);
    else
        phi = 0;
    end

    r = tan(theta / 2);
    x_stereo = r * cos(phi);
    y_stereo = r * sin(phi);

    % Swap for desired orientation
    x = y_stereo;
    y = x_stereo;
end

function [vertices, faces, colors] = createInterpolatedMesh(xData, yData, values, nGrid)
% CREATEINTERPOLATEDMESH Create interpolated mesh for smooth visualization

    if length(xData) < 3
        vertices = [];
        faces = [];
        colors = [];
        return;
    end

    % Create fine grid covering the triangle
    [x001, y001] = direction2stereoSingle(0, 0, 1);
    [x011, y011] = direction2stereoSingle(0, 1, 1);
    [x111, y111] = direction2stereoSingle(1, 1, 1);

    xMin = min([x001, x011, x111]) - 0.01;
    xMax = max([x001, x011, x111]) + 0.01;
    yMin = min([y001, y011, y111]) - 0.01;
    yMax = max([y001, y011, y111]) + 0.01;

    [xGrid, yGrid] = meshgrid(linspace(xMin, xMax, nGrid), linspace(yMin, yMax, nGrid));

    % Interpolate values onto grid
    try
        F = scatteredInterpolant(xData, yData, values, 'natural', 'none');
        vGrid = F(xGrid, yGrid);
    catch
        % Fallback to nearest neighbor
        vGrid = griddata(xData, yData, values, xGrid, yGrid, 'nearest');
    end

    % Mask points outside the triangle
    mask = isInsideTriangle(xGrid, yGrid);
    vGrid(~mask) = NaN;

    % Convert to patch format
    [faces, vertices, colors] = surf2patch(xGrid, yGrid, zeros(size(xGrid)), vGrid);

    % Remove faces with NaN colors
    validFaces = all(~isnan(colors(faces)), 2);
    faces = faces(validFaces, :);
end

function inside = isInsideTriangle(x, y)
% ISINSIDETRIANGLE Check if points are inside the standard triangle

    % Get triangle vertices in stereographic coordinates
    [x001, y001] = direction2stereoSingle(0, 0, 1);
    [x011, y011] = direction2stereoSingle(0, 1, 1);
    [x111, y111] = direction2stereoSingle(1, 1, 1);

    % Use barycentric coordinates
    v0 = [x111 - x001, y111 - y001];
    v1 = [x011 - x001, y011 - y001];

    inside = false(size(x));

    for i = 1:numel(x)
        v2 = [x(i) - x001, y(i) - y001];

        dot00 = dot(v0, v0);
        dot01 = dot(v0, v1);
        dot02 = dot(v0, v2);
        dot11 = dot(v1, v1);
        dot12 = dot(v1, v2);

        invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        u = (dot11 * dot02 - dot01 * dot12) * invDenom;
        v = (dot00 * dot12 - dot01 * dot02) * invDenom;

        inside(i) = (u >= -0.01) && (v >= -0.01) && (u + v <= 1.01);
    end
end
