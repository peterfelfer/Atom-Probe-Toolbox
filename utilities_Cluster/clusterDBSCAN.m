function [clusterIdx, clusterInfo] = clusterDBSCAN(pos, epsilon, minPts, options)
% CLUSTERDBSCAN Density-based spatial clustering for APT data
%
% [clusterIdx, clusterInfo] = clusterDBSCAN(pos, epsilon, minPts)
% [clusterIdx, clusterInfo] = clusterDBSCAN(pos, epsilon, minPts, 'useGPU', true)
%
% Performs DBSCAN (Density-Based Spatial Clustering of Applications with
% Noise) clustering on atom probe position data. Automatically uses GPU
% acceleration when available with CPU fallback.
%
% INPUT:
%   pos      - Nx3 array of atom positions [x, y, z] in nm
%              OR position table with x, y, z columns
%   epsilon  - Neighborhood radius in nm (typical: 0.3-1.0 nm for APT)
%   minPts   - Minimum points to form a cluster (typical: 5-20)
%
% OPTIONS:
%   'useGPU'      - Try to use GPU acceleration (default: true)
%   'useParallel' - Use parallel processing for CPU (default: true)
%   'soluteIdx'   - Logical or index vector of solute atoms to cluster
%                   (default: all atoms)
%   'chunkSize'   - Points per chunk for large datasets (default: 50000)
%   'showProgress'- Show progress indicator (default: true)
%   'algorithm'   - 'kdtree' or 'bruteforce' (default: 'kdtree')
%
% OUTPUT:
%   clusterIdx  - Nx1 array of cluster assignments
%                 -1 = noise, 0 = not analyzed, 1,2,3... = cluster ID
%   clusterInfo - Structure with cluster statistics:
%       .nClusters      - Number of clusters found
%       .nNoise         - Number of noise points
%       .nClustered     - Number of clustered points
%       .clusterSizes   - Number of atoms per cluster
%       .clusterCenters - Centroid of each cluster [nClusters x 3]
%       .clusterRadii   - Radius of gyration per cluster
%       .computeTime    - Time taken for clustering
%       .method         - 'GPU' or 'CPU'
%
% ALGORITHM:
%   DBSCAN groups points that are closely packed together, marking points
%   in low-density regions as outliers. Unlike k-means, it doesn't require
%   specifying the number of clusters and can find arbitrarily shaped clusters.
%
% EXAMPLES:
%   % Basic clustering of solute atoms
%   [idx, info] = clusterDBSCAN(pos, 0.5, 10);
%
%   % Cluster only specific solute atoms (e.g., Cu in Fe matrix)
%   soluteIdx = strcmp(posTable.ion, 'Cu');
%   [idx, info] = clusterDBSCAN(posTable, 0.6, 5, 'soluteIdx', soluteIdx);
%
% REFERENCES:
%   Ester et al. (1996) "A density-based algorithm for discovering clusters"
%   Marquis & Hyde (2010) "Applications of atom-probe tomography to the
%   characterisation of solute behaviours"
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos
    epsilon (1,1) double {mustBePositive}
    minPts (1,1) double {mustBePositive, mustBeInteger}
    options.useGPU (1,1) logical = true
    options.useParallel (1,1) logical = true
    options.soluteIdx = []
    options.chunkSize (1,1) double {mustBePositive} = 50000
    options.showProgress (1,1) logical = true
    options.algorithm (1,:) char {mustBeMember(options.algorithm, {'kdtree', 'bruteforce'})} = 'kdtree'
end

startTime = tic;

% Handle table input
if istable(pos)
    posArray = [pos.x, pos.y, pos.z];
else
    posArray = pos;
end

nPoints = size(posArray, 1);

% Handle solute selection
if ~isempty(options.soluteIdx)
    if islogical(options.soluteIdx)
        soluteIdx = options.soluteIdx;
    else
        soluteIdx = false(nPoints, 1);
        soluteIdx(options.soluteIdx) = true;
    end
    analyzeIdx = find(soluteIdx);
    posAnalyze = posArray(soluteIdx, :);
else
    analyzeIdx = (1:nPoints)';
    posAnalyze = posArray;
end

nAnalyze = size(posAnalyze, 1);

% Initialize output
clusterIdx = zeros(nPoints, 1);  % 0 = not analyzed

% Check GPU availability
useGPU = options.useGPU && isGPUAvailable();

if useGPU
    try
        [labels, method] = dbscanGPU(posAnalyze, epsilon, minPts, options);
    catch ME
        warning('clusterDBSCAN:gpuFailed', ...
            'GPU clustering failed: %s. Falling back to CPU.', ME.message);
        useGPU = false;
    end
end

if ~useGPU
    [labels, method] = dbscanCPU(posAnalyze, epsilon, minPts, options);
end

% Map results back to full position array
clusterIdx(analyzeIdx) = labels;

% Calculate cluster statistics
computeTime = toc(startTime);
clusterInfo = calculateClusterStats(posArray, clusterIdx, computeTime, method);

if options.showProgress
    fprintf('DBSCAN complete: %d clusters found, %d noise points (%.2f s, %s)\n', ...
        clusterInfo.nClusters, clusterInfo.nNoise, computeTime, method);
end

end

%% GPU Implementation
function [labels, method] = dbscanGPU(pos, epsilon, minPts, options)
    % GPU-accelerated DBSCAN using distance matrix computation
    method = 'GPU';

    nPoints = size(pos, 1);

    % Transfer to GPU
    posGPU = gpuArray(single(pos));

    % For large datasets, process in chunks
    if nPoints > options.chunkSize
        labels = dbscanChunkedGPU(posGPU, epsilon, minPts, options);
    else
        % Compute distance matrix on GPU
        D = pdist2GPU(posGPU, posGPU);

        % Find neighbors within epsilon
        neighbors = D <= epsilon;

        % Perform DBSCAN on GPU
        labels = gather(dbscanCore(neighbors, minPts));
    end
end

function D = pdist2GPU(X, Y)
    % Compute pairwise Euclidean distance on GPU
    % D(i,j) = ||X(i,:) - Y(j,:)||

    X2 = sum(X.^2, 2);
    Y2 = sum(Y.^2, 2)';
    D = sqrt(max(X2 + Y2 - 2 * (X * Y'), 0));
end

function labels = dbscanChunkedGPU(posGPU, epsilon, minPts, options)
    % Process large datasets in chunks on GPU

    nPoints = size(posGPU, 1);
    nChunks = ceil(nPoints / options.chunkSize);

    % First pass: find core points using chunked distance computation
    corePoints = false(nPoints, 1);

    if options.showProgress
        fprintf('Pass 1: Identifying core points...\n');
    end

    for i = 1:nChunks
        startIdx = (i-1) * options.chunkSize + 1;
        endIdx = min(i * options.chunkSize, nPoints);

        chunkPos = posGPU(startIdx:endIdx, :);

        % Count neighbors for this chunk against all points
        neighborCount = zeros(endIdx - startIdx + 1, 1, 'gpuArray');

        for j = 1:nChunks
            jStart = (j-1) * options.chunkSize + 1;
            jEnd = min(j * options.chunkSize, nPoints);

            D = pdist2GPU(chunkPos, posGPU(jStart:jEnd, :));
            neighborCount = neighborCount + sum(D <= epsilon, 2);
        end

        corePoints(startIdx:endIdx) = gather(neighborCount) >= minPts;
    end

    % Fall back to CPU for cluster assignment (complex graph operation)
    labels = assignClustersCPU(gather(posGPU), epsilon, corePoints, options);
end

%% CPU Implementation
function [labels, method] = dbscanCPU(pos, epsilon, minPts, options)
    method = 'CPU';

    nPoints = size(pos, 1);

    % Use MATLAB's built-in dbscan if available (R2019a+)
    if exist('dbscan', 'file') == 2
        labels = dbscan(pos, epsilon, minPts);
        return;
    end

    % Custom implementation using k-d tree
    if strcmp(options.algorithm, 'kdtree') && nPoints > 1000
        labels = dbscanKDTree(pos, epsilon, minPts, options);
    else
        labels = dbscanBruteForce(pos, epsilon, minPts, options);
    end
end

function labels = dbscanKDTree(pos, epsilon, minPts, options)
    % DBSCAN using k-d tree for efficient neighbor queries

    nPoints = size(pos, 1);

    % Build k-d tree
    tree = KDTreeSearcher(pos);

    % Find neighbors for all points
    neighborIdx = rangesearch(tree, pos, epsilon);

    % Identify core points
    neighborCounts = cellfun(@length, neighborIdx);
    corePoints = neighborCounts >= minPts;

    % Assign clusters
    labels = assignClustersCPU(pos, epsilon, corePoints, options, neighborIdx);
end

function labels = dbscanBruteForce(pos, epsilon, minPts, options)
    % Brute force DBSCAN for small datasets

    nPoints = size(pos, 1);

    % Compute full distance matrix
    D = pdist2(pos, pos);

    % Find neighbors
    neighbors = D <= epsilon;
    neighborCounts = sum(neighbors, 2);

    % Core points
    corePoints = neighborCounts >= minPts;

    % Initialize labels
    labels = zeros(nPoints, 1);
    clusterID = 0;
    visited = false(nPoints, 1);

    for i = 1:nPoints
        if visited(i) || ~corePoints(i)
            continue;
        end

        clusterID = clusterID + 1;
        [labels, visited] = expandCluster(i, neighbors, corePoints, labels, visited, clusterID);
    end

    % Mark remaining unvisited points as noise
    labels(~visited & labels == 0) = -1;
end

function labels = assignClustersCPU(pos, epsilon, corePoints, options, neighborIdx)
    % Assign clusters using pre-computed core points

    nPoints = size(pos, 1);

    if nargin < 5
        % Need to compute neighbors
        tree = KDTreeSearcher(pos);
        neighborIdx = rangesearch(tree, pos, epsilon);
    end

    % Initialize
    labels = zeros(nPoints, 1);
    clusterID = 0;
    visited = false(nPoints, 1);

    corePointsIdx = find(corePoints);

    if options.showProgress && length(corePointsIdx) > 1000
        prog = ProgressTracker(length(corePointsIdx), 'Assigning clusters');
    end

    processedCores = 0;
    for i = corePointsIdx'
        if visited(i)
            continue;
        end

        clusterID = clusterID + 1;

        % BFS to find all connected points
        queue = i;
        visited(i) = true;
        labels(i) = clusterID;

        while ~isempty(queue)
            current = queue(1);
            queue(1) = [];

            currentNeighbors = neighborIdx{current};

            for j = currentNeighbors
                if ~visited(j)
                    visited(j) = true;
                    labels(j) = clusterID;

                    if corePoints(j)
                        queue(end+1) = j;
                    end
                end
            end
        end

        processedCores = processedCores + 1;
        if options.showProgress && length(corePointsIdx) > 1000 && mod(processedCores, 100) == 0
            prog.update(processedCores);
        end
    end

    if options.showProgress && length(corePointsIdx) > 1000
        prog.finish();
    end

    % Mark noise
    labels(labels == 0) = -1;
end

function [labels, visited] = expandCluster(pointIdx, neighbors, corePoints, labels, visited, clusterID)
    % Expand cluster from a core point

    queue = pointIdx;
    visited(pointIdx) = true;
    labels(pointIdx) = clusterID;

    while ~isempty(queue)
        current = queue(1);
        queue(1) = [];

        currentNeighbors = find(neighbors(current, :));

        for j = currentNeighbors
            if ~visited(j)
                visited(j) = true;
                labels(j) = clusterID;

                if corePoints(j)
                    queue(end+1) = j;
                end
            end
        end
    end
end

function labels = dbscanCore(neighbors, minPts)
    % Core DBSCAN algorithm for GPU arrays

    nPoints = size(neighbors, 1);

    % Count neighbors
    neighborCounts = sum(neighbors, 2);
    corePoints = neighborCounts >= minPts;

    % Initialize labels on GPU
    labels = zeros(nPoints, 1, 'gpuArray');
    clusterID = 0;

    visited = false(nPoints, 1, 'gpuArray');

    for i = 1:nPoints
        if visited(i) || ~corePoints(i)
            continue;
        end

        clusterID = clusterID + 1;

        % Simple expansion (not optimal for GPU but works)
        stack = i;
        while ~isempty(stack)
            current = stack(end);
            stack(end) = [];

            if visited(current)
                continue;
            end

            visited(current) = true;
            labels(current) = clusterID;

            if corePoints(current)
                newNeighbors = find(gather(neighbors(current, :)) & ~gather(visited'));
                stack = [stack; newNeighbors(:)];
            end
        end
    end

    labels(~visited) = -1;
end

%% Helper Functions
function tf = isGPUAvailable()
    % Check if GPU is available for computation
    tf = false;
    try
        if exist('gpuDeviceCount', 'file') && gpuDeviceCount > 0
            gpu = gpuDevice;
            % Check for sufficient memory (at least 1 GB)
            if gpu.AvailableMemory > 1e9
                tf = true;
            end
        end
    catch
        tf = false;
    end
end

function clusterInfo = calculateClusterStats(pos, clusterIdx, computeTime, method)
    % Calculate cluster statistics

    clusterInfo = struct();
    clusterInfo.computeTime = computeTime;
    clusterInfo.method = method;

    % Basic counts
    uniqueClusters = unique(clusterIdx);
    uniqueClusters = uniqueClusters(uniqueClusters > 0);  % Exclude noise (-1) and unanalyzed (0)

    clusterInfo.nClusters = length(uniqueClusters);
    clusterInfo.nNoise = sum(clusterIdx == -1);
    clusterInfo.nClustered = sum(clusterIdx > 0);
    clusterInfo.nUnanalyzed = sum(clusterIdx == 0);

    % Per-cluster statistics
    if clusterInfo.nClusters > 0
        clusterInfo.clusterSizes = zeros(clusterInfo.nClusters, 1);
        clusterInfo.clusterCenters = zeros(clusterInfo.nClusters, 3);
        clusterInfo.clusterRadii = zeros(clusterInfo.nClusters, 1);

        for i = 1:clusterInfo.nClusters
            cID = uniqueClusters(i);
            mask = clusterIdx == cID;
            clusterPos = pos(mask, :);

            clusterInfo.clusterSizes(i) = sum(mask);
            clusterInfo.clusterCenters(i, :) = mean(clusterPos, 1);

            % Radius of gyration
            centered = clusterPos - clusterInfo.clusterCenters(i, :);
            clusterInfo.clusterRadii(i) = sqrt(mean(sum(centered.^2, 2)));
        end
    else
        clusterInfo.clusterSizes = [];
        clusterInfo.clusterCenters = [];
        clusterInfo.clusterRadii = [];
    end
end
