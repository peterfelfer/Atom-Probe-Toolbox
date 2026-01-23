function stats = spatialStatistics(pos, options)
% SPATIALSTATISTICS Compute spatial distribution statistics for APT data
%
% stats = spatialStatistics(pos)
% stats = spatialStatistics(pos, 'rdfMaxR', 5)
%
% Calculates various spatial statistics including radial distribution
% function (RDF/pair correlation), nearest neighbor distribution, and
% K-function (Ripley's K) for analyzing clustering and ordering.
%
% INPUT:
%   pos - Nx3 array of atom positions [x, y, z] in nm
%         OR position table with x, y, z columns
%
% OPTIONS:
%   'rdfMaxR'      - Maximum radius for RDF in nm (default: 3)
%   'rdfBinWidth'  - Bin width for RDF in nm (default: 0.05)
%   'nnMaxK'       - Maximum k for k-th nearest neighbor (default: 10)
%   'sampleSize'   - Points to sample for large datasets (default: 50000)
%   'reference'    - Reference positions for cross-RDF (default: same as pos)
%   'edgeCorrection' - Apply edge correction (default: true)
%   'showProgress' - Show progress bar (default: true)
%   'computeRDF'   - Compute RDF (default: true)
%   'computeNN'    - Compute nearest neighbor (default: true)
%   'computeK'     - Compute Ripley's K (default: true)
%
% OUTPUT:
%   stats - Structure containing:
%       .rdf - Radial distribution function
%           .r        - Radii (bin centers)
%           .g        - g(r) values
%           .gError   - Standard error of g(r)
%           .nPairs   - Number of pairs per bin
%       .nn  - Nearest neighbor statistics
%           .k        - Neighbor order (1st, 2nd, ..., k-th)
%           .meanDist - Mean distance to k-th neighbor
%           .stdDist  - Standard deviation
%           .histogram - Distribution of distances for each k
%       .ripley - Ripley's K function
%           .r        - Radii
%           .K        - K(r) values
%           .L        - L(r) = sqrt(K(r)/pi) - r (linearized)
%           .Ktheory  - Theoretical K for random (CSR)
%       .density - Local density estimate
%       .nPoints - Number of points analyzed
%
% THEORY:
%   RDF g(r): Probability of finding an atom at distance r relative to
%   a uniform random distribution. g(r) = 1 for random, >1 for clustering,
%   <1 for ordering/repulsion.
%
%   Nearest neighbor: In a random (Poisson) distribution, the mean
%   nearest neighbor distance is d_nn = 0.554 * rho^(-1/3)
%
% EXAMPLES:
%   % Basic analysis
%   stats = spatialStatistics(pos);
%   plot(stats.rdf.r, stats.rdf.g);
%
%   % Cross-RDF between two species
%   stats = spatialStatistics(posFe, 'reference', posCu);
%
% REFERENCES:
%   Moody et al. (2009) Ultramicroscopy
%   Philippe et al. (2009) Ultramicroscopy
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos
    options.rdfMaxR (1,1) double {mustBePositive} = 3
    options.rdfBinWidth (1,1) double {mustBePositive} = 0.05
    options.nnMaxK (1,1) double {mustBePositive, mustBeInteger} = 10
    options.sampleSize (1,1) double {mustBePositive} = 50000
    options.reference = []
    options.edgeCorrection (1,1) logical = true
    options.showProgress (1,1) logical = true
    options.computeRDF (1,1) logical = true
    options.computeNN (1,1) logical = true
    options.computeK (1,1) logical = true
end

% Handle table input
if istable(pos)
    posArray = [pos.x, pos.y, pos.z];
else
    posArray = pos;
end

nPoints = size(posArray, 1);

% Handle reference positions (for cross-RDF)
if isempty(options.reference)
    refArray = posArray;
    crossCorrelation = false;
else
    if istable(options.reference)
        refArray = [options.reference.x, options.reference.y, options.reference.z];
    else
        refArray = options.reference;
    end
    crossCorrelation = true;
end

% Sample if dataset is large
if nPoints > options.sampleSize
    sampleIdx = randperm(nPoints, options.sampleSize);
    posAnalyze = posArray(sampleIdx, :);
    if options.showProgress
        fprintf('Sampling %d of %d points for analysis\n', options.sampleSize, nPoints);
    end
else
    posAnalyze = posArray;
end
nAnalyze = size(posAnalyze, 1);

% Calculate bounding box and density
minPos = min(posArray, [], 1);
maxPos = max(posArray, [], 1);
boxSize = maxPos - minPos;
volume = prod(boxSize);
density = nPoints / volume;

% Initialize output
stats = struct();
stats.nPoints = nPoints;
stats.density = density;
stats.volume = volume;
stats.boxSize = boxSize;

% Radial Distribution Function
if options.computeRDF
    if options.showProgress
        fprintf('Computing radial distribution function...\n');
    end
    stats.rdf = computeRDF(posAnalyze, refArray, density, options);
end

% Nearest Neighbor Distribution
if options.computeNN
    if options.showProgress
        fprintf('Computing nearest neighbor statistics...\n');
    end
    stats.nn = computeNearestNeighbor(posAnalyze, options);
end

% Ripley's K Function
if options.computeK
    if options.showProgress
        fprintf('Computing Ripley''s K function...\n');
    end
    stats.ripley = computeRipleyK(posAnalyze, density, volume, boxSize, options);
end

if options.showProgress
    fprintf('Spatial statistics complete.\n');
end

end

%% Radial Distribution Function
function rdf = computeRDF(pos, ref, density, options)
    % Compute g(r) using histogram of pairwise distances

    nPos = size(pos, 1);
    nRef = size(ref, 1);

    % Define bins
    edges = 0:options.rdfBinWidth:options.rdfMaxR;
    r = edges(1:end-1) + options.rdfBinWidth/2;  % Bin centers
    nBins = length(r);

    % Compute pairwise distances in chunks to manage memory
    chunkSize = min(1000, nPos);
    nChunks = ceil(nPos / chunkSize);

    counts = zeros(1, nBins);
    nPairsTotal = 0;

    for i = 1:nChunks
        startIdx = (i-1)*chunkSize + 1;
        endIdx = min(i*chunkSize, nPos);
        chunkPos = pos(startIdx:endIdx, :);

        % Distance to reference points
        D = pdist2(chunkPos, ref);

        % Exclude self-distances if same dataset
        if isequal(pos, ref)
            selfIdx = startIdx:endIdx;
            for j = 1:size(chunkPos, 1)
                D(j, selfIdx(j)) = NaN;
            end
        end

        % Apply maximum distance cutoff
        D(D > options.rdfMaxR) = NaN;

        % Histogram
        validD = D(~isnan(D));
        counts = counts + histcounts(validD, edges);
        nPairsTotal = nPairsTotal + numel(validD);
    end

    % Normalize to get g(r)
    % Expected count in shell: N * rho * 4*pi*r^2*dr
    shellVolumes = 4/3 * pi * (edges(2:end).^3 - edges(1:end-1).^3);
    expectedCounts = nPos * density * shellVolumes;

    g = counts ./ expectedCounts;

    % Edge correction (simple approach - reduce by fraction of shell inside box)
    if options.edgeCorrection
        % Approximate correction factor
        for i = 1:nBins
            % Fraction of shell that might be outside
            edgeFrac = min(1, min(options.boxSize) / (2*r(i)));
            g(i) = g(i) / edgeFrac;
        end
    end

    % Standard error (assuming Poisson statistics)
    gError = g ./ sqrt(max(counts, 1));

    rdf = struct();
    rdf.r = r;
    rdf.g = g;
    rdf.gError = gError;
    rdf.nPairs = counts;
    rdf.binWidth = options.rdfBinWidth;
end

%% Nearest Neighbor Statistics
function nn = computeNearestNeighbor(pos, options)
    % Compute k-th nearest neighbor distances

    nPos = size(pos, 1);
    maxK = min(options.nnMaxK, nPos - 1);

    % Use k-d tree for efficient neighbor search
    tree = KDTreeSearcher(pos);

    % Find k+1 nearest neighbors (includes self)
    [~, D] = knnsearch(tree, pos, 'K', maxK + 1);

    % Remove self-distance (first column)
    D = D(:, 2:end);

    % Statistics for each k
    nn = struct();
    nn.k = (1:maxK)';
    nn.meanDist = mean(D, 1)';
    nn.stdDist = std(D, 0, 1)';
    nn.medianDist = median(D, 1)';
    nn.minDist = min(D, [], 1)';
    nn.maxDist = max(D, [], 1)';

    % Distribution histograms for first few neighbors
    nn.histogram = struct();
    for k = 1:min(5, maxK)
        edges = linspace(0, max(D(:,k))*1.1, 50);
        nn.histogram(k).edges = edges;
        nn.histogram(k).counts = histcounts(D(:,k), edges);
        nn.histogram(k).k = k;
    end

    % Theoretical values for random (Poisson) distribution
    % Mean k-th NN distance: Gamma(k+1/3) / (Gamma(k) * (4*pi*rho/3)^(1/3))
    density = nPos / prod(range(pos));
    nn.theoretical = struct();
    nn.theoretical.density = density;
    nn.theoretical.meanFirst = 0.554 * density^(-1/3);  % Mean 1st NN for random

    % Compare to random
    nn.clusteringIndex = nn.theoretical.meanFirst ./ nn.meanDist(1);
    % >1 suggests clustering, <1 suggests ordering
end

%% Ripley's K Function
function ripley = computeRipleyK(pos, density, volume, boxSize, options)
    % Compute Ripley's K function

    nPos = size(pos, 1);

    % Radii to evaluate
    r = linspace(0.1, options.rdfMaxR, 30);
    K = zeros(size(r));

    % Build k-d tree
    tree = KDTreeSearcher(pos);

    % Count neighbors within each radius
    for i = 1:length(r)
        [~, D] = rangesearch(tree, pos, r(i));
        neighborCounts = cellfun(@length, D) - 1;  % Subtract self

        % Edge correction (Ripley's isotropic correction approximation)
        if options.edgeCorrection
            % Simple correction: weight by fraction inside box
            weights = edgeCorrectionWeights(pos, r(i), boxSize);
            K(i) = sum(neighborCounts ./ weights) / nPos / density;
        else
            K(i) = mean(neighborCounts) / density;
        end
    end

    % Theoretical K for Complete Spatial Randomness (CSR)
    Ktheory = 4/3 * pi * r.^3;

    % L function (variance-stabilized)
    L = (3 * K / (4*pi)).^(1/3) - r;

    ripley = struct();
    ripley.r = r;
    ripley.K = K;
    ripley.L = L;
    ripley.Ktheory = Ktheory;

    % Deviation from CSR
    ripley.Lmax = max(abs(L));
    ripley.clusteringDetected = ripley.Lmax > 0.5;  % Rule of thumb threshold
end

function weights = edgeCorrectionWeights(pos, r, boxSize)
    % Compute edge correction weights for each point
    % Weight = 1 / (fraction of sphere inside box)

    nPos = size(pos, 1);
    weights = ones(nPos, 1);

    for i = 1:nPos
        % Distance to each face
        distToMin = pos(i, :);
        distToMax = boxSize - pos(i, :);

        % Minimum distance to any edge
        minDist = min([distToMin, distToMax]);

        if minDist < r
            % Approximate: fraction of sphere inside
            % Simple approximation using cap formula
            frac = 1;
            for d = 1:3
                if distToMin(d) < r
                    h = r - distToMin(d);
                    capFrac = (h^2 * (3*r - h)) / (4*r^3);
                    frac = frac * (1 - capFrac);
                end
                if distToMax(d) < r
                    h = r - distToMax(d);
                    capFrac = (h^2 * (3*r - h)) / (4*r^3);
                    frac = frac * (1 - capFrac);
                end
            end
            weights(i) = max(frac, 0.1);  % Avoid division by very small numbers
        end
    end
end

function r = range(x)
    r = max(x) - min(x);
end
