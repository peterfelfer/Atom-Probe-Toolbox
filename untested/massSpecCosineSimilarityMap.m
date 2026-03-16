function [simMap, info] = massSpecCosineSimilarityMap(pos, options)
% MASSSPECCOSINESIMLIARITYMAP Voxelised cosine similarity of local mass spectra.
%
% [simMap, info] = massSpecCosineSimilarityMap(pos)
% [simMap, info] = massSpecCosineSimilarityMap(pos, 'voxelSize', 2)
%
% Divides the dataset into a 3D voxel grid. For each voxel, computes the
% mass spectrum histogram and its cosine similarity to the global
% (whole-dataset) spectrum. Regions with different chemistry (e.g.
% precipitates, grain boundaries) appear as low-similarity voxels.
%
% INPUT
%   pos - pos table with x, y, z, mc columns
%
% OPTIONS
%   'voxelSize'  - voxel edge length in nm (default: 2)
%   'binWidth'   - mass spectrum bin width in Da (default: 0.1)
%   'mcMax'      - maximum m/c to include in Da (default: 120)
%   'minAtoms'   - minimum atoms per voxel to compute similarity (default: 20)
%   'reference'  - reference spectrum: 'global' (default) or a numeric
%                  vector (custom reference histogram)
%
% OUTPUT
%   simMap - 3D array (nx × ny × nz) of bias-corrected cosine similarity.
%            Corrected for the expected similarity of a random N-atom
%            sample (Poisson noise), so values are comparable across
%            voxels with different atom counts. NaN for voxels with fewer
%            than minAtoms atoms.
%
%   info   - struct with:
%       simMapRaw    - uncorrected cosine similarity (before bias correction)
%       countMap     - 3D array of atom counts per voxel
%       xEdges, yEdges, zEdges - voxel grid edges (nm)
%       globalHist   - the reference spectrum histogram
%       mcEdges      - mass spectrum bin edges (Da)
%       voxelSize    - voxel size used (nm)
%       binWidth     - bin width used (Da)
%       biasModel    - struct with a, b: E[cos_sim] = 1 - a/N^b
%
% EXAMPLE
%   pos = posLoad('dataset.epos');
%   [simMap, info] = massSpecCosineSimilarityMap(pos, 'voxelSize', 2);
%
%   % Visualise
%   [ix,iy,iz] = ind2sub(size(simMap), find(~isnan(simMap)));
%   xc = info.xEdges(ix) + info.voxelSize/2;
%   yc = info.yEdges(iy) + info.voxelSize/2;
%   zc = info.zEdges(iz) + info.voxelSize/2;
%   scatter3(xc, yc, zc, 8, simMap(~isnan(simMap)), 'filled');
%   colorbar; colormap(flipud(hot)); clim([0.9 1]);
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    pos table
    options.voxelSize (1,1) double {mustBePositive} = 2
    options.binWidth (1,1) double {mustBePositive} = 0.1
    options.mcMax (1,1) double {mustBePositive} = 120
    options.minAtoms (1,1) double {mustBePositive} = 20
    options.reference = 'global'
end

%% Setup
voxelSize = options.voxelSize;
binWidth = options.binWidth;
mcMax = options.mcMax;

mc = pos.mc;
mcEdges = 0:binWidth:mcMax;
nBins = numel(mcEdges) - 1;

%% Voxel grid
xEdges = min(pos.x):voxelSize:(max(pos.x) + voxelSize);
yEdges = min(pos.y):voxelSize:(max(pos.y) + voxelSize);
zEdges = min(pos.z):voxelSize:(max(pos.z) + voxelSize);
nx = numel(xEdges) - 1;
ny = numel(yEdges) - 1;
nz = numel(zEdges) - 1;

%% Assign atoms to voxels and mc bins
xi = discretize(pos.x, xEdges);
yi = discretize(pos.y, yEdges);
zi = discretize(pos.z, zEdges);
mi = discretize(mc, mcEdges);

valid = ~isnan(xi) & ~isnan(yi) & ~isnan(zi) & ~isnan(mi);
xi = xi(valid);
yi = yi(valid);
zi = zi(valid);
mi = mi(valid);

%% Reference spectrum
if ischar(options.reference) || isstring(options.reference)
    globalHist = accumarray(mi, 1, [nBins 1])';
else
    globalHist = options.reference(:)';
    if numel(globalHist) ~= nBins
        error('massSpecCosineSimilarityMap:refSize', ...
            'Reference histogram length (%d) must match number of bins (%d).', ...
            numel(globalHist), nBins);
    end
end

%% Build per-voxel histograms using sparse accumarray
% Linear voxel index for each atom
voxelIdx = sub2ind([nx ny nz], xi, yi, zi);
nVoxels = nx * ny * nz;

% Build a sparse (nVoxels × nBins) matrix: each row is one voxel's spectrum
% Using accumarray with a 2D subs: [voxelIdx, mcBinIdx] -> count
specMatrix = accumarray([voxelIdx, mi], 1, [nVoxels, nBins], [], [], true);

%% Atom counts per voxel
countVec = full(sum(specMatrix, 2));  % nVoxels × 1

%% Cosine similarity (vectorised)
% cos_sim = (A . B) / (|A| * |B|)
% A = each row of specMatrix, B = globalHist

globalNorm = norm(globalHist);
dotProducts = specMatrix * globalHist';       % nVoxels × 1
localNorms = sqrt(sum(specMatrix.^2, 2));     % nVoxels × 1

% Compute similarity (avoid division by zero)
simVec = NaN(nVoxels, 1);
validVoxel = countVec >= options.minAtoms & localNorms > 0;
simVec(validVoxel) = dotProducts(validVoxel) ./ (localNorms(validVoxel) * globalNorm);

%% Correct for atom-count bias
% A random draw of N atoms from the global distribution has an expected
% cosine similarity < 1 due to Poisson noise. This biases low-count voxels
% toward lower similarity. Correct by dividing by E[cos_sim(N)].
%
% The expected similarity follows: E[cos_sim] ≈ 1 - a / N^b
% where a and b are estimated by Monte Carlo sampling from the global
% distribution. This is dataset-specific because it depends on the shape
% of the global spectrum.
globalProb = globalHist / sum(globalHist);
nTrials = 200;
sampleN = [20 50 100 200 500 1000 5000];
expectedSim = zeros(size(sampleN));
for i = 1:numel(sampleN)
    sims = zeros(nTrials, 1);
    for t = 1:nTrials
        sample = mnrnd(sampleN(i), globalProb);
        sims(t) = dot(sample, globalHist) / (norm(sample) * globalNorm);
    end
    expectedSim(i) = mean(sims);
end

% Fit model: 1 - E[cos_sim] = a * N^(-b)  in log-log space
logDev = log(max(1 - expectedSim, eps));
logN = log(sampleN);
pFit = polyfit(logN, logDev, 1);
biasFitB = -pFit(1);
biasFitA = exp(pFit(2));

expectedSimVec = 1 - biasFitA ./ max(countVec, 1).^biasFitB;
expectedSimVec = max(expectedSimVec, 0.01);  % avoid division by near-zero

simVecCorrected = simVec ./ expectedSimVec;

%% Reshape to 3D
simMap = reshape(simVecCorrected, [nx ny nz]);
simMapRaw = reshape(simVec, [nx ny nz]);
countMap = reshape(countVec, [nx ny nz]);

%% Info struct
info = struct();
info.simMapRaw = simMapRaw;
info.countMap = countMap;
info.xEdges = xEdges;
info.yEdges = yEdges;
info.zEdges = zEdges;
info.globalHist = globalHist;
info.mcEdges = mcEdges;
info.voxelSize = voxelSize;
info.binWidth = binWidth;
info.biasModel = struct('a', biasFitA, 'b', biasFitB, ...
    'formula', 'E[cos_sim] = 1 - a / N^b');

end
