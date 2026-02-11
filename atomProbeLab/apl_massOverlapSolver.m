
function res = apl_massOverlapSolver(pos, rangeTable, varargin)
% apl_massOverlapSolver performs bulk MLE deconvolution

p = inputParser;
addParameter(p, 'isotopeTable', []);
addParameter(p, 'minAbundance', 0.21);
addParameter(p, 'ionTable', []);
addParameter(p, 'massSpec', []);
addParameter(p, 'verbose', false);
parse(p, varargin{:});

OP = apl_overlapProblemBuilder(rangeTable, 'isotopeTable', p.Results.isotopeTable, 'minAbundance', p.Results.minAbundance, 'ionTable', p.Results.ionTable, 'massSpec', p.Results.massSpec);
OP = apl_overlapProblem_getRangeCounts(OP, pos, []);

F = OP.fractionMatrix;
rangeCounts = OP.rangeCounts;
ions = OP.ions;
numIons = height(ions);

ionCounts = zeros(numIons,1);

% non-overlapped ranges
for r = 1:height(rangeTable)
    ionIdx = find(F(r,:) > 0);
    if numel(ionIdx) == 1
        ionCounts(ionIdx) = ionCounts(ionIdx) + rangeCounts(r);
    end
end

% overlapped groups
gof = struct('group', {}, 'chi2', {}, 'rmse', {});
for g = 1:numel(OP.overlapRangeIdx)
    rIdx = OP.overlapRangeIdx{g};
    iIdx = OP.overlapIons{g};
    C = F(rIdx, iIdx);
    n = rangeCounts(rIdx);
    if sum(n) == 0
        continue;
    end
    % initial guess from lsqnonneg
    try
        p0_full = lsqnonneg(C, n);
        if sum(p0_full) > 0
            p0_full = p0_full / sum(p0_full);
        else
            p0_full = ones(numel(iIdx),1) / numel(iIdx);
        end
        p0 = p0_full(1:end-1)';
    catch
        p0 = [];
    end

    [p, ~] = apl_mleSolver(n, C, p0);
    ionCounts(iIdx) = p(:) * sum(n);

    fit = C * ionCounts(iIdx);
    chi2 = sum((n - fit).^2 ./ max(fit, 1));
    rmse = sqrt(mean((n - fit).^2));
    gof(end+1) = struct('group', g, 'chi2', chi2, 'rmse', rmse); %#ok<AGROW>
end

fitCounts = F * ionCounts;

res = struct();
res.ionTable = ions;
res.ionCounts = ionCounts;
res.ionFractions = ionCounts / max(sum(ionCounts), 1);
res.rangeCounts = rangeCounts;
res.fitCounts = fitCounts;
res.overlapProblem = OP;
res.gof = gof;
end
