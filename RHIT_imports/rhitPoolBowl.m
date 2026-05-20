function pooledBowl = rhitPoolBowl(calibStructs, options)
% RHITPOOLBOWL Pool per-run 2D bowls into one global instrument bowl.
%
% pooledBowl = rhitPoolBowl({calib1, calib2, ...})
% pooledBowl = rhitPoolBowl(calibStructs, 'gridStep', 0.7, 'bowlDegree', 4)
%
% Strategy:
%   1. Sample each per-run bowl on a common (x, y) grid in absolute units.
%   2. Take the per-cell median across runs (robust to outliers).
%   3. Refit a single 2D polynomial of the requested degree to the
%      median grid.  The result is the global bowl.
%
% This works because, on a single LEAP, the bowl is stable to ~0.1% in
% C0 across runs (the entire variation is in t_offset).  Pooling a
% handful of trustworthy runs gives a calibration that reproduces EPOS
% mc to within 0.04 Da on the strong peaks for runs of any composition.
%
% INPUTS:
%   calibStructs - cell array of calibration structs from
%                  rhitCalibrateFromEpos.  Each must have a `.bowl`
%                  sub-struct of kind 'xy2d' (legacy radial structs are
%                  silently re-evaluated on the grid; that's fine).
%
% NAME-VALUE OPTIONS:
%   'gridStep'    - mm grid spacing (default 0.7)
%   'gridHalf'    - half-extent (default 17 mm)
%   'bowlDegree'  - degree of the refit polynomial (default 4)
%
% OUTPUT:
%   pooledBowl  - struct with kind='xy2d', degree, coeffs, terms,
%                 plus extras: pooledFrom, C0Global, fitResidualPctC0
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    calibStructs cell
    options.gridStep (1,1) double = 0.7
    options.gridHalf (1,1) double = 17
    options.bowlDegree (1,1) double {mustBeInteger, mustBePositive} = 4
end

if numel(calibStructs) < 2
    error('rhitPoolBowl:tooFew', 'Need at least 2 calibrations to pool.');
end

% Build grid
gx = -options.gridHalf:options.gridStep:options.gridHalf;
[X, Y] = meshgrid(gx, gx);
xs = X(:);
ys = Y(:);

% Sample each bowl on the grid
S = zeros(numel(calibStructs), numel(xs));
for k = 1:numel(calibStructs)
    if isfield(calibStructs{k}, 'bowl')
        bowl = calibStructs{k}.bowl;
    elseif isfield(calibStructs{k}, 'C_poly')
        bowl = struct('kind', 'radial', ...
            'degree', max(1, numel(calibStructs{k}.C_poly) - 1), ...
            'coeffs', calibStructs{k}.C_poly(:), ...
            'terms', (0:numel(calibStructs{k}.C_poly)-1)');
    else
        error('rhitPoolBowl:badStruct', ...
            'Entry %d has neither .bowl nor .C_poly', k);
    end
    S(k, :) = rhitEvaluateBowl(bowl, xs, ys)';
end

C_med = median(S, 1)';
C_std = std(S, 0, 1)';

% Refit polynomial to the median grid
[A, terms] = rhitDesignXY(xs, ys, options.bowlDegree);
coeffs = A \ C_med;
residStd = std(C_med - A * coeffs);

C0 = coeffs(1);

pooledBowl = struct();
pooledBowl.kind = 'xy2d';
pooledBowl.degree = options.bowlDegree;
pooledBowl.coeffs = coeffs;
pooledBowl.terms = terms;
pooledBowl.C0Global = C0;
pooledBowl.fitResidualPctC0 = 100 * residStd / C0;
pooledBowl.gridStdMeanPct = 100 * mean(C_std ./ C_med);
pooledBowl.gridStdMaxPct = 100 * max(C_std ./ C_med);
pooledBowl.nPooled = numel(calibStructs);
pooledBowl.gridStep = options.gridStep;
pooledBowl.gridHalf = options.gridHalf;

fprintf('Pooled %d bowls:\n', numel(calibStructs));
fprintf('  C(0,0)              = %.4e\n', C0);
fprintf('  per-cell std/median = mean %.4f%%, max %.4f%%\n', ...
    pooledBowl.gridStdMeanPct, pooledBowl.gridStdMaxPct);
fprintf('  refit residual      = %.4e (%.4f%% of C0)\n', ...
    residStd, pooledBowl.fitResidualPctC0);
end
