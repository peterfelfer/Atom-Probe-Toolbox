
function [p, exitflag] = apl_mleSolver(n, C, p0)
% apl_mleSolver MLE for multinomial with constrained composition

if nargin < 3 || isempty(p0)
    p0 = ones(1, size(C,2)-1) / size(C,2);
end

% Normalize C columns
C = bsxfun(@rdivide, C, sum(C,1));

n = n(:)';
if sum(n) <= 0 || all(n == 0)
    p = ones(1, size(C,2)) / size(C,2);
    exitflag = 0;
    return;
end

p0 = max(p0, 1e-8);

% Logistic transform for bounds
u0 = log(p0 ./ max(1 - p0, 1e-8));

L = @(u) negloglik(u, n, C);
opts = optimset('MaxFunEvals', 2000, 'MaxIter', 2000, 'Display','off');

try
    [u, ~, exitflag] = fminsearch(L, u0, opts);
catch
    u = u0;
    exitflag = -1;
end

praw = 1 ./ (1 + exp(-u));
p = compNorm(praw);
end

function val = negloglik(u, n, C)
    p = 1 ./ (1 + exp(-u));
    p = compNorm(p);
    q = p * C';
    q(q <= 0) = eps;
    val = -sum(n .* log(q));
end

function p1 = compNorm(p0)
    t = sum(p0);
    p1 = [p0, max(0, 1 - t)];
    t = sum(p1);
    if t == 0
        p1 = ones(1, numel(p1)) / numel(p1);
    else
        p1 = p1 ./ t;
    end
end
