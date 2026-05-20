function [A, terms] = rhitDesignXY(x, y, maxDeg)
% RHITDESIGNXY Design matrix for a 2D polynomial in (x, y).
%
% [A, terms] = rhitDesignXY(x, y, maxDeg)
%
% Returns the design matrix `A` whose columns are x.^i .* y.^j for every
% (i, j) with i + j <= maxDeg, listed in the same order as the matching
% `terms` (Ncols-by-2) array.
%
% INPUTS:
%   x, y   - column vectors of detector coordinates (mm)
%   maxDeg - polynomial degree (commonly 4 for a LEAP bowl)
%
% OUTPUTS:
%   A     - (numel(x), Ncols) design matrix
%   terms - (Ncols, 2) integer (i, j) exponent pairs
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    x (:,1) double
    y (:,1) double
    maxDeg (1,1) double {mustBeInteger, mustBeNonnegative}
end

if numel(x) ~= numel(y)
    error('rhitDesignXY:sizeMismatch', 'x and y must be the same length');
end

terms = zeros(0, 2);
for i = 0:maxDeg
    for j = 0:(maxDeg - i)
        terms(end+1, :) = [i, j]; %#ok<AGROW>
    end
end

n = numel(x);
A = zeros(n, size(terms, 1));
for k = 1:size(terms, 1)
    A(:, k) = (x .^ terms(k, 1)) .* (y .^ terms(k, 2));
end
end
