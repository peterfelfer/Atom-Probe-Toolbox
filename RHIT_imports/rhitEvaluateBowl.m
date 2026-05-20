function C = rhitEvaluateBowl(bowl, x, y)
% RHITEVALUATEBOWL Evaluate a calibration bowl polynomial at (x, y).
%
% C = rhitEvaluateBowl(bowl, x, y)
%
% INPUTS:
%   bowl - struct with fields:
%            kind    : 'xy2d' or 'radial'
%            degree  : polynomial degree
%            coeffs  : column vector of coefficients
%            terms   : exponent table
%                      'xy2d'   -> Ncols-by-2 (i, j) for x^i*y^j
%                      'radial' -> Ncols-by-1 [k] for (r^2)^k
%   x, y - detector coordinates (mm), same shape
%
% OUTPUT:
%   C    - bowl value at each (x, y), same shape as x
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    bowl struct
    x double
    y double
end

if ~isequal(size(x), size(y))
    error('rhitEvaluateBowl:sizeMismatch', 'x and y must be the same shape');
end

origShape = size(x);
xc = double(x(:));
yc = double(y(:));

switch lower(bowl.kind)
    case 'xy2d'
        A = rhitDesignXY(xc, yc, bowl.degree);
    case 'radial'
        r2 = xc .^ 2 + yc .^ 2;
        A = zeros(numel(xc), bowl.degree + 1);
        for k = 0:bowl.degree
            A(:, k + 1) = r2 .^ k;
        end
    otherwise
        error('rhitEvaluateBowl:unknownKind', 'Unknown bowl kind "%s"', bowl.kind);
end

C = reshape(A * bowl.coeffs(:), origShape);
end
