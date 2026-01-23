function [uncertainty, details] = concentrationUncertainty(counts, options)
% CONCENTRATIONUNCERTAINTY Calculate concentration uncertainty from counting statistics
%
% [uncertainty, details] = concentrationUncertainty(counts)
% [uncertainty, details] = concentrationUncertainty(counts, 'method', 'binomial')
%
% Calculates the uncertainty in composition measurements based on counting
% statistics. Several methods are available depending on the analysis context.
%
% INPUT:
%   counts - Vector or matrix of atom counts per element/ion
%            Each row represents a measurement, columns are elements
%            OR a table with ion counts
%
% OPTIONS:
%   'method'    - Uncertainty calculation method (default: 'binomial')
%                 'poisson'  - sqrt(N) counting statistics
%                 'binomial' - Full binomial uncertainty (recommended)
%                 'clopper'  - Clopper-Pearson exact confidence interval
%   'confidence' - Confidence level for intervals (default: 0.95)
%   'detectionEfficiency' - Detector efficiency for correction (default: 1)
%   'background' - Background counts to subtract (default: 0)
%
% OUTPUT:
%   uncertainty - Structure containing:
%       .concentration  - Calculated concentrations (atomic fraction)
%       .sigma          - Standard deviation of concentration
%       .ciLower        - Lower confidence interval
%       .ciUpper        - Upper confidence interval
%       .totalAtoms     - Total atom count per measurement
%
%   details - Additional calculation details
%
% THEORY:
%   For binomial statistics, the variance of concentration c = n/N is:
%   Var(c) = c(1-c)/N
%
%   This accounts for the constraint that concentrations must sum to 1.
%
% EXAMPLES:
%   % Single measurement with 3 elements
%   counts = [150, 40, 10];  % Fe, Cr, C atoms
%   unc = concentrationUncertainty(counts);
%
%   % Multiple measurements
%   counts = [150 40 10; 148 42 10; 152 38 10];
%   unc = concentrationUncertainty(counts, 'confidence', 0.99);
%
% REFERENCES:
%   Danoix, F. et al. (2007) Microscopy and Microanalysis
%   Miller, M.K. (2000) Atom Probe Tomography
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    counts {mustBeNumeric}
    options.method (1,:) char {mustBeMember(options.method, {'poisson', 'binomial', 'clopper'})} = 'binomial'
    options.confidence (1,1) double {mustBeInRange(options.confidence, 0, 1)} = 0.95
    options.detectionEfficiency (1,1) double {mustBePositive} = 1
    options.background = 0
end

% Handle table input
if istable(counts)
    counts = table2array(counts);
end

% Ensure counts is 2D (rows = measurements, columns = elements)
if isvector(counts)
    counts = counts(:)';  % Row vector
end

[nMeasurements, nElements] = size(counts);

% Subtract background if provided
if ~isscalar(options.background) || options.background ~= 0
    if isscalar(options.background)
        counts = counts - options.background;
    else
        counts = counts - options.background(:)';
    end
    counts = max(counts, 0);  % Ensure non-negative
end

% Calculate total atoms per measurement
totalAtoms = sum(counts, 2);

% Correct for detection efficiency
correctedCounts = counts / options.detectionEfficiency;
correctedTotal = totalAtoms / options.detectionEfficiency;

% Calculate concentrations
concentration = counts ./ totalAtoms;
concentration(isnan(concentration)) = 0;

% Initialize output
uncertainty = struct();
uncertainty.concentration = concentration;
uncertainty.totalAtoms = totalAtoms;

% Calculate uncertainty based on method
switch lower(options.method)
    case 'poisson'
        % Simple Poisson counting statistics
        % sigma(c) = sqrt(n) / N
        sigma = sqrt(counts) ./ totalAtoms;
        sigma(totalAtoms == 0, :) = 0;

        uncertainty.sigma = sigma;

        % Confidence interval using normal approximation
        z = norminv((1 + options.confidence) / 2);
        uncertainty.ciLower = max(0, concentration - z * sigma);
        uncertainty.ciUpper = min(1, concentration + z * sigma);

    case 'binomial'
        % Binomial uncertainty (recommended for APT)
        % Var(c) = c(1-c)/N
        variance = concentration .* (1 - concentration) ./ totalAtoms;
        variance(totalAtoms == 0, :) = 0;
        sigma = sqrt(variance);

        uncertainty.sigma = sigma;

        % Confidence interval using normal approximation
        z = norminv((1 + options.confidence) / 2);
        uncertainty.ciLower = max(0, concentration - z * sigma);
        uncertainty.ciUpper = min(1, concentration + z * sigma);

    case 'clopper'
        % Clopper-Pearson exact confidence interval
        % More accurate for small counts
        alpha = 1 - options.confidence;

        ciLower = zeros(nMeasurements, nElements);
        ciUpper = zeros(nMeasurements, nElements);

        for m = 1:nMeasurements
            N = totalAtoms(m);
            for e = 1:nElements
                n = counts(m, e);
                if N == 0
                    ciLower(m, e) = 0;
                    ciUpper(m, e) = 0;
                elseif n == 0
                    ciLower(m, e) = 0;
                    ciUpper(m, e) = 1 - (alpha/2)^(1/N);
                elseif n == N
                    ciLower(m, e) = (alpha/2)^(1/N);
                    ciUpper(m, e) = 1;
                else
                    ciLower(m, e) = betaincinv(alpha/2, n, N-n+1);
                    ciUpper(m, e) = betaincinv(1-alpha/2, n+1, N-n);
                end
            end
        end

        uncertainty.ciLower = ciLower;
        uncertainty.ciUpper = ciUpper;
        uncertainty.sigma = (ciUpper - ciLower) / (2 * norminv((1 + options.confidence) / 2));
end

% Additional details
details = struct();
details.method = options.method;
details.confidence = options.confidence;
details.detectionEfficiency = options.detectionEfficiency;
details.nMeasurements = nMeasurements;
details.nElements = nElements;
details.effectiveCounts = correctedCounts;

% If single measurement, simplify output
if nMeasurements == 1
    uncertainty.concentration = uncertainty.concentration(:)';
    uncertainty.sigma = uncertainty.sigma(:)';
    uncertainty.ciLower = uncertainty.ciLower(:)';
    uncertainty.ciUpper = uncertainty.ciUpper(:)';
end

end


function mustBeInRange(x, lo, hi)
    if x < lo || x > hi
        error('Value must be in range [%g, %g]', lo, hi);
    end
end
