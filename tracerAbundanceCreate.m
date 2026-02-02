function tracerAbund = tracerAbundanceCreate(element, isotopes, abundances)
% tracerAbundanceCreate creates a tracer isotope abundance table.
%
% tracerAbund = tracerAbundanceCreate(element, isotopes, abundances)
%
% This function creates a table defining the isotopic composition of a
% tracer material, which can be used by posCalculateConcentrationDeconvolved
% to perform multi-source deconvolution (separating natural vs. tracer
% contributions).
%
% INPUT
% element:    Element symbol as string or char (e.g., 'H', 'O', 'Fe')
%
% isotopes:   Vector of mass numbers for the isotopes present in the tracer
%             (e.g., [16, 18] for oxygen with 16O and 18O)
%             For a single isotope, can be a scalar (e.g., 2 for deuterium)
%
% abundances: Vector of relative abundances (fractions 0-1) for each isotope
%             Must have the same length as isotopes
%             Values are automatically normalized to sum to 1
%             For a single isotope at 100%, use 1 or 100
%
% OUTPUT
% tracerAbund: Table with columns:
%              - element: element symbol (categorical)
%              - isotope: mass number
%              - abundance: normalized abundance (0-1)
%
% EXAMPLES
%
% % Deuterium tracer (100% 2H)
% tracerAbund = tracerAbundanceCreate('H', 2, 1.0);
%
% % Oxygen tracer with 50% 16O, 50% 18O
% tracerAbund = tracerAbundanceCreate('O', [16, 18], [0.5, 0.5]);
%
% % Oxygen-18 enriched tracer (90% 18O, 10% 16O)
% tracerAbund = tracerAbundanceCreate('O', [16, 18], [0.1, 0.9]);
%
% % Carbon-13 enriched tracer (99% 13C, 1% 12C)
% tracerAbund = tracerAbundanceCreate('C', [12, 13], [0.01, 0.99]);
%
% % Multiple elements can be combined by stacking tables:
% tracerH = tracerAbundanceCreate('H', 2, 1.0);
% tracerO = tracerAbundanceCreate('O', [16, 18], [0.5, 0.5]);
% tracerAbund = [tracerH; tracerO];
%
% USAGE WITH CONCENTRATION CALCULATION
%
% % After defining ions with tracer suffix in mass spectrum:
% ionAdd(spec, 'H+');        % Natural hydrogen
% ionAdd(spec, 'H+tracer');  % Hydrogen tracer (will use deuterium abundances)
%
% % Define tracer composition
% tracerAbund = tracerAbundanceCreate('H', 2, 1.0);  % 100% deuterium
%
% % Calculate concentration with multi-source deconvolution
% [conc, info] = posCalculateConcentrationDeconvolved(pos, detEff, {}, 'ROI', ...
%     'massSpec', spec, 'tracerAbundance', tracerAbund);
%
% % Results will separate "H" (natural) from "H (tracer)" (deuterium)
%
% NOTES
% - The tracer abundance table defines ONLY the isotope ratios in the tracer
%   material, not the natural abundances
% - For single-isotope tracers (e.g., pure deuterium), the abundance is 1.0
% - For mixed-isotope tracers (e.g., enriched 18O), provide all isotopes
%   present in the tracer with their relative abundances
% - Abundances are normalized automatically, so [1, 1] is equivalent to [0.5, 0.5]
%
% SEE ALSO
% posCalculateConcentrationDeconvolved, ionAdd, ionConvertName
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Input validation
if nargin < 3
    error('tracerAbundanceCreate:insufficientInputs', ...
        'Requires element, isotopes, and abundances as inputs.');
end

% Convert element to string if needed
if ischar(element)
    element = string(element);
end

% Validate element
if ~isstring(element) && ~ischar(element)
    error('tracerAbundanceCreate:invalidElement', ...
        'Element must be a string or char (e.g., ''H'', ''O'').');
end

% Ensure isotopes is a vector
isotopes = isotopes(:);
if ~isnumeric(isotopes) || any(isotopes <= 0) || any(mod(isotopes, 1) ~= 0)
    error('tracerAbundanceCreate:invalidIsotopes', ...
        'Isotopes must be positive integers (mass numbers).');
end

% Ensure abundances is a vector with same length as isotopes
abundances = abundances(:);
if ~isnumeric(abundances) || any(abundances < 0)
    error('tracerAbundanceCreate:invalidAbundances', ...
        'Abundances must be non-negative numbers.');
end

if numel(isotopes) ~= numel(abundances)
    error('tracerAbundanceCreate:sizeMismatch', ...
        'Isotopes and abundances must have the same number of elements.');
end

%% Normalize abundances to sum to 1
totalAbund = sum(abundances);
if totalAbund <= 0
    error('tracerAbundanceCreate:zeroAbundance', ...
        'Total abundance must be greater than zero.');
end

% Handle percentage vs fraction input
if totalAbund > 1 + 1e-6  % Allow small tolerance for floating point
    % Assume percentages if total > 1
    abundances = abundances / 100;
    totalAbund = sum(abundances);
end

% Normalize
abundances = abundances / totalAbund;

%% Create output table
numIsotopes = numel(isotopes);
elementCol = repmat(categorical(element), numIsotopes, 1);

tracerAbund = table(elementCol, isotopes, abundances, ...
    'VariableNames', {'element', 'isotope', 'abundance'});

end
