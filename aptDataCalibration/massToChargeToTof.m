function tof = massToChargeToTof(mc, VDC, detx, dety, flightPathLength)
% MASSTOCHARGETOTOF Convert mass-to-charge ratio to time-of-flight.
%
% tof = massToChargeToTof(mc, VDC, detx, dety, flightPathLength)
%
% INPUT
%   mc               - mass-to-charge ratio (Da), scalar or column vector
%   VDC              - DC specimen voltage (V), scalar or column vector
%   detx             - detector x position (mm), scalar or column vector
%   dety             - detector y position (mm), scalar or column vector
%   flightPathLength - nominal flight path (mm), scalar
%
% OUTPUT
%   tof - time-of-flight (ns), same size as mc
%
% NOTE: This is the simple inverse (laser mode, no t0 offset).
%       Add t0 to the result if needed.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    mc (:,1) double
    VDC (:,1) double
    detx (:,1) double
    dety (:,1) double
    flightPathLength (1,1) double
end

e   = 1.602176634e-19;
amu = 1.66053906660e-27;

L = sqrt((detx*1e-3).^2 + (dety*1e-3).^2 + (flightPathLength*1e-3).^2);

tof = sqrt((mc .* amu .* L.^2) ./ (2 .* e .* VDC)) * 1e9;

end
