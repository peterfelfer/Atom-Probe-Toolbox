function mc = tofToMassToCharge(tof, VDC, detx, dety, flightPathLength, t0, options)
% TOFTOMASSTOCHARGE Convert time-of-flight to mass-to-charge ratio.
%
% mc = tofToMassToCharge(tof, VDC, detx, dety, flightPathLength, t0)
% mc = tofToMassToCharge(..., 'mode', 'voltage', 'VP', VP)
%
% INPUT
%   tof              - time-of-flight (ns), column vector
%   VDC              - DC specimen voltage (V), column vector
%   detx             - detector x position (mm), column vector
%   dety             - detector y position (mm), column vector
%   flightPathLength - nominal flight path (mm), scalar
%   t0               - propagation delay (ns), scalar
%
% OPTIONS
%   'VP'    - pulse voltage (V), column vector (required for voltage mode)
%   'mode'  - 'laser' (default) or 'voltage'
%   'alpha' - amplification factor (default 1.015)
%   'beta'  - averaging factor (default 0.7)
%
% OUTPUT
%   mc - mass-to-charge ratio (Da), column vector
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    tof (:,1) double
    VDC (:,1) double
    detx (:,1) double
    dety (:,1) double
    flightPathLength (1,1) double
    t0 (1,1) double
    options.VP (:,1) double = []
    options.mode (1,:) char {mustBeMember(options.mode, {'laser','voltage'})} = 'laser'
    options.alpha (1,1) double = 1.015
    options.beta (1,1) double = 0.7
end

% Physical constants
e   = 1.602176634e-19;     % elementary charge (C)
amu = 1.66053906660e-27;   % atomic mass unit (kg)

% Unit conversions: all to SI
t_s = (tof - t0) * 1e-9;                   % ns  -> s
L   = sqrt((detx*1e-3).^2 + ...            % mm  -> m
            (dety*1e-3).^2 + ...
            (flightPathLength*1e-3).^2);

switch options.mode
    case 'laser'
        mc = 2 .* VDC .* e .* (t_s ./ L).^2;

    case 'voltage'
        if isempty(options.VP)
            error('tofToMassToCharge:missingVP', ...
                'VP (pulse voltage) is required for voltage mode.');
        end
        mc = 2 .* options.alpha .* (VDC + options.beta .* options.VP) ...
             .* e .* (t_s ./ L).^2;
end

mc = mc ./ amu;   % kg -> Da

end
