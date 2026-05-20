function hits = rhitApplyCalibration(hits, calib)
% RHITAPPLYCALIBRATION Apply a saved calibration to RHIT data.
%
% hits = rhitApplyCalibration(hits, calib)
%
% Computes  mc = C(x, y) · VDC · (tof - t_offset)²  and writes it into
% hits.mc.  Accepts both new (2D bowl, struct field `bowl`) and legacy
% (radial bowl, field `C_poly`) calibration structs.
%
% INPUTS:
%   hits  - Table from rhitLoad
%   calib - Calibration struct.  Either:
%             new format: fields tOffsetNs + bowl (struct with kind,
%                         degree, coeffs, terms)
%             legacy:     fields t_offset (ns) + C_poly = [C0, C1, C2]
%                         for C(r²) = C0 + C1·r² + C2·r⁴
%
% See also: rhitCalibrateFromEpos, rhitLoad, rhitEvaluateBowl, rhitPoolBowl
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    calib struct
end

if isfield(calib, 'bowl') && isstruct(calib.bowl)
    bowl = calib.bowl;
    tOffset = calib.tOffsetNs;
elseif isfield(calib, 'C_poly')
    % Legacy radial polynomial
    bowl = struct( ...
        'kind', 'radial', ...
        'degree', max(1, numel(calib.C_poly) - 1), ...
        'coeffs', calib.C_poly(:), ...
        'terms', (0:numel(calib.C_poly)-1)');
    if isfield(calib, 't_offset')
        tOffset = calib.t_offset;
    else
        tOffset = calib.tOffsetNs;
    end
else
    error('rhitApplyCalibration:badStruct', ...
        'Calibration struct must have either .bowl or .C_poly field.');
end

C = rhitEvaluateBowl(bowl, double(hits.detx), double(hits.dety));
tofCorr = double(hits.tof) - tOffset;
hits.mc = C .* double(hits.VDC) .* tofCorr .^ 2;

if isfield(calib, 'fitResidualPctC0')
    fprintf('Applied calibration: t_offset=%.2f ns, bowl=%s deg=%d (resid %.3f%%)\n', ...
        tOffset, bowl.kind, bowl.degree, calib.fitResidualPctC0);
elseif isfield(calib, 'fit_residual_std')
    fprintf('Applied calibration: t_offset=%.2f ns (legacy radial)\n', tOffset);
else
    fprintf('Applied calibration: t_offset=%.2f ns\n', tOffset);
end
end
