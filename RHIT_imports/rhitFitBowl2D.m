function [tOffset, bowl, info] = rhitFitBowl2D(hits, epos, matchedE, matchedR, options)
% RHITFITBOWL2D Fit a 2D xy polynomial bowl from matched RHIT/EPOS events.
%
% [tOffset, bowl, info] = rhitFitBowl2D(hits, epos, matchedE, matchedR)
% [tOffset, bowl, info] = rhitFitBowl2D(..., 'bowlKind', 'xy2d', 'bowlDegree', 4)
%
% Determines:
%   1. t_offset = median(tof_rhit - tof_epos) over the matched events
%   2. C(x, y) = sum_{i+j<=deg} a_ij x^i y^j   (bowlKind = 'xy2d')
%      such that  mc_epos ≈ C(x, y) · VDC · (tof - t_offset)²
%      on the matched events.
%
% A radial fallback (bowlKind = 'radial', polynomial in r²) is also
% supported.
%
% INPUTS:
%   hits      - table from rhitLoad (detx, dety in mm; tof ns; VDC V)
%   epos      - table from posLoad  (detx, dety in mm; tof ns; VDC V; mc Da)
%   matchedE  - row indices into epos
%   matchedR  - row indices into hits, same length as matchedE
%
% NAME-VALUE:
%   'bowlKind'      - 'xy2d' (default) or 'radial'
%   'bowlDegree'    - polynomial degree (default 4)
%   'mcMin','mcMax' - mc range used for fitting (default 0.5..200 Da)
%   'tofMin'        - minimum EPOS tof for events to enter the fit (default 10 ns)
%   'cOutlierFrac'  - reject events with |C - median(C)| > frac*median(C)
%                     (default 0.2)
%
% OUTPUTS:
%   tOffset - TOF offset in ns
%   bowl    - struct with fields kind, degree, coeffs, terms (consumed by
%             rhitEvaluateBowl / rhitApplyCalibration)
%   info    - diagnostic struct (residual std, n_used, L_eff at radii etc.)
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    epos table
    matchedE (:,1) double
    matchedR (:,1) double
    options.bowlKind (1,:) char {mustBeMember(options.bowlKind, {'xy2d','radial'})} = 'xy2d'
    options.bowlDegree (1,1) double {mustBeInteger, mustBePositive} = 4
    options.mcMin (1,1) double = 0.5
    options.mcMax (1,1) double = 200
    options.tofMin (1,1) double = 10
    options.cOutlierFrac (1,1) double = 0.2
end

if numel(matchedE) ~= numel(matchedR)
    error('rhitFitBowl2D:sizeMismatch', 'matchedE and matchedR must match in length.');
end

% --- t_offset ---
tofR = double(hits.tof(matchedR));
tofE = double(epos.tof(matchedE));
tOffset = median(tofR - tofE);
tOffsetStd = std(tofR - tofE);

% --- Per-event calibration constant C_i = mc_e / (VDC * tof_e^2) ---
mcE = double(epos.mc(matchedE));
vdcM = double(epos.VDC(matchedE));
valid = mcE > options.mcMin & mcE < options.mcMax & ...
        tofE > options.tofMin & vdcM > 0;
nValid = sum(valid);

C_perEvent = mcE(valid) ./ (vdcM(valid) .* tofE(valid).^2);
detxM = double(hits.detx(matchedR(valid)));
detyM = double(hits.dety(matchedR(valid)));

C_med = median(C_perEvent);
good = abs(C_perEvent - C_med) < options.cOutlierFrac * C_med;
nGood = sum(good);

% --- Build design matrix ---
switch options.bowlKind
    case 'xy2d'
        [A, terms] = rhitDesignXY(detxM(good), detyM(good), options.bowlDegree);
    case 'radial'
        r2 = detxM(good).^2 + detyM(good).^2;
        terms = (0:options.bowlDegree)';
        A = zeros(numel(r2), options.bowlDegree + 1);
        for k = 0:options.bowlDegree
            A(:, k + 1) = r2 .^ k;
        end
end

coeffs = A \ C_perEvent(good);
resid = C_perEvent(good) - A * coeffs;
residStd = std(resid);
C0 = coeffs(1);

bowl = struct();
bowl.kind = options.bowlKind;
bowl.degree = options.bowlDegree;
bowl.coeffs = coeffs;
bowl.terms = terms;

% --- Effective flight path at a few radii (sanity) ---
e_const = 1.602176634e-19;
u_const = 1.66053906660e-27;
LeffByR = struct();
for r = [0, 5, 10, 15]
    if strcmp(options.bowlKind, 'xy2d')
        ang = linspace(0, 2*pi, 64);
        ang(end) = [];
        Cv = rhitEvaluateBowl(bowl, r * cos(ang), r * sin(ang));
        Cmean = mean(Cv);
    else
        Cmean = rhitEvaluateBowl(bowl, r, 0);
    end
    if Cmean > 0
        LeffByR.(sprintf('r%d', r)) = ...
            sqrt(2 * e_const / (u_const * Cmean) * 1e-18) * 1000;
    end
end

% --- Info struct ---
info = struct();
info.tOffsetNs = tOffset;
info.tOffsetStdNs = tOffsetStd;
info.nMatched = numel(matchedE);
info.nValidForBowl = nValid;
info.nGoodForFit = nGood;
info.fitResidualStd = residStd;
info.fitResidualPctC0 = 100 * residStd / C0;
info.LeffMmByRMm = LeffByR;
info.C0 = C0;
end
