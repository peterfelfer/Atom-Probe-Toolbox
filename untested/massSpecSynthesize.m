function [mc, gt, info] = massSpecSynthesize(species, options)
% MASSSPECSYNTHESIZE Generate a synthetic atom probe mass spectrum.
%
% Builds a synthetic mass-to-charge spectrum with known ground truth,
% following the methodology of Hudson, Smith and Gault (2011): generate
% in time-of-flight (TOF) space, apply an exponentially modified Gaussian
% peak shape, add a constant TOF background that becomes A/sqrt(mc) under
% the TOF -> m/q transformation, and sample one entry per ion. Returns the
% list of synthesized m/q values, a ground-truth peak table, and an info
% struct with the parameters used.
%
% [mc, gt, info] = massSpecSynthesize(species)
% [mc, gt, info] = massSpecSynthesize(species, 'tofSigma_ns', 0.5, ...)
%
% INPUT
% species   table, one row per ion species with columns:
%             element       string/char element symbol, e.g. "Fe"
%             chargeState   positive integer, e.g. 2 for Fe2+
%             counts        nonnegative integer, total ions for the species
%             isotopes      (optional) numeric vector of mass numbers to
%                           restrict synthesis to (e.g. [56 54]); abundances
%                           are renormalized over the kept isotopes. When
%                           omitted or empty, all natural isotopes are used.
%
% OPTIONS
%   isotopeTable     isotope table; default loads
%                    isotopeTable_naturalAbundances.mat
%   t0_ns            TOF zero offset (default: 50 ns)
%   tofScale         TOF scale factor s.t. tof = t0 + tofScale*sqrt(mc),
%                    in ns/sqrt(Da) (default: 1000)
%   tofSigma_ns      Gaussian sigma of peak in TOF (default: 0.6 ns)
%                    Approximate mass resolution m/dm at FWHM is
%                       m/dm = tof_peak / (2 * 2.355 * sigma).
%   tofTau_ns        exponential decay time of EMG tail in TOF
%                    (default: 1.0 ns; set 0 for pure Gaussian).
%                    For laser-pulsed-mode realism, push to 3-8 ns.
%   tailPowerExp     power-law exponent p for an additional heavy
%                    right-side kick on top of the EMG core (default: NaN,
%                    disabled). When set to a finite value > 1, a fraction
%                    tailPowerFrac of each species' ions receive an extra
%                    Lomax-distributed offset on top of their EMG sample,
%                    smoothly continuous with the peak (kick = 0 at U = 1),
%                    decaying as dt^(-p) in TOF for large dt. Smaller p =
%                    heavier tail; p = 2 has divergent mean and gives
%                    laser-mode-like very long tails.
%   tailPowerFrac    fraction of each species' ions that receive the heavy
%                    kick (default: 0). For real laser-pulsed data
%                    0.10-0.30 is realistic. The other (1-tailPowerFrac)
%                    of ions still go through the EMG core, so the peak
%                    stays single-modal.
%   tailStart_ns     scale parameter of the Lomax kick (default: NaN,
%                    falls back to tofSigma_ns). Larger values stretch the
%                    heavy tail out further. Voltage-pulsed data: leave at
%                    default. Laser-pulsed data: 2-10 ns.
%   bgIons           total number of background ions, distributed uniformly
%                    in TOF over the synthesis range (default: 0)
%   mcRange          [mcMin mcMax] for synthesis and background, in Da
%                    (default: [0 100])
%   randomSeed       rng seed for reproducibility (default: [], no reset)
%
% OUTPUT
% mc        column vector of synthesized mass-to-charge values [Da].
% gt        ground-truth table, one row per (isotope, chargeState):
%             element, isotope, chargeState, peakMc, peakTof_ns,
%             isotopeAbundance, trueCounts, ionLabel
% info      struct with all parameters, the TOF<->m/q mappings, expected
%           and realized background counts, total ion count, and a
%           backgroundCoefA field (the analytic A in A/sqrt(mc) per Da).
%
% TOF/m-q mapping
%   tof = t0_ns + tofScale * sqrt(mc),           mc = ((tof - t0_ns)/tofScale)^2
% A constant TOF density rho_t (ions/ns) maps to a m/q density
%   rho_mc(mc) = rho_t * tofScale / (2*sqrt(mc)) = A/sqrt(mc),  A = rho_t*tofScale/2.
%
% Generation algorithm (after Hudson et al 2011, Sect. 2.1)
%   1. Resolve per-isotope true peak positions:  peakMc = isoMass / cs.
%   2. Distribute species counts across isotopes by natural abundance.
%      Sub-counts are deterministic round() so ground truth is exact.
%   3. For each ion, draw observed TOF as
%        tof_obs = peakTof + sigma*randn() - tau*log(rand())
%      i.e. exponentially modified Gaussian (EMG) with right-side tail.
%   4. Convert tof_obs to mc via the inverse mapping.
%   5. Background: draw bgIons uniform samples in [tofMin, tofMax],
%      convert to mc.  This produces the experimentally observed
%      A/sqrt(mc) background shape in the m/q domain.
%
% LIMITATIONS
%   - Monatomic ions only. Molecular ions and isotope-mixed clusters
%     are not produced. Add by extending species rows or by combining
%     outputs from multiple calls.
%   - No multi-hit / detector-deadtime / pulse-by-pulse drift effects.
%   - Single global PSF (sigma, tau) shared across all peaks.
%
% EXAMPLE
%   species = table( ...
%       ["Fe"; "Cr"; "C"], ...
%       [2; 2; 1], ...
%       [1e5; 5e3; 1e3], ...
%       'VariableNames', {'element','chargeState','counts'});
%   [mc, gt] = massSpecSynthesize(species, ...
%       'tofSigma_ns', 0.6, 'tofTau_ns', 1.0, ...
%       'bgIons', 5e3, 'mcRange', [0 60], 'randomSeed', 1);
%   disp(gt);
%   massSpecPlot(mc, 0.01, 'count');
%
% References
%   Hudson, Smith, Gault. Ultramicroscopy 111 (2011) 480-486.
%       doi:10.1016/j.ultramic.2010.11.007
%   Johnson et al. Ultramicroscopy 132 (2013) 60-67.
%       doi:10.1016/j.ultramic.2013.03.015
%   Coakley, Sanford. Ultramicroscopy 240 (2022) 113521.
%       doi:10.1016/j.ultramic.2022.113521
%
% See also: rangeAutoEER, massSpecPlot, massSpecFindPeaks
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    species table
    options.isotopeTable table = table()
    options.t0_ns (1,1) double = 50
    options.tofScale (1,1) double {mustBePositive} = 1000
    options.tofSigma_ns (1,1) double {mustBeNonnegative} = 0.6
    options.tofTau_ns (1,1) double {mustBeNonnegative} = 1.0
    options.tailPowerExp (1,1) double = NaN
    options.tailPowerFrac (1,1) double {mustBeInRange(options.tailPowerFrac, 0, 1)} = 0
    options.tailStart_ns (1,1) double = NaN
    options.bgIons (1,1) double {mustBeNonnegative} = 0
    options.mcRange (1,2) double {mustBeNonnegative} = [0 100]
    options.randomSeed = []
end

if ~isempty(options.randomSeed)
    rng(options.randomSeed);
end

%% Resolve isotope table and detect column names
if isempty(options.isotopeTable) || height(options.isotopeTable) == 0
    s = load('isotopeTable_naturalAbundances.mat');
    isoT = s.isotopeTable;
else
    isoT = options.isotopeTable;
end

names = isoT.Properties.VariableNames;
elemColName  = pickColName(names, {'element','symbol','Element','Symbol'});
isoColName   = pickColName(names, {'isotope','massnumber','mass_number','A','MassNumber'});
massColName  = pickColName(names, {'weight','mass','atomicmass','Mass','AtomicMass'});
abundColName = pickColName(names, {'abundance','naturalabundance','relativeabundance','Abundance'});

elemCol  = string(isoT.(elemColName));
isoCol   = double(isoT.(isoColName));
massCol  = double(isoT.(massColName));
abundCol = double(isoT.(abundColName));
if any(abundCol > 1)
    abundCol = abundCol / 100;
end

%% TOF mapping
t0 = options.t0_ns;
scale = options.tofScale;
tofOfMc = @(mcv) t0 + scale .* sqrt(max(mcv, 0));
mcOfTof = @(tofv) ((max(tofv, t0) - t0) ./ scale).^2;

mcMin = max(options.mcRange(1), 0);
mcMax = options.mcRange(2);
if mcMax <= mcMin
    error('massSpecSynthesize:badRange', 'mcRange must be increasing.');
end
tofMin = tofOfMc(mcMin);
tofMax = tofOfMc(mcMax);

%% Validate species table
required = ["element","chargeState","counts"];
for r = 1:numel(required)
    if ~ismember(required(r), species.Properties.VariableNames)
        error('massSpecSynthesize:missingColumn', ...
            'species table must contain column "%s".', required(r));
    end
end
hasIsoFilter = ismember("isotopes", species.Properties.VariableNames);

%% Synthesize signal ions species by species
gtRows = cell(0, 8);
mcChunks = cell(0, 1);

for s = 1:height(species)
    elem = string(species.element(s));
    cs   = double(species.chargeState(s));
    nIon = round(double(species.counts(s)));
    if cs <= 0
        error('massSpecSynthesize:badChargeState', ...
            'Row %d: chargeState must be a positive integer.', s);
    end
    if nIon <= 0
        continue;
    end

    elemMask = strcmpi(elemCol, elem);
    if ~any(elemMask)
        warning('massSpecSynthesize:elementNotFound', ...
            'Element "%s" not in isotope table; skipping row %d.', elem, s);
        continue;
    end

    isoNumbers = isoCol(elemMask);
    isoMasses  = massCol(elemMask);
    isoAbund   = abundCol(elemMask);

    if hasIsoFilter
        wanted = species.isotopes{s};
        if ~isempty(wanted)
            keep = ismember(isoNumbers, wanted);
            if ~any(keep)
                warning('massSpecSynthesize:isotopeNotFound', ...
                    'No requested isotopes of %s found; skipping row %d.', elem, s);
                continue;
            end
            isoNumbers = isoNumbers(keep);
            isoMasses  = isoMasses(keep);
            isoAbund   = isoAbund(keep);
        end
    end

    sumA = sum(isoAbund);
    if sumA <= 0
        warning('massSpecSynthesize:zeroAbundance', ...
            'Row %d: isotope abundances sum to zero; skipping.', s);
        continue;
    end
    isoAbund = isoAbund / sumA;

    perIso = round(nIon * isoAbund);

    for k = 1:numel(isoNumbers)
        n = perIso(k);
        if n <= 0
            continue;
        end
        peakMc  = isoMasses(k) / cs;
        peakTof = tofOfMc(peakMc);

        % EMG core for ALL ions (single peak position, no second bump).
        gaussPart = options.tofSigma_ns * randn(n, 1);
        if options.tofTau_ns > 0
            expPart = -options.tofTau_ns * log(rand(n, 1));
        else
            expPart = zeros(n, 1);
        end

        % Optional heavy power-law kick added to a fraction of ions.
        % Lomax (shifted Pareto): dt = tailStart * (U^(-1/(p-1)) - 1).
        % At U=1 the kick is 0 (continuous with the EMG core), at U=0 the
        % kick diverges -> heavy right tail. Smaller p = heavier tail.
        usePowerTail = isfinite(options.tailPowerExp) && options.tailPowerExp > 1 ...
                       && options.tailPowerFrac > 0;
        heavyKick = zeros(n, 1);
        if usePowerTail
            p = options.tailPowerExp;
            tailStart = options.tailStart_ns;
            if ~isfinite(tailStart) || tailStart <= 0
                tailStart = options.tofSigma_ns;
            end
            heavyMask = rand(n, 1) < options.tailPowerFrac;
            nHeavy = sum(heavyMask);
            if nHeavy > 0
                u = rand(nHeavy, 1);
                heavyKick(heavyMask) = tailStart * (u.^(-1/(p - 1)) - 1);
            end
        end

        tofIon = peakTof + gaussPart + expPart + heavyKick;

        % Reject ions outside the synthesis range
        keep = tofIon > t0 & tofIon >= tofMin & tofIon <= tofMax;
        tofIon = tofIon(keep);
        mcIon = mcOfTof(tofIon);

        mcChunks{end+1, 1} = mcIon; %#ok<AGROW>

        ionLabel = sprintf('%d%s%s', isoNumbers(k), char(elem), repmat('+', 1, cs));
        gtRows(end+1, :) = {char(elem), isoNumbers(k), cs, peakMc, peakTof, ...
            isoAbund(k), n, ionLabel}; %#ok<AGROW>
    end
end

%% Background: uniform in TOF over [tofMin, tofMax]
nBg = round(options.bgIons);
if nBg > 0
    tofBg = tofMin + (tofMax - tofMin) * rand(nBg, 1);
    mcBg  = mcOfTof(tofBg);
    mcChunks{end+1, 1} = mcBg;
end

if isempty(mcChunks)
    mc = zeros(0, 1);
else
    mc = vertcat(mcChunks{:});
end
mc = mc(mc >= mcMin & mc <= mcMax);

%% Ground-truth table
if isempty(gtRows)
    gt = table();
else
    gt = cell2table(gtRows, ...
        'VariableNames', {'element','isotope','chargeState','peakMc', ...
                          'peakTof_ns','isotopeAbundance','trueCounts', ...
                          'ionLabel'});
    gt = sortrows(gt, 'peakMc');
end

%% Info struct
nSignal = 0;
if ~isempty(gt)
    nSignal = sum(gt.trueCounts);
end

% Analytic A in A/sqrt(mc): bg density per ns is nBg/(tofMax-tofMin),
% pulled through the Jacobian gives A = rho_t * tofScale / 2.
if nBg > 0 && tofMax > tofMin
    bgPerNs = nBg / (tofMax - tofMin);
else
    bgPerNs = 0;
end
backgroundCoefA = bgPerNs * scale / 2;

info = struct();
info.t0_ns = t0;
info.tofScale = scale;
info.tofSigma_ns = options.tofSigma_ns;
info.tofTau_ns = options.tofTau_ns;
info.tofRange_ns = [tofMin tofMax];
info.mcRange = [mcMin mcMax];
info.bgIons = nBg;
info.bgPerNs = bgPerNs;
info.backgroundCoefA = backgroundCoefA;
info.nIonSignal = nSignal;
info.nIonBackground = nBg;
info.nIonTotal = numel(mc);
info.tofOfMc = tofOfMc;
info.mcOfTof = mcOfTof;
info.options = options;

end


function col = pickColName(names, candidates)
    col = '';
    for i = 1:numel(candidates)
        idx = find(strcmpi(names, candidates{i}), 1, 'first');
        if ~isempty(idx)
            col = names{idx};
            return;
        end
    end
    error('massSpecSynthesize:missingIsotopeColumn', ...
        'Isotope table missing one of: %s', strjoin(candidates, ', '));
end
