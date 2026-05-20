function hits = strCalibrateFromRhit(hits, rhitHits, rhitHistograms, rhitMetadata)
% STRCALIBRATEFROMRHIT Calibrate STR data using a matching RHIT file.
%
% hits = strCalibrateFromRhit(hits, rhitHits, rhitHistograms, rhitMetadata)
%
% Uses the matching RHIT file (same run) to transfer per-event voltage and
% determine the TDC clock period, then computes mass-to-charge ratio.
%
% Added columns:
%   VDC      - specimen DC voltage (V, transferred from RHIT)
%   tof_ns   - time-of-flight in nanoseconds
%   mc       - mass-to-charge ratio (Da, voltage-corrected)
%
% The voltage is transferred by mapping fractional run position (STR event
% index / total) to the RHIT voltage curve. This is accurate to ~0.03%
% because voltage changes slowly (~1 mV/event).
%
% The TDC clock period is determined by matching the dominant TOF peak
% position between RHIT (in ns) and STR (in TDC counts).
%
% INPUTS:
%   hits            - Table from strLoad + strCalculatePositions
%   rhitHits        - Table from rhitLoad (same run)
%   rhitHistograms  - Histograms from rhitLoad
%   rhitMetadata    - Metadata from rhitLoad
%
% EXAMPLE:
%   [strH, strM] = strLoad('data.STR');
%   strH = strCalculatePositions(strH);
%   [rhitH, rhitHist, rhitMeta] = rhitLoad('data.RHIT');
%   strH = strCalibrateFromRhit(strH, rhitH, rhitHist, rhitMeta);
%   massSpecPlot(strH, 0.01, 'normalised');
%
% See also: strLoad, strCalculatePositions, rhitLoad
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hits table
    rhitHits table
    rhitHistograms struct
    rhitMetadata struct
end

nSTR = height(hits);
nRHIT = height(rhitHits);
fprintf('Calibrating %d STR events from %d RHIT events\n', nSTR, nRHIT);

%% Transfer voltage via fractional run position
strFrac = (1:nSTR)' / nSTR;
rhitIdx = max(1, min(nRHIT, round(strFrac * nRHIT)));
hits.VDC = rhitHits.VDC(rhitIdx);

fprintf('Voltage transferred: [%.1f, %.1f] V\n', min(hits.VDC), max(hits.VDC));

%% Determine TDC clock and t0 by optimizing mc spectrum match
% mc = C * V * ((tof_TDC - t0) * clock)^2
% Optimize clock and t0 so the STR mc spectrum matches the RHIT mc spectrum.

L = rhitMetadata.instrumentParams.flight_path_mm / 1000;
kf = double(rhitHits.Vref(1)) / double(rhitHits.VDC(1));
e_const = 1.602176634e-19;
u_const = 1.66053906660e-27;
C = 2 * e_const / (u_const * L^2) * 1e-18 / sqrt(kf);

% Reference: RHIT mc spectrum
mcEdges = 0:0.05:120;
hRef = double(histcounts(rhitHits.mc, mcEdges))';

% Subsample STR for speed during optimization
stride = 10;
tofSub = hits.tof(1:stride:end);
vdcSub = hits.VDC(1:stride:end);
valid = ~isnan(tofSub);
tofSub = tofSub(valid);
vdcSub = vdcSub(valid);

% Cost function: negative correlation between STR and RHIT mc spectra
costFn = @(p) -corr( ...
    double(histcounts(C .* vdcSub .* ((tofSub - p(1)) * p(2)).^2, mcEdges))', ...
    hRef);

% Initial guess: clock ~30 ps, t0 ~ 0 (start broad)
p0 = [0, 0.030];
opts = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-14, 'MaxIter', 500);

% Multi-start to avoid local minima
bestCost = inf;
bestP = p0;
for t0_try = linspace(-5000, 30000, 20)
    for clk_try = [0.025, 0.030, 0.035, 0.040, 0.050]
        try
            [pOpt, fval] = fminsearch(costFn, [t0_try, clk_try], opts);
            if fval < bestCost && pOpt(2) > 0.01 && pOpt(2) < 0.1
                bestCost = fval;
                bestP = pOpt;
            end
        catch
        end
    end
end

t0TDC = bestP(1);
clockNs = bestP(2);

fprintf('Optimized: t0=%.1f TDC (%.1f ns), clock=%.4f ns (%.2f ps)\n', ...
    t0TDC, t0TDC * clockNs, clockNs, clockNs * 1000);
fprintf('Spectrum correlation: %.4f\n', -bestCost);

%% Convert TOF to nanoseconds and compute mc
hits.tof_ns = (hits.tof - t0TDC) * clockNs;
hits.mc = C .* double(hits.VDC) .* hits.tof_ns.^2;

%% Set units
hits.Properties.VariableUnits{strcmp(hits.Properties.VariableNames, 'VDC')} = 'V';
hits.Properties.VariableUnits{strcmp(hits.Properties.VariableNames, 'tof_ns')} = 'ns';
hits.Properties.VariableUnits{strcmp(hits.Properties.VariableNames, 'mc')} = 'Da';

%% Report peak accuracy
ms = rhitHistograms.massSpectrum;
mcC = (double(ms.edges(1:end-1)) + double(ms.edges(2:end))) / 2;
[~, pkRef] = max(double(ms.values));
mcRef = mcC(pkRef);

mcEdges = mcRef-3:0.01:mcRef+3;
mcCC = mcEdges(1:end-1) + 0.005;
hStr = histcounts(hits.mc, mcEdges);
[~, pkCalc] = max(hStr);

fprintf('Main peak: pMass=%.3f Da, STR=%.3f Da (%+.3f Da)\n', ...
    mcRef, mcCC(pkCalc), mcCC(pkCalc) - mcRef);
fprintf('Calibration complete: %d events with mc\n', sum(~isnan(hits.mc)));

end
