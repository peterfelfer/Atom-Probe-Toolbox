function [hits, calib, diagnostics] = rhitCalibrateFromStoredMass(hits, histograms, metadata, options)
% RHITCALIBRATEFROMSTOREDMASS RHIT-only mass calibration from stored pMass.
%
% [hits, calib, diagnostics] = rhitCalibrateFromStoredMass(hits, histograms, metadata)
%
% Uses the calibrated mass histogram stored inside the RHIT file (pMass) as
% the spectral target. A detector-position dependent C(x,y) is estimated by
% aligning local detector-cell spectra to pMass, then represented as a 2D
% polynomial bowl usable by rhitApplyCalibration.

arguments
    hits table
    histograms struct
    metadata struct
    options.mcHi (1,1) double = 100
    options.binWidth (1,1) double = 0.02
    options.gridHalf (1,1) double = 16
    options.gridN (1,1) double {mustBeInteger, mustBeGreaterThan(options.gridN, 2)} = 9
    options.bowlDegree (1,1) double {mustBeInteger, mustBePositive} = 4
    options.minCellHits (1,1) double = 6000
    options.scaleHalfRange (1,1) double = 0.12
    options.scaleStepCoarse (1,1) double = 0.01
    options.scaleStepFine (1,1) double = 0.0015
    options.refineTOffset (1,1) logical = true
    options.verbose (1,1) logical = true
end

if ~isfield(histograms, 'massSpectrum')
    error('rhitCalibrateFromStoredMass:noMassSpectrum', ...
        'histograms.massSpectrum is required.');
end

edges = 0:options.binWidth:options.mcHi;
refValues = rebinStoredSpectrum(histograms.massSpectrum, edges);
refValues = double(refValues(:));
refLog = log1p(refValues);
refNorm = refLog ./ max(sum(refLog), eps);
mask = edges(1:end-1)' >= 0.5 & edges(1:end-1)' <= options.mcHi;

tOffset = initialTOffset(metadata);
C0 = initialC0(hits, tOffset);

if options.verbose
    fprintf('RHIT-only pMass calibration: t0 seed %.4f ns, C0 %.6e\n', ...
        tOffset, C0);
end

[bowl, cellTable] = fitBowlFromCells(hits, C0, tOffset, edges, refNorm, mask, options);

if options.refineTOffset
    refSpectrum = struct('edges', edges, 'values', refValues');
    [tOffset, tInfo] = rhitRefineTOffset(hits, refSpectrum, bowl, tOffset, ...
        'mcHi', options.mcHi, ...
        'coarseHalf', 20, ...
        'coarseStep', 1, ...
        'fineHalf', 2, ...
        'fineStep', 0.1);
    if options.verbose
        fprintf('Refined t_offset against stored pMass: %.4f ns\n', tOffset);
    end
    [bowl, cellTable] = fitBowlFromCells(hits, C0, tOffset, edges, refNorm, mask, options);
else
    tInfo = struct();
end

calib = struct();
calib.tOffsetNs = tOffset;
calib.bowl = bowl;
calib.reference = 'RHIT histograms.massSpectrum (pMass)';
calib.method = 'local detector-cell spectrum alignment to stored pMass';
calib.nCells = height(cellTable);
calib.fitResidualPctC0 = bowl.fitResidualPctC0;

rawSimilarity = spectrumSimilarity(hits.mc, edges, refNorm, mask);
hits = rhitApplyCalibration(hits, calib);
calibratedSimilarity = spectrumSimilarity(hits.mc, edges, refNorm, mask);

diagnostics = struct();
diagnostics.cellTable = cellTable;
diagnostics.tOffsetInfo = tInfo;
diagnostics.referenceEdges = edges;
diagnostics.referenceValues = refValues;
diagnostics.rawSimilarity = rawSimilarity;
diagnostics.calibratedSimilarity = calibratedSimilarity;
diagnostics.calibration = calib;
end


function t0 = initialTOffset(metadata)
if isfield(metadata, 'instrumentParams') && isfield(metadata.instrumentParams, 't0_ns')
    t0 = double(metadata.instrumentParams.t0_ns);
else
    t0 = 45;
end
end


function C0 = initialC0(hits, tOffset)
base = double(hits.VDC) .* (double(hits.tof) - tOffset).^2;
ok = isfinite(base) & base > 0 & isfinite(hits.mc) & hits.mc > 0;
C0 = median(double(hits.mc(ok)) ./ base(ok));
end


function [bowl, cellTable] = fitBowlFromCells(hits, C0, tOffset, edges, refNorm, mask, opt)
xEdges = linspace(-opt.gridHalf, opt.gridHalf, opt.gridN + 1);
yEdges = linspace(-opt.gridHalf, opt.gridHalf, opt.gridN + 1);

cellX = [];
cellY = [];
cellC = [];
cellN = [];
cellScore = [];

baseAll = double(hits.VDC) .* (double(hits.tof) - tOffset).^2;
xAll = double(hits.detx);
yAll = double(hits.dety);
validAll = isfinite(baseAll) & baseAll > 0 & isfinite(xAll) & isfinite(yAll);

for ix = 1:opt.gridN
    for iy = 1:opt.gridN
        inCell = validAll & xAll >= xEdges(ix) & xAll < xEdges(ix+1) & ...
            yAll >= yEdges(iy) & yAll < yEdges(iy+1);
        nCell = sum(inCell);
        if nCell < opt.minCellHits
            continue
        end

        base = baseAll(inCell);
        [scale, score] = bestScaleForCell(base, C0, edges, refNorm, mask, opt);
        xc = 0.5 * (xEdges(ix) + xEdges(ix+1));
        yc = 0.5 * (yEdges(iy) + yEdges(iy+1));

        cellX(end+1, 1) = xc; %#ok<AGROW>
        cellY(end+1, 1) = yc; %#ok<AGROW>
        cellC(end+1, 1) = C0 * scale; %#ok<AGROW>
        cellN(end+1, 1) = nCell; %#ok<AGROW>
        cellScore(end+1, 1) = score; %#ok<AGROW>
    end
end

if numel(cellC) < nchoosek(opt.bowlDegree + 2, 2)
    error('rhitCalibrateFromStoredMass:tooFewCells', ...
        'Only %d usable detector cells for degree-%d fit.', numel(cellC), opt.bowlDegree);
end

[A, terms] = rhitDesignXY(cellX, cellY, opt.bowlDegree);
w = sqrt(cellN ./ median(cellN));
coeffs = (A .* w) \ (cellC .* w);
resid = cellC - A * coeffs;
C0Fit = coeffs(1);

bowl = struct();
bowl.kind = 'xy2d';
bowl.degree = opt.bowlDegree;
bowl.coeffs = coeffs;
bowl.terms = terms;
bowl.C0Global = C0Fit;
bowl.fitResidualPctC0 = 100 * std(resid) / C0Fit;

cellTable = table(cellX, cellY, cellC, cellN, cellScore, resid, ...
    'VariableNames', {'detx','dety','C','nHits','score','residual'});

if opt.verbose
    fprintf('Fitted pMass bowl from %d cells: C0 %.6e, residual %.4f%%\n', ...
        height(cellTable), C0Fit, bowl.fitResidualPctC0);
end
end


function [bestScale, bestScore] = bestScaleForCell(base, C0, edges, refNorm, mask, opt)
coarse = (1 - opt.scaleHalfRange):opt.scaleStepCoarse:(1 + opt.scaleHalfRange);
[bestScale, bestScore] = scanScale(base, C0, edges, refNorm, mask, coarse);
fine = (bestScale - 3 * opt.scaleStepCoarse):opt.scaleStepFine:(bestScale + 3 * opt.scaleStepCoarse);
[bestScale, bestScore] = scanScale(base, C0, edges, refNorm, mask, fine);
end


function [bestScale, bestScore] = scanScale(base, C0, edges, refNorm, mask, scales)
bestScale = scales(1);
bestScore = -Inf;
for k = 1:numel(scales)
    mc = C0 * scales(k) .* base;
    score = spectrumSimilarity(mc, edges, refNorm, mask);
    if score > bestScore
        bestScore = score;
        bestScale = scales(k);
    end
end
end


function score = spectrumSimilarity(mc, edges, refNorm, mask)
h = histcounts(mc(mc > 0 & mc < edges(end)), edges);
hLog = log1p(double(h(:)));
hNorm = hLog ./ max(sum(hLog), eps);
score = sum(hNorm(mask) .* refNorm(mask));
end


function values = rebinStoredSpectrum(ms, edges)
oldEdges = double(ms.edges(:)');
oldValues = double(ms.values(:)');
oldCenters = 0.5 * (oldEdges(1:end-1) + oldEdges(2:end));
targetBin = discretize(oldCenters, edges);
values = accumarray(targetBin(~isnan(targetBin))', oldValues(~isnan(targetBin))', ...
    [numel(edges)-1, 1], @sum, 0);
end
