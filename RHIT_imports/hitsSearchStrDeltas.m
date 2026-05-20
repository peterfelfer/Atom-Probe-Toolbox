function results = hitsSearchStrDeltas(hitsFile, strFile, options)
% HITSSEARCHSTRDELTAS Test HITS payload bytes against STR previous-ion deltas.
%
% results = hitsSearchStrDeltas(hitsFile, strFile)
%
% This is a targeted diagnostic for predictive/residual HITS encoding. STR
% provides decoded delay-line timings; the target variables are
% value(i)-value(i-1) along the loaded STR ion stream. HITS contributes raw
% byte/word predictors from the 8-byte payload.
%
% OPTIONS:
%   eventOffset  - HITS eventIdx -> STR eventIdx offset (default: 0)
%   sampleStride - stride through candidate HITS rows (default: 500)
%   maxRows      - maximum aligned rows after filtering (default: 50000)
%   parentType   - "SINGLE", "STATUS", "MULTI", "OTHER", or "ALL"
%                  (default: "SINGLE")
%   verbose      - print top correlations (default: true)
%
% See also: hitsSearchStrBitfields, hitsAuditStrAlignment
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.eventOffset (1,1) double = 0
    options.sampleStride (1,1) double = 500
    options.maxRows (1,1) double = 50000
    options.parentType (1,1) string {mustBeMember(options.parentType, ...
        ["SINGLE","STATUS","MULTI","OTHER","ALL"])} = "SINGLE"
    options.verbose (1,1) logical = true
end

[hitsSample, hitSummary] = sampleHitsFromEvents(hitsFile, options);
[strHits, ~] = strLoad(strFile);
strHits = strCalculatePositions(strHits);

maxEvent = max([double(hitsSample.eventIdx) + abs(options.eventOffset); ...
                double(strHits.eventIdx)]);
strRowByEvent = zeros(maxEvent + 2, 1, 'uint32');
strRowByEvent(double(strHits.eventIdx)) = uint32((1:height(strHits))');

eventQuery = double(hitsSample.eventIdx) + options.eventOffset;
inBounds = eventQuery >= 1 & eventQuery <= numel(strRowByEvent);
rows = zeros(height(hitsSample), 1, 'uint32');
rows(inBounds) = strRowByEvent(eventQuery(inBounds));
valid = rows > 1;

hitsSample = hitsSample(valid, :);
strRows = double(rows(valid));
strSample = strHits(strRows, :);
strPrev = strHits(strRows - 1, :);

predictors = buildPredictors(hitsSample);
targets = buildDeltaTargets(strSample, strPrev);

predictorOut = strings(0, 1);
targetOut = strings(0, 1);
rhoOut = zeros(0, 1);
rhoDetrendedOut = zeros(0, 1);
nOut = zeros(0, 1);

for pi = 1:numel(predictors.names)
    x = predictors.values{pi};
    for ti = 1:numel(targets.names)
        y = targets.values{ti};
        ok = isfinite(x) & isfinite(y);
        if nnz(ok) < 100
            continue;
        end

        rho = corr(x(ok), y(ok), 'Rows', 'complete');
        xd = detrendLocal(x(ok));
        yd = detrendLocal(y(ok));
        rhoDetrended = corr(xd, yd, 'Rows', 'complete');

        predictorOut(end+1, 1) = predictors.names(pi); %#ok<AGROW>
        targetOut(end+1, 1) = targets.names(ti); %#ok<AGROW>
        rhoOut(end+1, 1) = rho; %#ok<AGROW>
        rhoDetrendedOut(end+1, 1) = rhoDetrended; %#ok<AGROW>
        nOut(end+1, 1) = nnz(ok); %#ok<AGROW>
    end
end

results = table(predictorOut, targetOut, rhoOut, abs(rhoOut), ...
    rhoDetrendedOut, abs(rhoDetrendedOut), nOut, ...
    'VariableNames', {'predictor','target','correlation', ...
                      'absCorrelation','detrendedCorrelation', ...
                      'absDetrendedCorrelation','n'});

sortScore = results.absDetrendedCorrelation;
sortScore(isnan(sortScore)) = -Inf;
[~, ord] = sort(sortScore, 'descend');
results = results(ord, :);

if options.verbose
    fprintf('\n--- HITS/STR previous-ion delta search ---\n');
    fprintf('Parent type: %s\n', options.parentType);
    fprintf('Event offset: %+d\n', options.eventOffset);
    fprintf('Candidate payload rows before sampling: %d\n', hitSummary.nCandidateRows);
    fprintf('Sampled aligned rows: %d\n', height(hitsSample));
    disp(results(1:min(30, height(results)), :));
end
end

function [hitsSample, summary] = sampleHitsFromEvents(hitsFile, options)
[events, ~, raw] = hitsScanEvents(hitsFile);

eventFields = double(events.nFields);
completeBlocks = floor(max(0, eventFields - 1) / 2);
payloadCounts = zeros(size(eventFields));

payloadCounts((events.type == 'SINGLE' | events.type == 'STATUS' | ...
    events.type == 'OTHER') & completeBlocks > 0) = 1;
isMultiPayload = events.type == 'MULTI' & completeBlocks > 0;
payloadCounts(isMultiPayload) = completeBlocks(isMultiPayload);

typeSel = payloadCounts > 0;
if options.parentType ~= "ALL"
    typeSel = typeSel & events.type == options.parentType;
end

eventIdxAll = find(typeSel);
nCandidateRows = sum(payloadCounts(eventIdxAll));
targetRows = 1:max(1, round(options.sampleStride)):nCandidateRows;
if numel(targetRows) > options.maxRows
    targetRows = targetRows(1:options.maxRows);
end

cumRows = cumsum(payloadCounts(eventIdxAll));
eventPos = discretize(targetRows, [0; cumRows(:)]);
sampleEventIdx = eventIdxAll(eventPos);
prevCum = [0; cumRows(1:end-1)];
hitInEvent = targetRows(:) - prevCum(eventPos);

n = numel(sampleEventIdx);
B = zeros(n, 8, 'uint8');
isDelta = false(n, 1);
parentType = strings(n, 1);

for k = 1:n
    ei = sampleEventIdx(k);
    payloadOffset = double(events.byteOffset(ei)) + 4 + ...
        (double(hitInEvent(k)) - 1) * 8;
    B(k, :) = raw.bytes(payloadOffset+1:payloadOffset+8);
    isDelta(k) = bitand(B(k, 4), uint8(1)) ~= 0;
    parentType(k) = string(events.type(ei));
end

keep = ~isDelta;
B = B(keep, :);
sampleEventIdx = sampleEventIdx(keep);
hitInEvent = hitInEvent(keep);
parentType = parentType(keep);

hitsSample = table(uint32(sampleEventIdx(:)), categorical(parentType(:)), ...
    uint16(hitInEvent(:)), B(:,1), B(:,2), B(:,3), B(:,4), ...
    B(:,5), B(:,6), B(:,7), B(:,8), ...
    'VariableNames', {'eventIdx','parentType','hitInEvent', ...
                      'b0','b1','b2','b3','b4','b5','b6','b7'});

summary = struct();
summary.nEvents = numel(events.fieldIdx);
summary.nCandidateRows = nCandidateRows;
summary.nSampledRows = height(hitsSample);
summary.nDeltaSkipped = nnz(~keep);
end

function predictors = buildPredictors(hits)
byteNames = "b" + string(0:7);
names = strings(0, 1);
values = {};

for k = 1:numel(byteNames)
    x = double(hits.(char(byteNames(k))));
    names(end+1, 1) = byteNames(k); %#ok<AGROW>
    values{end+1, 1} = x; %#ok<AGROW>
    names(end+1, 1) = byteNames(k) + "_s8"; %#ok<AGROW>
    values{end+1, 1} = signedUnsigned(x, 8); %#ok<AGROW>
end

for k = 0:6
    lo = double(hits.(sprintf('b%d', k)));
    hi = double(hits.(sprintf('b%d', k+1)));
    u16 = lo + 256 * hi;
    names(end+1, 1) = "u16_b" + k + "b" + (k+1); %#ok<AGROW>
    values{end+1, 1} = u16; %#ok<AGROW>
    names(end+1, 1) = "s16_b" + k + "b" + (k+1); %#ok<AGROW>
    values{end+1, 1} = signedUnsigned(u16, 16); %#ok<AGROW>
end

for k = 0:5
    b0 = double(hits.(sprintf('b%d', k)));
    b1 = double(hits.(sprintf('b%d', k+1)));
    b2 = double(hits.(sprintf('b%d', k+2)));
    u24 = b0 + 256 * b1 + 65536 * b2;
    names(end+1, 1) = "u24_b" + k + "b" + (k+2); %#ok<AGROW>
    values{end+1, 1} = u24; %#ok<AGROW>
    names(end+1, 1) = "s24_b" + k + "b" + (k+2); %#ok<AGROW>
    values{end+1, 1} = signedUnsigned(u24, 24); %#ok<AGROW>
end

predictors = struct('names', names, 'values', {values});
end

function targets = buildDeltaTargets(strNow, strPrev)
sumXNow = double(strNow.detxt1) + double(strNow.detxt2);
sumYNow = double(strNow.detyt1) + double(strNow.detyt2);
sumWNow = double(strNow.detwt1) + double(strNow.detwt2);
sumXPrev = double(strPrev.detxt1) + double(strPrev.detxt2);
sumYPrev = double(strPrev.detyt1) + double(strPrev.detyt2);
sumWPrev = double(strPrev.detwt1) + double(strPrev.detwt2);

names = ["d_detxt1","d_detxt2","d_detyt1","d_detyt2","d_detwt1", ...
    "d_detwt2","d_detxRaw","d_detyRaw","d_detwRaw","d_sumX", ...
    "d_sumY","d_sumW","d_tof","d_quality"];
values = {
    double(strNow.detxt1) - double(strPrev.detxt1)
    double(strNow.detxt2) - double(strPrev.detxt2)
    double(strNow.detyt1) - double(strPrev.detyt1)
    double(strNow.detyt2) - double(strPrev.detyt2)
    double(strNow.detwt1) - double(strPrev.detwt1)
    double(strNow.detwt2) - double(strPrev.detwt2)
    double(strNow.detxRaw) - double(strPrev.detxRaw)
    double(strNow.detyRaw) - double(strPrev.detyRaw)
    double(strNow.detwRaw) - double(strPrev.detwRaw)
    sumXNow - sumXPrev
    sumYNow - sumYPrev
    sumWNow - sumWPrev
    double(strNow.tof) - double(strPrev.tof)
    double(strNow.quality) - double(strPrev.quality)
};
targets = struct('names', names, 'values', {values});
end

function y = signedUnsigned(x, nBits)
signBit = 2^(nBits - 1);
fullScale = 2^nBits;
y = x;
y(x >= signBit) = y(x >= signBit) - fullScale;
end

function y = detrendLocal(x)
x = double(x(:));
w = min(101, max(11, 2*floor(numel(x) / 1000) + 1));
y = x - movmean(x, w, 'omitnan');
end
