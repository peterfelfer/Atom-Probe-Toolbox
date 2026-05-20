function results = hitsSearchStrMultiPairs(hitsFile, strFile, options)
% HITSSEARCHSTRMULTIPAIRS Test HITS multi-hit pair deltas against STR pairs.
%
% results = hitsSearchStrMultiPairs(hitsFile, strFile)
%
% HITS MULTI events can contain multiple 8-byte payload blocks in a single
% a4-delimited event. STR provides one decoded channel-bearing row per event
% after filtering. This diagnostic tests several explicit STR pairing
% hypotheses around each HITS MULTI event and correlates HITS hit2-hit1 byte
% differences with STR channel differences.
%
% OPTIONS:
%   eventOffset  - HITS eventIdx -> STR eventIdx offset (default: 0)
%   sampleStride - stride through candidate HITS MULTI events (default: 1)
%   maxRows      - maximum sampled MULTI pairs (default: 50000)
%   verbose      - print top correlations (default: true)
%
% See also: hitsSearchStrDeltas, hitsStrEventMap
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.eventOffset (1,1) double = 0
    options.sampleStride (1,1) double = 1
    options.maxRows (1,1) double = 50000
    options.verbose (1,1) logical = true
end

[pairs, hitSummary] = sampleMultiPairs(hitsFile, options);
[strHits, ~] = strLoad(strFile);
strHits = strCalculatePositions(strHits);

strRowByEvent = zeros(max(double(strHits.eventIdx)) + ...
    abs(options.eventOffset) + 3, 1, 'uint32');
strRowByEvent(double(strHits.eventIdx)) = uint32((1:height(strHits))');

predictors = buildPairPredictors(pairs);
[pairTargets, pairSummary] = buildPairTargets(pairs.eventIdx, strHits, ...
    strRowByEvent, options.eventOffset);

predictorOut = strings(0, 1);
pairingOut = strings(0, 1);
targetOut = strings(0, 1);
rhoOut = zeros(0, 1);
rhoDetrendedOut = zeros(0, 1);
nOut = zeros(0, 1);

for pi = 1:numel(predictors.names)
    x = predictors.values{pi};
    for gi = 1:numel(pairTargets.names)
        targetGroup = pairTargets.groups{gi};
        for ti = 1:numel(targetGroup.names)
            y = targetGroup.values{ti};
            ok = isfinite(x) & isfinite(y);
            if nnz(ok) < 100
                continue;
            end
            rho = corr(x(ok), y(ok), 'Rows', 'complete');
            rhoDetrended = corr(detrendLocal(x(ok)), detrendLocal(y(ok)), ...
                'Rows', 'complete');

            predictorOut(end+1, 1) = predictors.names(pi); %#ok<AGROW>
            pairingOut(end+1, 1) = pairTargets.names(gi); %#ok<AGROW>
            targetOut(end+1, 1) = targetGroup.names(ti); %#ok<AGROW>
            rhoOut(end+1, 1) = rho; %#ok<AGROW>
            rhoDetrendedOut(end+1, 1) = rhoDetrended; %#ok<AGROW>
            nOut(end+1, 1) = nnz(ok); %#ok<AGROW>
        end
    end
end

results = table(predictorOut, pairingOut, targetOut, rhoOut, abs(rhoOut), ...
    rhoDetrendedOut, abs(rhoDetrendedOut), nOut, ...
    'VariableNames', {'predictor','strPairing','target','correlation', ...
                      'absCorrelation','detrendedCorrelation', ...
                      'absDetrendedCorrelation','n'});

sortScore = results.absDetrendedCorrelation;
sortScore(isnan(sortScore)) = -Inf;
[~, ord] = sort(sortScore, 'descend');
results = results(ord, :);

if options.verbose
    fprintf('\n--- HITS/STR multi-hit pair search ---\n');
    fprintf('Event offset: %+d\n', options.eventOffset);
    fprintf('Candidate HITS events with >=2 payload blocks: %d\n', ...
        hitSummary.nCandidateEvents);
    fprintf('Sampled HITS pairs: %d\n', height(pairs));
    disp(pairSummary);
    disp(results(1:min(30, height(results)), :));
end
end

function [pairs, summary] = sampleMultiPairs(hitsFile, options)
[events, ~, raw] = hitsScanEvents(hitsFile);
eventFields = double(events.nFields);
completeBlocks = floor(max(0, eventFields - 1) / 2);
candidate = events.type == 'MULTI' & completeBlocks >= 2;
eventIdxAll = find(candidate);

eventIdxAll = eventIdxAll(1:max(1, round(options.sampleStride)):end);
if numel(eventIdxAll) > options.maxRows
    eventIdxAll = eventIdxAll(1:options.maxRows);
end

n = numel(eventIdxAll);
B1 = zeros(n, 8, 'uint8');
B2 = zeros(n, 8, 'uint8');
for k = 1:n
    ei = eventIdxAll(k);
    payloadOffset = double(events.byteOffset(ei)) + 4;
    B1(k, :) = raw.bytes(payloadOffset+1:payloadOffset+8);
    B2(k, :) = raw.bytes(payloadOffset+9:payloadOffset+16);
end

pairs = table(uint32(eventIdxAll(:)), B1(:,1), B1(:,2), B1(:,3), B1(:,4), ...
    B1(:,5), B1(:,6), B1(:,7), B1(:,8), B2(:,1), B2(:,2), B2(:,3), ...
    B2(:,4), B2(:,5), B2(:,6), B2(:,7), B2(:,8), ...
    'VariableNames', {'eventIdx','h1b0','h1b1','h1b2','h1b3', ...
                      'h1b4','h1b5','h1b6','h1b7','h2b0','h2b1', ...
                      'h2b2','h2b3','h2b4','h2b5','h2b6','h2b7'});

summary = struct();
summary.nCandidateEvents = nnz(candidate);
summary.nSampledPairs = height(pairs);
end

function predictors = buildPairPredictors(pairs)
names = strings(0, 1);
values = {};

for k = 0:7
    h1 = double(pairs.(sprintf('h1b%d', k)));
    h2 = double(pairs.(sprintf('h2b%d', k)));
    d = h2 - h1;
    dm = mod(d + 128, 256) - 128;

    names(end+1, 1) = "d_b" + k; %#ok<AGROW>
    values{end+1, 1} = d; %#ok<AGROW>
    names(end+1, 1) = "dmod_b" + k; %#ok<AGROW>
    values{end+1, 1} = dm; %#ok<AGROW>
    names(end+1, 1) = "h2_b" + k; %#ok<AGROW>
    values{end+1, 1} = h2; %#ok<AGROW>
end

for k = 0:6
    h1 = double(pairs.(sprintf('h1b%d', k))) + ...
         256 * double(pairs.(sprintf('h1b%d', k+1)));
    h2 = double(pairs.(sprintf('h2b%d', k))) + ...
         256 * double(pairs.(sprintf('h2b%d', k+1)));
    names(end+1, 1) = "d_u16_b" + k + "b" + (k+1); %#ok<AGROW>
    values{end+1, 1} = h2 - h1; %#ok<AGROW>
    names(end+1, 1) = "dmod_u16_b" + k + "b" + (k+1); %#ok<AGROW>
    values{end+1, 1} = mod(h2 - h1 + 32768, 65536) - 32768; %#ok<AGROW>
end

predictors = struct('names', names, 'values', {values});
end

function [pairTargets, summary] = buildPairTargets(eventIdx, strHits, rowByEvent, eventOffset)
baseEvent = double(eventIdx) + eventOffset;
rowSame = eventRows(rowByEvent, baseEvent);
rowEventPlus1 = eventRows(rowByEvent, baseEvent + 1);
rowEventMinus1 = eventRows(rowByEvent, baseEvent - 1);
rowLoadedPlus1 = rowSame + 1;
rowLoadedMinus1 = rowSame - 1;

pairNames = ["same_to_loadedPlus1","same_to_eventPlus1", ...
    "loadedMinus1_to_same","eventMinus1_to_same"];
rowA = {rowSame, rowSame, rowLoadedMinus1, rowEventMinus1};
rowB = {rowLoadedPlus1, rowEventPlus1, rowSame, rowSame};

groups = cell(numel(pairNames), 1);
validCounts = zeros(numel(pairNames), 1);
for k = 1:numel(pairNames)
    valid = rowA{k} >= 1 & rowA{k} <= height(strHits) & ...
            rowB{k} >= 1 & rowB{k} <= height(strHits);
    validCounts(k) = nnz(valid);
    groups{k} = buildStrDeltaGroup(strHits, rowA{k}, rowB{k}, valid);
end

pairTargets = struct('names', pairNames, 'groups', {groups});
summary = table(pairNames(:), validCounts, ...
    'VariableNames', {'strPairing','validPairs'});
end

function rows = eventRows(rowByEvent, eventIdx)
rows = zeros(size(eventIdx));
inBounds = eventIdx >= 1 & eventIdx <= numel(rowByEvent);
rows(inBounds) = double(rowByEvent(eventIdx(inBounds)));
end

function group = buildStrDeltaGroup(strHits, rowA, rowB, valid)
names = ["d_detxt1","d_detxt2","d_detyt1","d_detyt2","d_detwt1", ...
    "d_detwt2","d_detxRaw","d_detyRaw","d_detwRaw","d_sumX", ...
    "d_sumY","d_sumW","d_tof","d_quality"];
values = cell(numel(names), 1);
for k = 1:numel(values)
    values{k} = nan(numel(rowA), 1);
end

if any(valid)
    A = strHits(rowA(valid), :);
    B = strHits(rowB(valid), :);
    sumXA = double(A.detxt1) + double(A.detxt2);
    sumYA = double(A.detyt1) + double(A.detyt2);
    sumWA = double(A.detwt1) + double(A.detwt2);
    sumXB = double(B.detxt1) + double(B.detxt2);
    sumYB = double(B.detyt1) + double(B.detyt2);
    sumWB = double(B.detwt1) + double(B.detwt2);

    deltas = {
        double(B.detxt1) - double(A.detxt1)
        double(B.detxt2) - double(A.detxt2)
        double(B.detyt1) - double(A.detyt1)
        double(B.detyt2) - double(A.detyt2)
        double(B.detwt1) - double(A.detwt1)
        double(B.detwt2) - double(A.detwt2)
        double(B.detxRaw) - double(A.detxRaw)
        double(B.detyRaw) - double(A.detyRaw)
        double(B.detwRaw) - double(A.detwRaw)
        sumXB - sumXA
        sumYB - sumYA
        sumWB - sumWA
        double(B.tof) - double(A.tof)
        double(B.quality) - double(A.quality)
    };
    for k = 1:numel(values)
        values{k}(valid) = deltas{k};
    end
end

group = struct('names', names, 'values', {values});
end

function y = detrendLocal(x)
x = double(x(:));
w = min(101, max(11, 2*floor(numel(x) / 1000) + 1));
y = x - movmean(x, w, 'omitnan');
end
