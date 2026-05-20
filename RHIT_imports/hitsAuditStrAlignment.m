function report = hitsAuditStrAlignment(hitsFile, strFile, options)
% HITSAUDITSTRALIGNMENT Audit same-run HITS and STR raw detector streams.
%
% report = hitsAuditStrAlignment(hitsFile, strFile)
%
% This is the gatekeeper for HITS bit-packing reverse engineering when a
% same-run STR companion exists. STR provides decoded delay-line timings;
% HITS provides the packed b0..b7 payload. Before fitting any bit mapping,
% the two event streams must be aligned.
%
% The function:
%   1. loads HITS payloads and the full HITS event index,
%   2. loads STR TLV-channel events while preserving their original eventIdx,
%   3. tests simple event-index offsets on sampled HITS payload rows, and
%   4. reports raw-byte correlations for the best offset as diagnostics.
%
% It does not claim a HITS decoder. Low raw-byte correlations are expected
% if the payload is bit-packed or conditional on mode bytes.
%
% OPTIONS:
%   sampleStride       - stride through HITS payload rows (default: 1000)
%   offsetCandidates   - event-index offsets to test (default: -2000:2000)
%   calculatePositions - add STR detxRaw/detyRaw/tof targets (default: false)
%   verbose            - print report (default: true)
%
% See also: hitsLoad, strLoad, strCalculatePositions
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.sampleStride (1,1) double = 1000
    options.offsetCandidates (:,1) double = (-2000:2000)'
    options.calculatePositions (1,1) logical = false
    options.verbose (1,1) logical = true
end

if ~exist(hitsFile, 'file')
    error('hitsAuditStrAlignment:hitsMissing', ...
        'HITS file does not exist: %s', hitsFile);
end
if ~exist(strFile, 'file')
    error('hitsAuditStrAlignment:strMissing', ...
        'STR file does not exist: %s', strFile);
end

report = struct();
report.hitsFile = hitsFile;
report.strFile = strFile;

[hits, hitsHeader, hitEvents] = hitsLoad(hitsFile);
[strHits, strMeta] = strLoad(strFile);
if options.calculatePositions
    strHits = strCalculatePositions(strHits);
end

report.hitsHeader = hitsHeader;
report.strMetadata = strMeta;
report.nHitsPayloads = height(hits);
report.nHitsEvents = numel(hitEvents.fieldIdx);
report.nStrEventsLoaded = height(strHits);
report.nStrOriginalEvents = max(double(strHits.eventIdx));
report.nEventCountDiff = report.nHitsEvents - report.nStrOriginalEvents;
report.nPayloadCountDiff = report.nStrEventsLoaded - report.nHitsPayloads;

idx = 1:max(1, round(options.sampleStride)):height(hits);
sample = hits(idx, :);
report.nSampledPayloads = height(sample);

maxEvent = max([double(sample.eventIdx); double(strHits.eventIdx)]);
strRowByEvent = zeros(maxEvent + max(abs(options.offsetCandidates)) + 2, 1, 'uint32');
strRowByEvent(double(strHits.eventIdx)) = uint32((1:height(strHits))');

offsetReport = table();
offsets = options.offsetCandidates(:);
validFraction = zeros(numel(offsets), 1);
medianEventDistance = zeros(numel(offsets), 1);
maxAbsByteCorr = zeros(numel(offsets), 1);

targetNames = availableTargetNames(strHits);
for oi = 1:numel(offsets)
    off = offsets(oi);
    eventQuery = double(sample.eventIdx) + off;
    inBounds = eventQuery >= 1 & eventQuery <= numel(strRowByEvent);
    rows = zeros(height(sample), 1, 'uint32');
    rows(inBounds) = strRowByEvent(eventQuery(inBounds));
    valid = rows > 0;
    validFraction(oi) = mean(valid);
    if any(valid)
        medianEventDistance(oi) = median(double(sample.eventIdx(valid)) + off - ...
            double(strHits.eventIdx(double(rows(valid)))));
        corrVals = abs(byteTargetCorrelations( ...
            sample(valid, :), strHits(double(rows(valid)), :), targetNames));
        maxAbsByteCorr(oi) = max(corrVals, [], 'omitnan');
    else
        medianEventDistance(oi) = NaN;
        maxAbsByteCorr(oi) = NaN;
    end
end

offsetReport.offset = offsets;
offsetReport.validFraction = validFraction;
offsetReport.medianEventDistance = medianEventDistance;
offsetReport.maxAbsByteCorr = maxAbsByteCorr;
[~, bestIdx] = sortrows([validFraction, maxAbsByteCorr], [-1 -2]);
bestOffset = offsets(bestIdx(1));

report.offsetReport = offsetReport;
report.bestOffset = bestOffset;
report.bestValidFraction = validFraction(bestIdx(1));
report.bestMaxAbsByteCorr = maxAbsByteCorr(bestIdx(1));

eventQuery = double(sample.eventIdx) + bestOffset;
inBounds = eventQuery >= 1 & eventQuery <= numel(strRowByEvent);
rows = zeros(height(sample), 1, 'uint32');
rows(inBounds) = strRowByEvent(eventQuery(inBounds));
valid = rows > 0;
report.rawByteCorrelations = quickByteCorrelations( ...
    sample(valid, :), strHits(double(rows(valid)), :), targetNames);

if options.verbose
    printReport(report);
end
end

function names = availableTargetNames(strHits)
candidateNames = ["detxt1","detxt2","detyt1","detyt2","detwt1","detwt2", ...
                  "quality","detxRaw","detyRaw","detwRaw","tof","hitType"];
names = candidateNames(ismember(candidateNames, string(strHits.Properties.VariableNames)));
end

function rho = byteTargetCorrelations(hits, strHits, targetNames)
byteNames = "b" + string(0:7);
rho = zeros(numel(byteNames) * numel(targetNames), 1);
ri = 0;
for bi = 1:numel(byteNames)
    x = double(hits.(char(byteNames(bi))));
    for ti = 1:numel(targetNames)
        y = double(strHits.(char(targetNames(ti))));
        ri = ri + 1;
        rho(ri) = corr(x, y, 'Rows', 'complete');
    end
end
end

function C = quickByteCorrelations(hits, strHits, targetNames)
byteNames = "b" + string(0:7);
rawByte = strings(0, 1);
target = strings(0, 1);
rho = zeros(0, 1);

for bi = 1:numel(byteNames)
    x = double(hits.(char(byteNames(bi))));
    for ti = 1:numel(targetNames)
        y = double(strHits.(char(targetNames(ti))));
        rawByte(end+1, 1) = byteNames(bi); %#ok<AGROW>
        target(end+1, 1) = targetNames(ti); %#ok<AGROW>
        rho(end+1, 1) = corr(x, y, 'Rows', 'complete'); %#ok<AGROW>
    end
end

C = table(rawByte, target, rho, ...
    'VariableNames', {'rawByte', 'target', 'correlation'});
end

function printReport(report)
fprintf('\n--- HITS/STR alignment audit ---\n');
fprintf('HITS events: %d\n', report.nHitsEvents);
fprintf('STR original events: %d\n', report.nStrOriginalEvents);
fprintf('Event count diff HITS - STR: %+d\n', report.nEventCountDiff);
fprintf('HITS payload rows: %d\n', report.nHitsPayloads);
fprintf('STR loaded events: %d\n', report.nStrEventsLoaded);
fprintf('Payload count diff STR - HITS: %+d\n', report.nPayloadCountDiff);
fprintf('Sampled HITS payload rows: %d\n', report.nSampledPayloads);
fprintf('Best simple event offset: %+d\n', report.bestOffset);
fprintf('Best valid sample fraction: %.4f\n', report.bestValidFraction);
fprintf('Best max |raw-byte correlation|: %.4f\n', report.bestMaxAbsByteCorr);

fprintf('\nTop simple-offset candidates:\n');
T = report.offsetReport;
[~, ord] = sortrows([T.validFraction, T.maxAbsByteCorr], [-1 -2]);
disp(T(ord(1:min(10, height(T))), :));

if ~isempty(report.rawByteCorrelations)
    fprintf('\nRaw-byte correlations at best offset, sorted by |rho|:\n');
    C = report.rawByteCorrelations;
    [~, ord] = sort(abs(C.correlation), 'descend');
    disp(C(ord(1:min(16, height(C))), :));
end
end
