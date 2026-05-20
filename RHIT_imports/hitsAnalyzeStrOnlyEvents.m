function report = hitsAnalyzeStrOnlyEvents(hitsFile, strFile, options)
% HITSANALYZESTRONLYEVENTS Inspect HITS events where STR has channel data.
%
% report = hitsAnalyzeStrOnlyEvents(hitsFile, strFile)
%
% Finds event indices where STR has detector-channel data but HITS has no
% complete 8-byte ion payload according to hitsLoad/hitsStrEventMap. These
% events are the most likely hiding place for the remaining HITS-encoded
% ions. The function summarizes HITS event type, event size, raw payload
% length, tag distribution, and a byte-position fingerprint.
%
% OPTIONS:
%   maxEventsPerGroup - maximum events sampled per type/size group
%                       (default: 50000)
%   verbose           - print summaries (default: true)
%
% See also: hitsStrEventMap, hitsScanEvents, strScanEvents
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.maxEventsPerGroup (1,1) double = 50000
    options.verbose (1,1) logical = true
end

[hitEvents, hitHeader, hitRaw] = hitsScanEvents(hitsFile);
[strEvents, strHeader] = strScanEvents(strFile);

n = min(numel(hitEvents.fieldIdx), numel(strEvents.eventIdx));
payloadCount = hitsPayloadCounts(hitEvents);
strOnly = payloadCount(1:n) == 0 & strEvents.hasAnyChannel(1:n);

report = struct();
report.hitsFile = hitsFile;
report.strFile = strFile;
report.hitsHeader = hitHeader;
report.strHeader = strHeader;
report.nStrOnlyEvents = sum(strOnly);

typeStr = string(hitEvents.type(1:n));
nFields = double(hitEvents.nFields(1:n));
nCompleteDL = double(strEvents.nCompleteDelayLines(1:n));

[groupKey, groupType, groupFields, groupCompleteDL] = makeGroupKey( ...
    typeStr(strOnly), nFields(strOnly), nCompleteDL(strOnly));
[uKey, ~, ic] = unique(groupKey);
groupCount = accumarray(ic, 1);

summary = table();
summary.group = uKey;
summary.hitsType = groupType;
summary.nFields = groupFields;
summary.nCompleteDelayLines = groupCompleteDL;
summary.count = groupCount;
[~, ord] = sort(summary.count, 'descend');
summary = summary(ord, :);
report.groupSummary = summary;

details = struct();
for gi = 1:height(summary)
    sel = strOnly & typeStr == summary.hitsType(gi) & ...
        nFields == summary.nFields(gi) & ...
        nCompleteDL == summary.nCompleteDelayLines(gi);
    eventIdx = find(sel);
    if numel(eventIdx) > options.maxEventsPerGroup
        eventIdx = eventIdx(round(linspace(1, numel(eventIdx), options.maxEventsPerGroup)));
    end
    name = matlab.lang.makeValidName(sprintf('%s_nf%d_dl%d', ...
        summary.hitsType(gi), summary.nFields(gi), summary.nCompleteDelayLines(gi)));
    details.(name) = analyzeGroup(eventIdx, hitEvents, hitRaw, strEvents);
end
report.details = details;

if options.verbose
    printReport(report);
end
end

function counts = hitsPayloadCounts(hitEvents)
nf = double(hitEvents.nFields);
completeBlocks = floor(max(0, nf - 1) / 2);
counts = zeros(size(nf));
counts((hitEvents.type == 'SINGLE' | hitEvents.type == 'STATUS' | ...
    hitEvents.type == 'OTHER') & completeBlocks > 0) = 1;
isMultiPayload = hitEvents.type == 'MULTI' & completeBlocks > 0;
counts(isMultiPayload) = completeBlocks(isMultiPayload);
end

function [key, typeOut, fieldsOut, dlOut] = makeGroupKey(typeStr, nFields, nCompleteDL)
key = typeStr + "_nf" + string(nFields) + "_dl" + string(nCompleteDL);
[uKey, ia] = unique(key, 'stable');
typeOut = typeStr(ia);
fieldsOut = nFields(ia);
dlOut = nCompleteDL(ia);
key = key(:);
typeOut = typeOut(:);
fieldsOut = fieldsOut(:);
dlOut = dlOut(:);

% Return group metadata aligned to unique(key) order used by unique later.
[~, loc] = ismember(unique(key), uKey);
typeOut = typeOut(loc);
fieldsOut = fieldsOut(loc);
dlOut = dlOut(loc);
end

function detail = analyzeGroup(eventIdx, hitEvents, hitRaw, strEvents)
n = numel(eventIdx);
maxPayloadLength = max(max(0, double(hitEvents.nFields(eventIdx)) - 1) * 4);
payload = nan(n, maxPayloadLength);
payloadTag = nan(n, max(0, double(max(hitEvents.nFields(eventIdx))) - 1));
payloadVal = nan(n, max(0, double(max(hitEvents.nFields(eventIdx))) - 1));

for k = 1:n
    ei = eventIdx(k);
    startField = double(hitEvents.fieldIdx(ei)) + 1;
    nPayloadFields = double(hitEvents.nFields(ei)) - 1;
    if nPayloadFields <= 0
        continue;
    end
    byteStart = double(hitEvents.byteOffset(ei)) + 4;
    nBytes = nPayloadFields * 4;
    payload(k, 1:nBytes) = double(hitRaw.bytes(byteStart+1:byteStart+nBytes));
    payloadTag(k, 1:nPayloadFields) = double(hitRaw.tag(startField:startField+nPayloadFields-1));
    payloadVal(k, 1:nPayloadFields) = double(hitRaw.val(startField:startField+nPayloadFields-1));
end

byteStats = summarizeBytes(payload);
tagStats = summarizeTags(payloadTag);

detail = struct();
detail.nSampled = n;
detail.eventIdx = uint32(eventIdx(:));
detail.byteStats = byteStats;
detail.tagStats = tagStats;
detail.strQuality = summarizeVector(strEvents.quality(eventIdx));
detail.strCompleteDelayLines = summarizeVector(double(strEvents.nCompleteDelayLines(eventIdx)));
detail.strChannelMaskCounts = channelMaskCounts(strEvents.hasChannel(eventIdx, :));
detail.firstPayloadBytes = payload(1:min(10, n), :);
detail.firstPayloadTags = payloadTag(1:min(10, n), :);
detail.firstPayloadVals = payloadVal(1:min(10, n), :);
end

function T = summarizeBytes(payload)
if isempty(payload)
    T = table();
    return;
end

nPos = size(payload, 2);
bytePos = (1:nPos)';
nFinite = zeros(nPos, 1);
nUnique = zeros(nPos, 1);
meanVal = nan(nPos, 1);
stdVal = nan(nPos, 1);
minVal = nan(nPos, 1);
maxVal = nan(nPos, 1);

for p = 1:nPos
    x = payload(:, p);
    x = x(isfinite(x));
    nFinite(p) = numel(x);
    if isempty(x)
        continue;
    end
    nUnique(p) = numel(unique(x));
    meanVal(p) = mean(x);
    stdVal(p) = std(x);
    minVal(p) = min(x);
    maxVal(p) = max(x);
end

T = table(bytePos, nFinite, nUnique, meanVal, stdVal, minVal, maxVal);
end

function T = summarizeTags(payloadTag)
x = payloadTag(:);
x = x(isfinite(x));
if isempty(x)
    T = table();
    return;
end

[u, ~, ic] = unique(x);
count = accumarray(ic, 1);
[count, ord] = sort(count, 'descend');
tag = u(ord);
fraction = count / sum(count);
T = table(tag, count, fraction);
end

function s = summarizeVector(x)
x = double(x);
x = x(isfinite(x));
s = struct('n', numel(x), 'mean', mean(x), 'std', std(x), ...
    'min', min(x), 'median', median(x), 'max', max(x));
end

function T = channelMaskCounts(mask)
key = string(mask(:,1)) + string(mask(:,2)) + string(mask(:,3)) + ...
      string(mask(:,4)) + string(mask(:,5)) + string(mask(:,6));
[u, ~, ic] = unique(key);
count = accumarray(ic, 1);
[count, ord] = sort(count, 'descend');
maskKey = u(ord);
T = table(maskKey, count);
end

function printReport(report)
fprintf('\n--- HITS events with STR channel data but no complete 8-byte HITS payload ---\n');
fprintf('Total events: %d\n', report.nStrOnlyEvents);
fprintf('\nTop groups by HITS type / nFields / STR complete delay lines:\n');
disp(report.groupSummary(1:min(20, height(report.groupSummary)), :));

names = string(fieldnames(report.details));
for k = 1:min(8, numel(names))
    d = report.details.(char(names(k)));
    fprintf('\nGroup %s, sampled %d events\n', names(k), d.nSampled);
    fprintf('Top payload tags:\n');
    disp(d.tagStats(1:min(12, height(d.tagStats)), :));
    fprintf('Byte fingerprint:\n');
    disp(d.byteStats(1:min(16, height(d.byteStats)), :));
    fprintf('STR channel masks:\n');
    disp(d.strChannelMaskCounts(1:min(8, height(d.strChannelMaskCounts)), :));
end
end
