function eventMap = hitsStrEventMap(hitsFile, strFile, options)
% HITSSTREVENTMAP Build a shared event-axis map for HITS and STR files.
%
% eventMap = hitsStrEventMap(hitsFile, strFile)
%
% The map is designed to diagnose stream alignment before bit-packing
% inference. It compares HITS 0xa4-delimited event classes with STR
% 0x18-delimited channel-bearing events and extracts selected status
% counter anchors from both streams.
%
% OPTIONS:
%   'verbose' - print summary (default: true)
%
% See also: hitsScanEvents, strScanEvents
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.verbose (1,1) logical = true
end

[hitEvents, hitHeader, hitRaw] = hitsScanEvents(hitsFile);
[strEvents, strHeader] = strScanEvents(strFile);

nHitsEvents = numel(hitEvents.fieldIdx);
nStrEvents = numel(strEvents.eventIdx);
n = max(nHitsEvents, nStrEvents);

eventIdx = uint32((1:n)');
hitsType = categorical(repmat("MISSING", n, 1), ...
    ["STATUS","SINGLE","MULTI","OTHER","MISSING"]);
hitsType(1:nHitsEvents) = addcats(hitEvents.type, "MISSING");

hitsPayloadCount = zeros(n, 1, 'uint16');
hitsPayloadCount(1:nHitsEvents) = uint16(hitsPayloadCounts(hitEvents));
hitsHasPayload = hitsPayloadCount > 0;
hitsNFields = zeros(n, 1, 'uint16');
hitsNFields(1:nHitsEvents) = uint16(hitEvents.nFields);

strHasAnyChannel = false(n, 1);
strNCompleteDelayLines = zeros(n, 1, 'uint8');
strNFields = zeros(n, 1, 'uint16');
strQuality = nan(n, 1);
if nStrEvents > 0
    strHasAnyChannel(1:nStrEvents) = strEvents.hasAnyChannel(:);
    strNCompleteDelayLines(1:nStrEvents) = strEvents.nCompleteDelayLines(:);
    strNFields(1:nStrEvents) = strEvents.nFields(:);
    strQuality(1:nStrEvents) = strEvents.quality(:);
end

statusTags = ["tag05","tag0b","tag49","tag1e","tag20"];
hitsStatus = extractHitsStatus(hitEvents, hitRaw, statusTags);

strStatus = struct();
for k = 1:numel(statusTags)
    vals = nan(n, 1);
    vals(1:nStrEvents) = strEvents.(char(statusTags(k)));
    strStatus.(char(statusTags(k))) = vals;
end

eventMap = struct();
eventMap.hitsFile = hitsFile;
eventMap.strFile = strFile;
eventMap.hitsHeader = hitHeader;
eventMap.strHeader = strHeader;
eventMap.eventIdx = eventIdx;
eventMap.hitsType = hitsType;
eventMap.hitsNFields = hitsNFields;
eventMap.hitsPayloadCount = hitsPayloadCount;
eventMap.hitsHasPayload = hitsHasPayload;
eventMap.strNFields = strNFields;
eventMap.strHasAnyChannel = strHasAnyChannel;
eventMap.strNCompleteDelayLines = strNCompleteDelayLines;
eventMap.strQuality = strQuality;
eventMap.hitsStatus = hitsStatus;
eventMap.strStatus = strStatus;
eventMap.summary = summarizeMap(eventMap, nHitsEvents, nStrEvents);

if options.verbose
    printSummary(eventMap.summary);
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

function status = extractHitsStatus(hitEvents, raw, statusTags)
n = numel(hitEvents.fieldIdx);
status = struct();
for k = 1:numel(statusTags)
    status.(char(statusTags(k))) = nan(n, 1);
end

tagCodes = uint8([0x05 0x0b 0x49 0x1e 0x20]);
for ei = 1:n
    startField = double(hitEvents.fieldIdx(ei)) + 1;
    endField = startField + double(hitEvents.nFields(ei)) - 2;
    if endField < startField || startField > raw.nFields
        continue;
    end
    endField = min(endField, raw.nFields);

    tags = raw.tag(startField:endField);
    vals = raw.val(startField:endField);
    for k = 1:numel(tagCodes)
        idx = find(tags == tagCodes(k), 1, 'last');
        if ~isempty(idx)
            status.(char(statusTags(k)))(ei) = double(vals(idx));
        end
    end
end
end

function summary = summarizeMap(eventMap, nHitsEvents, nStrEvents)
commonN = min(nHitsEvents, nStrEvents);
summary = struct();
summary.nHitsEvents = nHitsEvents;
summary.nStrEvents = nStrEvents;
summary.eventCountDiff = nHitsEvents - nStrEvents;
summary.nHitsPayloadEvents = sum(eventMap.hitsHasPayload);
summary.nHitsPayloadRows = sum(double(eventMap.hitsPayloadCount));
summary.nStrChannelEvents = sum(eventMap.strHasAnyChannel);
summary.nStrCompleteSixChannelEvents = sum(eventMap.strNCompleteDelayLines == 3);
summary.strMinusHitsPayloadRows = summary.nStrChannelEvents - summary.nHitsPayloadRows;

summary.commonEvents = commonN;
summary.bothPayloadAndChannel = sum(eventMap.hitsHasPayload(1:commonN) & ...
    eventMap.strHasAnyChannel(1:commonN));
summary.hitsPayloadOnly = sum(eventMap.hitsHasPayload(1:commonN) & ...
    ~eventMap.strHasAnyChannel(1:commonN));
summary.strChannelOnly = sum(~eventMap.hitsHasPayload(1:commonN) & ...
    eventMap.strHasAnyChannel(1:commonN));
summary.neitherPayloadNorChannel = sum(~eventMap.hitsHasPayload(1:commonN) & ...
    ~eventMap.strHasAnyChannel(1:commonN));

summary.hitsTypeCounts = table(categories(eventMap.hitsType), ...
    countcats(eventMap.hitsType), ...
    'VariableNames', {'hitsType','count'});
summary.strCompleteDelayLineCounts = table((0:3)', ...
    arrayfun(@(k) sum(eventMap.strNCompleteDelayLines == k), (0:3)'), ...
    'VariableNames', {'nCompleteDelayLines','count'});
summary.hitsTypeVsStrChannel = hitsTypeVsStrChannel(eventMap);
summary.strChannelOnlyByHitsType = strChannelOnlyByHitsType(eventMap, commonN);
summary.strChannelOnlyByCompleteDL = strChannelOnlyByCompleteDL(eventMap, commonN);

summary.statusAnchorMatches = statusAnchorSummary(eventMap, commonN);
end

function T = hitsTypeVsStrChannel(eventMap)
cats = categories(eventMap.hitsType);
hitsType = strings(numel(cats), 1);
nEvents = zeros(numel(cats), 1);
nStrChannel = zeros(numel(cats), 1);
nHitsPayload = zeros(numel(cats), 1);
nBoth = zeros(numel(cats), 1);

for k = 1:numel(cats)
    sel = eventMap.hitsType == cats{k};
    hitsType(k) = string(cats{k});
    nEvents(k) = sum(sel);
    nStrChannel(k) = sum(sel & eventMap.strHasAnyChannel);
    nHitsPayload(k) = sum(sel & eventMap.hitsHasPayload);
    nBoth(k) = sum(sel & eventMap.strHasAnyChannel & eventMap.hitsHasPayload);
end

T = table(hitsType, nEvents, nStrChannel, nHitsPayload, nBoth);
end

function T = strChannelOnlyByHitsType(eventMap, commonN)
cats = categories(eventMap.hitsType);
hitsType = strings(numel(cats), 1);
count = zeros(numel(cats), 1);
selBase = ~eventMap.hitsHasPayload(1:commonN) & eventMap.strHasAnyChannel(1:commonN);

for k = 1:numel(cats)
    hitsType(k) = string(cats{k});
    count(k) = sum(selBase & eventMap.hitsType(1:commonN) == cats{k});
end

T = table(hitsType, count);
end

function T = strChannelOnlyByCompleteDL(eventMap, commonN)
nCompleteDelayLines = (0:3)';
count = zeros(4, 1);
selBase = ~eventMap.hitsHasPayload(1:commonN) & eventMap.strHasAnyChannel(1:commonN);

for k = 0:3
    count(k+1) = sum(selBase & eventMap.strNCompleteDelayLines(1:commonN) == k);
end

T = table(nCompleteDelayLines, count);
end

function T = statusAnchorSummary(eventMap, commonN)
names = string(fieldnames(eventMap.hitsStatus));
nameCol = strings(numel(names), 1);
nBoth = zeros(numel(names), 1);
fractionEqual = nan(numel(names), 1);
correlation = nan(numel(names), 1);

for k = 1:numel(names)
    h = eventMap.hitsStatus.(char(names(k)))(1:commonN);
    s = eventMap.strStatus.(char(names(k)))(1:commonN);
    both = isfinite(h) & isfinite(s);
    nameCol(k) = names(k);
    nBoth(k) = sum(both);
    if nBoth(k) > 0
        fractionEqual(k) = mean(h(both) == s(both));
    end
    if nBoth(k) > 2
        correlation(k) = corr(h(both), s(both), 'Rows', 'complete');
    end
end

T = table(nameCol, nBoth, fractionEqual, correlation, ...
    'VariableNames', {'tag','nBoth','fractionEqual','correlation'});
end

function printSummary(summary)
fprintf('\n--- HITS/STR event map ---\n');
fprintf('HITS events: %d\n', summary.nHitsEvents);
fprintf('STR events: %d\n', summary.nStrEvents);
fprintf('Event count diff HITS - STR: %+d\n', summary.eventCountDiff);
fprintf('HITS payload events: %d\n', summary.nHitsPayloadEvents);
fprintf('HITS payload rows: %d\n', summary.nHitsPayloadRows);
fprintf('STR channel events: %d\n', summary.nStrChannelEvents);
fprintf('STR six-channel events: %d\n', summary.nStrCompleteSixChannelEvents);
fprintf('STR channel events - HITS payload rows: %+d\n', ...
    summary.strMinusHitsPayloadRows);

fprintf('\nCommon event-axis classification:\n');
fprintf('  both HITS payload and STR channel: %d\n', summary.bothPayloadAndChannel);
fprintf('  HITS payload only: %d\n', summary.hitsPayloadOnly);
fprintf('  STR channel only: %d\n', summary.strChannelOnly);
fprintf('  neither: %d\n', summary.neitherPayloadNorChannel);

fprintf('\nHITS type counts:\n');
disp(summary.hitsTypeCounts);
fprintf('HITS type vs STR channel/payload:\n');
disp(summary.hitsTypeVsStrChannel);
fprintf('STR-channel-only events by HITS type:\n');
disp(summary.strChannelOnlyByHitsType);
fprintf('STR-channel-only events by complete delay-line count:\n');
disp(summary.strChannelOnlyByCompleteDL);
fprintf('STR complete-delay-line counts:\n');
disp(summary.strCompleteDelayLineCounts);
fprintf('Status anchor overlap at same event index:\n');
disp(summary.statusAnchorMatches);
end
