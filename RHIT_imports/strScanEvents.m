function [events, header, raw] = strScanEvents(fileName)
% STRSCANEVENTS Scan STR TLV records into event-level fingerprints.
%
% [events, header, raw] = strScanEvents(fileName)
%
% STR files use 4-byte TLV records and tag 0x18 as the event delimiter.
% Some AP Suite-era STR files carry a v3 header but retain STR-style bulk
% TLV channel records. This scanner preserves the original event index and
% summarizes channel presence plus selected status/counter tags.
%
% OUTPUTS:
%   events - struct with one row per original STR event:
%              eventIdx, fieldIdxEnd, byteOffsetEnd, nFields
%              hasChannel [nEvents x 6] for tags 0x01,0x02,0x03,0x04,0x21,0x22
%              completeDelayLines [nEvents x 3] for x,y,w pairs
%              nCompleteDelayLines, hasAnyChannel, quality
%              selected last-value tags: tag05, tag0b, tag49, tag1e, tag20
%   header - parsed header metadata
%   raw    - parsed TLV arrays and byte buffer
%
% See also: strLoad, hitsScanEvents, hitsStrEventMap
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
end

if isempty(fileName)
    [file, path] = uigetfile( ...
        {'*.STR;*.str', 'STR files (*.STR)'}, 'Select STR file');
    if isequal(file, 0)
        events = struct(); header = struct(); raw = struct();
        return;
    end
    fileName = fullfile(path, file);
end

fid = fopen(fileName, 'r', 'l');
if fid < 0
    error('strScanEvents:fileOpen', 'Cannot open %s', fileName);
end
cleanupObj = onCleanup(@() fclose(fid));

bytes = fread(fid, Inf, '*uint8');
nBytes = numel(bytes);
nFields = floor(nBytes / 4);
tail = nBytes - nFields * 4;
fprintf('strScanEvents: %s (%.1f MB, %d records, %d trailing bytes)\n', ...
    fileName, nBytes/1e6, nFields, tail);

rec = reshape(bytes(1:nFields*4), 4, [])';
b0 = rec(:, 1);
b1 = rec(:, 2);
b2 = rec(:, 3);
tag = rec(:, 4);

val = int32(b0) + int32(b1)*256 + int32(b2)*65536;
neg = val >= 8388608;
val(neg) = val(neg) - 16777216;

raw = struct('b0', b0, 'b1', b1, 'b2', b2, 'tag', tag, 'val', val, ...
             'bytes', bytes, 'nBytes', nBytes, 'tail', tail, ...
             'nFields', nFields);
header = parseHeader(tag, val);

isEnd = tag == 0x18;
endPos = find(isEnd);
if isempty(endPos)
    events = struct();
    warning('strScanEvents:noEvents', 'No STR 0x18 event delimiters found.');
    return;
end

eventId = cumsum(isEnd);
eventId = eventId - isEnd + 1;
nEvents = max(eventId);

prevEnd = [0; endPos(1:end-1)];
nFieldsPerEvent = uint16(endPos - prevEnd);

channelTags = uint8([0x01 0x02 0x03 0x04 0x21 0x22]);
hasChannel = false(nEvents, 6);
for ci = 1:6
    ev = eventId(tag == channelTags(ci));
    hasChannel(ev, ci) = true;
end

completeDelayLines = [hasChannel(:,1) & hasChannel(:,2), ...
                      hasChannel(:,3) & hasChannel(:,4), ...
                      hasChannel(:,5) & hasChannel(:,6)];
nCompleteDelayLines = uint8(sum(completeDelayLines, 2));
hasAnyChannel = any(hasChannel, 2);

events = struct();
events.eventIdx = uint32((1:nEvents)');
fieldIdxEnd = endPos;
if nEvents > numel(endPos)
    fieldIdxEnd(end+1:nEvents, 1) = nFields;
    nFieldsPerEvent(end+1:nEvents, 1) = uint16(nFields - endPos(end));
end
events.fieldIdxEnd = uint32(fieldIdxEnd);
events.byteOffsetEnd = uint64(fieldIdxEnd - 1) * uint64(4);
events.nFields = nFieldsPerEvent;
events.hasChannel = hasChannel;
events.completeDelayLines = completeDelayLines;
events.nCompleteDelayLines = nCompleteDelayLines;
events.hasAnyChannel = hasAnyChannel;
quality = nan(nEvents, 1);
quality(1:numel(endPos)) = double(val(isEnd));
events.quality = quality;

statusTags = uint8([0x05 0x0b 0x49 0x1e 0x20]);
statusNames = {'tag05', 'tag0b', 'tag49', 'tag1e', 'tag20'};
for si = 1:numel(statusTags)
    out = nan(nEvents, 1);
    mask = tag == statusTags(si);
    out(eventId(mask)) = double(val(mask));
    events.(statusNames{si}) = out;
end

fprintf(['strScanEvents: %d events, %d channel-bearing, ' ...
         '%d complete six-channel\n'], ...
    nEvents, sum(hasAnyChannel), sum(nCompleteDelayLines == 3));
end

function header = parseHeader(tags, vals)
header = struct();
if ~isempty(tags) && tags(1) == 0xa0
    header.version = bitand(vals(1), 255);
else
    header.version = 0;
end

labels = char(zeros(1, 6));
for i = 1:6
    idx = find(tags(1:min(50,numel(tags))) == uint8(0x26 + i), 1);
    if ~isempty(idx)
        labels(i) = char(vals(idx));
    end
end
header.detectorLabels = labels;

idx = find(tags(1:min(100,numel(tags))) == 0x1b, 1);
if ~isempty(idx)
    header.voltage = double(vals(idx));
end
end
