function [events, header, raw] = hitsScanEvents(fileName)
% HITSSCANEVENTS Scan a HITS v3 file and index its events.
%
% [events, header, raw] = hitsScanEvents(fileName)
%
% Walks the file, locates all event markers (tag 0xa4), classifies each
% event by byte2 of the marker (0x00 STATUS, 0x20 SINGLE, 0x40 MULTI, else
% OTHER), and reports the size of each event in 4-byte records.
%
% INPUTS:
%   fileName - path to .HITS file. If empty, opens a file dialog.
%
% OUTPUTS:
%   events - struct with arrays of length nEvents:
%              fieldIdx    - record index of the 0xa4 marker (1-based)
%              byteOffset  - byte offset of the marker (0-based)
%              typeByte    - byte2 of marker (uint8: 0x00/0x20/0x40/...)
%              type        - categorical: STATUS / SINGLE / MULTI / OTHER
%              nFields     - record gap to next event (event size in records)
%              markerVal16 - uint16 from bytes 0-1 of the marker
%   header - struct with parsed header fields (version, labels, voltage,
%            thresholds, walkCorrections, timingChannels, ...).
%   raw    - struct with parsed TLV arrays for downstream payload work:
%              b0, b1, b2  - per-record bytes 0..2 (uint8)
%              tag         - per-record tag byte (uint8)
%              val         - signed 24-bit value per record (int32)
%              nBytes, tail - file size and trailing-byte count
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
end

if isempty(fileName)
    [file, path] = uigetfile( ...
        {'*.HITS;*.hits', 'HITS files (*.HITS)'}, 'Select HITS file');
    if isequal(file, 0)
        events = struct(); header = struct(); raw = struct();
        return;
    end
    fileName = fullfile(path, file);
end

%% Read bytes
fid = fopen(fileName, 'r', 'l');
if fid < 0
    error('hitsScanEvents:fileOpen', 'Cannot open %s', fileName);
end
cleanupObj = onCleanup(@() fclose(fid));

bytes   = fread(fid, Inf, '*uint8');
nBytes  = numel(bytes);
nFields = floor(nBytes / 4);
tail    = nBytes - nFields*4;
fprintf('hitsScanEvents: %s (%.1f MB, %d records, %d trailing bytes)\n', ...
    fileName, nBytes/1e6, nFields, tail);

%% Parse 4-byte TLV records: [b0 b1 b2 tag]
rec = reshape(bytes(1:nFields*4), 4, [])';   % nFields x 4
b0  = rec(:,1);
b1  = rec(:,2);
b2  = rec(:,3);
tag = rec(:,4);

% Signed 24-bit value (sign-extend bit 23)
val = int32(b0) + int32(b1)*256 + int32(b2)*65536;
neg = val >= 8388608;
val(neg) = val(neg) - 16777216;

raw = struct('b0', b0, 'b1', b1, 'b2', b2, 'tag', tag, 'val', val, ...
             'bytes', bytes, ...
             'nBytes', nBytes, 'tail', tail, 'nFields', nFields);

%% Verify magic (tag=0xa0; byte0 = version, expect 3 for HITS v3)
if tag(1) ~= 0xa0 || b0(1) ~= 3
    warning('hitsScanEvents:magic', ...
        'Magic mismatch: expected tag=0xa0 byte0=3, got tag=0x%02x byte0=%d', ...
        tag(1), b0(1));
end

%% Locate event markers
isA4   = (tag == 0xa4);
a4Pos  = find(isA4);
nEvents = numel(a4Pos);
fprintf('hitsScanEvents: %d 0xa4 markers found\n', nEvents);

if nEvents == 0
    events = struct(); header = parseHeader(tag, val, nFields); return;
end

%% Header = everything before the first 0xa4 marker
headerEnd = a4Pos(1) - 1;
header    = parseHeader(tag(1:headerEnd), val(1:headerEnd), headerEnd);

%% Per-event index
typeByte    = b2(a4Pos);
markerVal16 = uint16(b0(a4Pos)) + bitshift(uint16(b1(a4Pos)), 8);

% Gap to next event (= event size in records). Last event extends to EOF.
nextPos   = [a4Pos(2:end); int32(nFields) + 1];
gapFields = int32(nextPos) - int32(a4Pos);

% Classify by byte2 of the marker
typeName = repmat("OTHER", nEvents, 1);
typeName(typeByte == 0x00) = "STATUS";
typeName(typeByte == 0x20) = "SINGLE";
typeName(typeByte == 0x40) = "MULTI";
typeCat  = categorical(typeName, ["STATUS","SINGLE","MULTI","OTHER"]);

events = struct();
events.fieldIdx    = a4Pos;
events.byteOffset  = (a4Pos - 1) * 4;
events.typeByte    = typeByte;
events.type        = typeCat;
events.nFields     = gapFields;
events.markerVal16 = markerVal16;

%% Console summary
nStatus = sum(typeCat == "STATUS");
nSingle = sum(typeCat == "SINGLE");
nMulti  = sum(typeCat == "MULTI");
nOther  = sum(typeCat == "OTHER");
fprintf('hitsScanEvents: STATUS=%d  SINGLE=%d  MULTI=%d  OTHER=%d  TOTAL=%d\n', ...
    nStatus, nSingle, nMulti, nOther, nEvents);
end


function header = parseHeader(tag, val, nHeaderFields)
% Parse the header region (records before the first 0xa4 marker).
header = struct();

% 0xa0 = version. Only byte0 encodes the version number; bytes 1-2 form
% an additional 16-bit signature (typically 0x8010) that we expose as well.
idx = find(tag == 0xa0, 1);
if ~isempty(idx)
    rawVal = double(val(idx));
    if rawVal < 0, rawVal = rawVal + 16777216; end   % unsigned 24-bit
    header.version    = mod(rawVal, 256);
    header.versionSig = floor(rawVal / 256);          % bytes 1-2
end

% Detector channel labels: tags 0x27..0x2c -> 6 ASCII chars
labels = char(zeros(1, 6));
for k = 1:6
    idx = find(tag == uint8(0x26 + k), 1);
    if ~isempty(idx)
        labels(k) = char(mod(double(val(idx)), 256));
    end
end
header.detectorLabels = labels;

% Detection thresholds: 0x2d..0x2f
thr = nan(1, 3);
for k = 1:3
    idx = find(tag == uint8(0x2c + k), 1);
    if ~isempty(idx), thr(k) = double(val(idx)); end
end
header.thresholds = thr;

% Walk corrections: 0x33..0x38
walk = nan(1, 6);
for k = 1:6
    idx = find(tag == uint8(0x32 + k), 1);
    if ~isempty(idx), walk(k) = double(val(idx)); end
end
header.walkCorrections = walk;

% Five 48-bit timing channels: lo=0x3b..0x3f, hi=0x4c..0x50
loTags = uint8([0x3b 0x3c 0x3d 0x3e 0x3f]);
hiTags = uint8([0x4c 0x4d 0x4e 0x4f 0x50]);
chan = nan(1, 5);
for k = 1:5
    iLo = find(tag == loTags(k), 1);
    if ~isempty(iLo)
        lo = double(val(iLo));
        iHi = find(tag == hiTags(k), 1);
        if ~isempty(iHi)
            hi = double(val(iHi));
            chan(k) = hi * 16777216 + lo;
        else
            chan(k) = lo;
        end
    end
end
header.timingChannels = chan;

% Initial voltage (0x1b) and multi-hit indicator (0xa8)
idx = find(tag == 0x1b, 1);
if ~isempty(idx), header.voltage = double(val(idx)); end
idx = find(tag == 0xa8, 1);
if ~isempty(idx), header.multiHitIndicator = double(val(idx)); end

% Other documented unknown header fields
extraTags = uint8([0xa1 0xa2 0xa3 0x09 0x1e 0x49 0x40 0x13 0x14 0x15 ...
                   0x20 0x23 0x24 0x25 0x47 0x4b]);
for k = 1:numel(extraTags)
    idx = find(tag == extraTags(k), 1);
    if ~isempty(idx)
        header.(sprintf('tag_%02x', extraTags(k))) = double(val(idx));
    end
end

header.nHeaderFields = nHeaderFields;
end
