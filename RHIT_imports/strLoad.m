function [hits, metadata] = strLoad(fileName)
% STRLOAD Load raw detector data from a Cameca STR or HITS file.
%
% [hits, metadata] = strLoad(fileName)
%
% The STR format (version 2) and HITS format (version 3) are compact
% binary formats used by Cameca LEAP atom probes. They store raw TDC
% timing data from the delay-line detector (6 channels). This function
% decodes the binary data directly in MATLAB — no Python required.
%
% INPUTS:
%   fileName - Path to the .STR or .HITS file. If empty, opens a dialog.
%
% OUTPUTS:
%   hits     - Table with per-event raw delay-line timing data:
%                ionIdx   - sequential loaded-hit index (1-based)
%                eventIdx - original STR event index before empty-marker
%                           filtering
%                detxt1, detxt2 - DL-x end timings (TDC counts)
%                detyt1, detyt2 - DL-y end timings (TDC counts)
%                detwt1, detwt2 - DL-w end timings (45° diagonal)
%                quality  - hit finding quality / chi2
%              Events with missing channels have NaN for those fields.
%              Use strCalculatePositions to compute detxRaw, detyRaw, tof.
%
%   metadata - Struct with file header information:
%                version        - format version (2=STR, 3=HITS)
%                detectorLabels - 6-char string of channel labels
%                thresholds     - [3x1] detection thresholds
%                walkCorrections - [6x1] walk correction values
%                voltage        - initial specimen voltage (V)
%                nEvents        - total event count
%
% FORMAT:
%   The file consists of 4-byte records: [int24_LE value][uint8 tag].
%   Per-event records use tags 0x01-0x04, 0x21-0x22 (6 timing channels)
%   and tag 0x18 (quality, marks end of event).
%
% EXAMPLE:
%   [hits, meta] = strLoad('data.STR');
%   histogram(hits.tof, 1000); xlabel('TOF (TDC counts)');
%
% See also: rhitLoad, posLoad
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
end

if isempty(fileName)
    [file, path] = uigetfile( ...
        {'*.STR;*.str;*.HITS;*.hits', 'STR/HITS files (*.STR, *.HITS)'}, ...
        'Select STR or HITS file');
    if isequal(file, 0)
        hits = table();
        metadata = struct();
        return;
    end
    fileName = fullfile(path, file);
end

%% Read file
fid = fopen(fileName, 'r', 'l');  % little-endian
if fid < 0
    error('strLoad:fileOpen', 'Cannot open %s', fileName);
end
cleanupObj = onCleanup(@() fclose(fid));

raw = fread(fid, Inf, '*uint8');
nBytes = numel(raw);
nFields = floor(nBytes / 4);
fprintf('Reading %s (%.1f MB, %d fields)\n', fileName, nBytes/1e6, nFields);

%% Parse 4-byte TLV records: [int24_LE][uint8_tag]
raw4 = reshape(raw(1:nFields*4), 4, [])';
tagVec  = raw4(:, 4);
valVec  = int32(raw4(:,1)) + int32(raw4(:,2)) * 256 + int32(raw4(:,3)) * 65536;
valVec(valVec >= 8388608) = valVec(valVec >= 8388608) - 16777216; % sign extend 24-bit

%% Parse header
metadata = parseHeader(tagVec, valVec);

%% Extract events — format depends on version
% STR (v2): bulk data uses compact tags 0x01-0x04, 0x21-0x22
%           tag 0x18 = quality, marks end of each event
% HITS (v3): different bulk encoding (delimited by 0xa4 markers).
%            Header + per-event index are decoded; the 8-byte hit-data
%            block is returned raw (bit-level packing into TOF/detx/dety
%            is not reverse-engineered). Delegate to hitsLoad.

[~, ~, ext] = fileparts(fileName);
isHitsFile = any(strcmpi(ext, {'.HITS', '.hits'}));

if metadata.version == 3 && isHitsFile
    warning('strLoad:hitsFormat', ...
        ['HITS format (v3) detected. strLoad delegates to hitsLoad. ' ...
         'Returned table contains the raw 8-byte hit-data block per ion ' ...
         '(b0..b7), NOT decoded TOF / detx / dety. See ' ...
         'RHIT_imports/HITS_format_analysis.md.']);
    [hits, hdr] = hitsLoad(fileName);
    metadata = hdr;
    return;
elseif metadata.version == 3
    fprintf(['STR file with v3 header detected. Continuing with STR-style ' ...
             'TLV channel parser.\n']);
end

%% Vectorized event extraction
% Strategy: assign each field to its event using cumulative sum of 0x18 tags.
% Tag 0x18 marks the END of each event, so fields between consecutive
% 0x18 markers (plus the 0x18 itself) belong to the same event.

fprintf('Extracting events (vectorized)...\n');

% Event ID for each field: increment at each 0x18 (end-of-event marker)
isEnd = (tagVec == 0x18);
% Event numbering: fields BEFORE the first 0x18 = event 1, etc.
eventId = cumsum(isEnd);
% Shift so that the 0x18 field itself belongs to its own event
eventId = eventId - isEnd + 1;

nEvents = max(eventId);
metadata.nEvents = nEvents;
fprintf('Found %d events\n', nEvents);

% For each channel tag, build a sparse (event → value) mapping.
% Use the LAST occurrence of each tag within an event (closest to 0x18).
channelTags = uint8([0x01 0x02 0x03 0x04 0x21 0x22]);
chData = nan(nEvents, 6);

for ci = 1:6
    mask = (tagVec == channelTags(ci));
    idx = find(mask);
    ev  = eventId(idx);
    v   = double(valVec(idx));
    % Keep last occurrence per event (accumarray with @last)
    % Since idx is sorted, the last value written wins with direct indexing
    chData(ev, ci) = v;
end

% Quality (tag 0x18)
qualData = double(valVec(isEnd));
% Pad if fewer 0x18 than events (shouldn't happen, but be safe)
if numel(qualData) < nEvents
    qualData(end+1:nEvents) = 0;
end

%% Filter: keep only events with at least one channel (discard segment markers)
% Tag 0x18 with value=5000 between pulse counters is a segment marker,
% not a real hit event. Real events have at least one channel tag.
hasAnyChannel = any(~isnan(chData), 2);
sourceEventIdx = uint32(find(hasAnyChannel));
nKept = sum(hasAnyChannel);

chData = chData(hasAnyChannel, :);
qualData = qualData(hasAnyChannel);

nAll6 = sum(~any(isnan(chData), 2));
fprintf('Events: %d with data (of %d 0x18 markers), %d with all 6 channels\n', ...
    nKept, nEvents, nAll6);
metadata.nEvents = nKept;

%% Build output table
hits = table();
hits.ionIdx  = (1:nKept)';
hits.eventIdx = sourceEventIdx;
hits.detxt1  = chData(:,1);     % DL-x end 1 (TDC counts)
hits.detxt2  = chData(:,2);     % DL-x end 2
hits.detyt1  = chData(:,3);     % DL-y end 1
hits.detyt2  = chData(:,4);     % DL-y end 2
hits.detwt1  = chData(:,5);     % DL-w end 1 (45° diagonal)
hits.detwt2  = chData(:,6);     % DL-w end 2
hits.quality = qualData;

hits.Properties.VariableUnits = {'1','1','TDC','TDC','TDC','TDC','TDC','TDC','1'};

fprintf('Loaded %d events from %s\n', nKept, fileName);
end


function metadata = parseHeader(tags, vals)
% Parse the file header from the initial TLV fields

metadata = struct();

% Version (first field, tag 0xa0)
if tags(1) == 0xa0
    metadata.version = bitand(vals(1), 255);
else
    metadata.version = 0;
end

% Detector channel labels (tags 0x27-0x2c)
labels = char(zeros(1, 6));
labelTags = uint8([0x27 0x28 0x29 0x2a 0x2b 0x2c]);
for i = 1:6
    idx = find(tags(1:min(50,numel(tags))) == labelTags(i), 1);
    if ~isempty(idx)
        labels(i) = char(vals(idx));
    end
end
metadata.detectorLabels = labels;

% Thresholds (tags 0x2d-0x2f)
thresh = zeros(3, 1);
for i = 1:3
    idx = find(tags(1:min(50,numel(tags))) == uint8(0x2c + i), 1);
    if ~isempty(idx)
        thresh(i) = vals(idx);
    end
end
metadata.thresholds = thresh;

% Walk corrections (tags 0x33-0x38)
walk = zeros(6, 1);
for i = 1:6
    idx = find(tags(1:min(50,numel(tags))) == uint8(0x32 + i), 1);
    if ~isempty(idx)
        walk(i) = vals(idx);
    end
end
metadata.walkCorrections = walk;

% Initial voltage (tag 0x1b)
idx = find(tags(1:min(100,numel(tags))) == 0x1b, 1);
if ~isempty(idx)
    metadata.voltage = double(vals(idx));
end
end
