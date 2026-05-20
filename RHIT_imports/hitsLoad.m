function [hits, header, events] = hitsLoad(fileName, options)
% HITSLOAD Canonical loader for Cameca HITS (v3) detector data files.
%
% [hits, header, events] = hitsLoad(fileName)
% [hits, header, events] = hitsLoad(fileName, 'payloadMode', mode)
%
% HITS files (Cameca AP Suite / IVAS 3.8+) store raw pre-reconstruction
% detector data. The file format is partially decoded: the header and the
% per-event index are fully parsed, but the bit-level packing of the
% 8-byte hit-data block (mapping to TOF / detx / dety) has not been
% reverse-engineered. This loader returns the RAW 8 bytes per ion plus
% the marker timing word, leaving downstream decoding to future work.
%
% INPUTS:
%   fileName - path to .HITS file. If empty, opens a file dialog.
%
% OUTPUTS:
%   hits - table with one row per ion-bearing payload. Any event with at
%          least one complete leading 8-byte payload contributes a row.
%          In the default 'recovered' payload mode, SINGLE events contribute
%          one ion, MULTI events are split into complete payload pairs, and
%          STATUS/OTHER events with at least one complete 8-byte leading
%          block contribute one recovered ion. Parent type is preserved
%          because HITS type bytes are not a reliable no-hit/hit
%          discriminator. Columns:
%            ionIdx       - sequential index (1-based)
%            eventIdx     - parent event index in the events struct
%            byteOffset   - 0-based byte offset of this ion's 8-byte block
%            markerOffset - 0-based byte offset of the parent 0xa4 marker
%            markerVal16  - inter-event timing field (uint16, ~ns)
%            parentType   - STATUS / SINGLE / MULTI / OTHER
%            hitInEvent   - 1-based ion number inside the parent event
%            nHitsInEvent - number of ion payloads in the parent event
%            b0..b7       - the 8 raw hit-data bytes (uint8)
%            eventSize    - record gap to next event (=3 bare hit, =5 hit
%                           + pulse-counter trailer, etc.)
%            hasEventTrailer - true if payload bytes remain after the
%                              complete 8-byte ion blocks
%            isDelta      - true when d1.b3 bit 0 is set (= 0x05 family).
%                           These ~0.1% events are delta-encoded
%                           continuation hits with a different byte layout.
%   header - struct with parsed file header (version, detector channel
%            labels, voltage, thresholds, walk corrections, timing
%            channels, ...).
%   events - struct from hitsScanEvents with the full event index for
%            ALL events (STATUS / SINGLE / MULTI / OTHER), if the caller
%            wants to walk trailers or status records.
%
% LIMITATIONS (HITS v3, as of the current decode):
%   - Bit-level packing of d1+d2 into TOF / detx / dety is unknown.
%     Multi-hit byte-role assay shows bytes 0,1,4,5 vary independently
%     per ion (position-like) and bytes 2,6 are shared/delta within a
%     multi (TOF-like), but no direct decoder works. See
%     RHIT_imports/HITS_format_analysis.md sections 5-7 for the full
%     analysis.
%   - 'recovered' mode is based on R5124 same-run STR evidence: many HITS
%     STATUS events carry one hit-like leading 8-byte block followed by
%     trailer/status records. Use 'allBlocks' for forensic extraction of
%     every complete 8-byte pair, including likely trailers.
%   - MULTI ion events are split into per-ion rows structurally. This
%     recovers the payload stream for alignment, but does not decode the
%     raw 8 bytes into physical detector quantities.
%
% EXAMPLE:
%   [hits, header, events] = hitsLoad('R5121_12823.HITS');
%   fprintf('Detector labels: %s\n', header.detectorLabels);
%   fprintf('%d ion payloads extracted\n', height(hits));
%   % Example: histogram one of the raw bytes
%   histogram(hits.b2, 256);
%
% See also: hitsScanEvents, hitsExtractPayloads, posLoad, strLoad
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    fileName (1,:) char = ''
    options.payloadMode (1,1) string {mustBeMember(options.payloadMode, ...
        ["recovered","standard","allBlocks"])} = "recovered"
end

if isempty(fileName)
    [file, path] = uigetfile( ...
        {'*.HITS;*.hits', 'HITS files (*.HITS)'}, 'Select HITS file');
    if isequal(file, 0)
        hits = table(); header = struct(); events = struct();
        return;
    end
    fileName = fullfile(path, file);
end

%% Walk the file
[events, header, raw] = hitsScanEvents(fileName);

if ~isfield(events, 'fieldIdx') || isempty(events.fieldIdx)
    hits = table();
    return;
end

%% Build the per-ion payload table from leading 8-byte blocks
% Same-run STR comparison shows that many HITS STATUS and marker-like events
% have real STR channel data and carry the same first 8-byte payload
% fingerprint as SINGLE/MULTI events. For SINGLE/STATUS/OTHER, only the
% leading block is treated as an ion by default; later complete pairs are
% often TLV trailers, not additional detector hits.
eventFields = double(events.nFields);
payloadCounts = inferPayloadCounts(events, options.payloadMode);

isSinglePayload = events.type == 'SINGLE' & payloadCounts > 0;
isMarkerOnlySingle = events.type == 'SINGLE' & double(events.nFields) == 1;

isMulti = events.type == 'MULTI';
isMultiPayload = isMulti & payloadCounts > 0;
isEmptyMulti = isMulti & eventFields == 1;
isIncompleteMulti = isMulti & eventFields == 2;
isStatusPayload = events.type == 'STATUS' & payloadCounts > 0;
isOtherPayload = events.type == 'OTHER' & payloadCounts > 0;

nSingleIons = sum(payloadCounts(isSinglePayload));
nMultiIons = sum(payloadCounts(isMultiPayload));
nStatusIons = sum(payloadCounts(isStatusPayload));
nOtherIons = sum(payloadCounts(isOtherPayload));
nHits = sum(payloadCounts);

ionIdx = (1:nHits)';
eventIdx = zeros(nHits, 1, 'uint32');
byteOffset = zeros(nHits, 1, 'uint32');
markerOffset = zeros(nHits, 1, 'uint32');
markerVal = zeros(nHits, 1, 'uint16');
parentType = strings(nHits, 1);
hitInEvent = zeros(nHits, 1, 'uint16');
nHitsInEvent = zeros(nHits, 1, 'uint16');
eventSize = zeros(nHits, 1, 'uint16');
hasEventTrailer = false(nHits, 1);
B = zeros(nHits, 8, 'uint8');

row = 0;
payloadEventIdx = find(payloadCounts > 0);
for k = 1:numel(payloadEventIdx)
    ei = payloadEventIdx(k);
    marker = double(events.byteOffset(ei));
    nThis = payloadCounts(ei);
    typeThis = string(events.type(ei));

    for hi = 1:nThis
        row = row + 1;
        payloadOffset = marker + 4 + (hi - 1) * 8;

        eventIdx(row) = uint32(ei);
        byteOffset(row) = uint32(payloadOffset);
        markerOffset(row) = uint32(marker);
        markerVal(row) = events.markerVal16(ei);
        parentType(row) = typeThis;
        hitInEvent(row) = uint16(hi);
        nHitsInEvent(row) = uint16(nThis);
        eventSize(row) = uint16(events.nFields(ei));
        hasEventTrailer(row) = double(events.nFields(ei)) > (1 + 2*nThis);
        B(row, :) = raw.bytes(payloadOffset+1 : payloadOffset+8);
    end
end

if row ~= nHits
    error('hitsLoad:internalCountMismatch', ...
        'Expected %d payload rows but filled %d rows.', nHits, row);
end

isDelta = bitand(B(:, 4), uint8(1)) ~= 0;
isMultiIon = parentType == "MULTI";

hits = table( ...
    ionIdx, eventIdx, byteOffset, markerOffset, markerVal, ...
    categorical(parentType), hitInEvent, nHitsInEvent, ...
    B(:,1), B(:,2), B(:,3), B(:,4), ...
    B(:,5), B(:,6), B(:,7), B(:,8), ...
    eventSize, hasEventTrailer, isDelta, isMultiIon, ...
    'VariableNames', {'ionIdx', 'eventIdx', 'byteOffset', 'markerOffset', ...
                      'markerVal16', 'parentType', 'hitInEvent', ...
                      'nHitsInEvent', 'b0', 'b1', 'b2', 'b3', ...
                      'b4', 'b5', 'b6', 'b7', ...
                      'eventSize', 'hasEventTrailer', 'isDelta', 'isMulti'});

hits.Properties.VariableUnits = {'1', '1', 'B', 'B', 'ns', '1', '1', '1', ...
    'raw', 'raw', 'raw', 'raw', 'raw', 'raw', 'raw', 'raw', ...
    'records', '1', '1', '1'};

header.format = 'HITSv3';
header.decodeStatus = 'header/event-index/raw-ion-payloads';
header.nEvents = numel(events.fieldIdx);
header.nSingleEvents = sum(events.type == 'SINGLE');
header.nSingleIonPayloads = nSingleIons;
header.nStatusIonPayloads = nStatusIons;
header.nOtherIonPayloads = nOtherIons;
header.payloadMode = char(options.payloadMode);
header.nMarkerOnlySingleEvents = sum(isMarkerOnlySingle);
header.nMultiEvents = sum(isMulti);
header.nMultiIonPayloads = nMultiIons;
header.nEmptyMultiEvents = sum(isEmptyMulti);
header.nIncompleteMultiEvents = sum(isIncompleteMulti);
header.nEventsWithTrailers = sum(payloadCounts > 0 & eventFields > (1 + 2*payloadCounts));
header.nMultiEventsWithTrailers = sum(isMultiPayload & eventFields > (1 + 2*payloadCounts));
header.nIonPayloads = nHits;

if header.nIncompleteMultiEvents > 0
    warning('hitsLoad:incompleteMultiEvents', ...
        '%d MULTI events have only one payload record and were not split.', ...
        header.nIncompleteMultiEvents);
end

fprintf(['hitsLoad: %d ion payloads returned ' ...
         '(%d STATUS + %d SINGLE + %d MULTI + %d OTHER; mode=%s; header: v%d, %s)\n'], ...
    nHits, nStatusIons, nSingleIons, nMultiIons, nOtherIons, ...
    options.payloadMode, header.version, header.detectorLabels);
end

function payloadCounts = inferPayloadCounts(events, payloadMode)
eventFields = double(events.nFields);
completeBlocks = floor(max(0, eventFields - 1) / 2);
payloadCounts = zeros(size(eventFields));

switch payloadMode
    case "standard"
        payloadCounts(events.type == 'SINGLE' & completeBlocks > 0) = 1;
        isMultiPayload = events.type == 'MULTI' & completeBlocks > 0;
        payloadCounts(isMultiPayload) = completeBlocks(isMultiPayload);

    case "recovered"
        isSinglePayload = events.type == 'SINGLE' & completeBlocks > 0;
        isStatusPayload = events.type == 'STATUS' & completeBlocks > 0;
        isOtherPayload = events.type == 'OTHER' & completeBlocks > 0;
        isMultiPayload = events.type == 'MULTI' & completeBlocks > 0;

        payloadCounts(isSinglePayload | isStatusPayload | isOtherPayload) = 1;
        payloadCounts(isMultiPayload) = completeBlocks(isMultiPayload);

    case "allBlocks"
        payloadCounts = completeBlocks;
end
end
