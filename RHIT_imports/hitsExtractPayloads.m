function payloads = hitsExtractPayloads(events, raw)
% HITSEXTRACTPAYLOADS Per-event payload byte ranges from a HITS scan.
%
% payloads = hitsExtractPayloads(events, raw)
%
% Given the events struct and raw byte buffer from hitsScanEvents,
% returns per-event payload offsets and lengths. The payload of event i
% is the bytes immediately following the 0xa4 marker, up to (but not
% including) the next 0xa4 marker.
%
% Indexing convention: payloadOffset is the 0-based byte offset of the
% first payload byte; payloadLength is the byte count. To grab a payload:
%
%   bytesI = raw.bytes(payloads.offset(i)+1 : ...
%                      payloads.offset(i)+payloads.length(i));
%
% OUTPUTS:
%   payloads.offset    - nEvents x 1 (uint32) 0-based byte offset
%   payloads.length    - nEvents x 1 (uint32) byte count
%   payloads.endOffset - nEvents x 1 (uint32) 0-based offset just past payload
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    events (1,1) struct
    raw    (1,1) struct
end

% Marker is at byteOffset .. byteOffset+3.
% Payload runs from byteOffset+4 to byteOffset + nFields*4 - 1.
markerOffset = uint32(events.byteOffset);
nFields      = uint32(events.nFields);

payloads = struct();
payloads.offset    = markerOffset + uint32(4);
payloads.length    = (nFields - uint32(1)) * uint32(4);
payloads.endOffset = payloads.offset + payloads.length;

% Sanity: never exceed the file
endsOk = payloads.endOffset <= uint32(raw.nBytes);
if ~all(endsOk)
    nBad = sum(~endsOk);
    warning('hitsExtractPayloads:overflow', ...
        '%d events have payloads extending past EOF (clipping)', nBad);
    payloads.endOffset(~endsOk) = uint32(raw.nBytes);
    payloads.length(~endsOk) = payloads.endOffset(~endsOk) - payloads.offset(~endsOk);
end
end
