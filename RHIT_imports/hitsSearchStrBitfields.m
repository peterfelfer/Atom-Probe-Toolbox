function results = hitsSearchStrBitfields(hitsFile, strFile, options)
% HITSSEARCHSTRBITFIELDS Search HITS payload bitfields against STR timings.
%
% results = hitsSearchStrBitfields(hitsFile, strFile)
%
% Uses a same-run HITS+STR pair. STR supplies decoded delay-line timing
% targets; HITS supplies packed b0..b7 payloads. The search tests contiguous
% unsigned and signed bitfields in both little- and big-endian byte order
% against STR timing-derived targets.
%
% This is a hypothesis generator, not a decoder. Candidate fields must still
% pass cross-file validation and physical constraints.
%
% OPTIONS:
%   eventOffset  - tentative HITS eventIdx -> STR eventIdx offset (default: 0)
%   sampleStride - stride through HITS payload rows (default: 1000)
%   maxRows      - maximum sampled rows after filtering (default: 50000)
%   minBits      - minimum bitfield length (default: 4)
%   maxBits      - maximum bitfield length (default: 24)
%   singleOnly   - restrict to SINGLE, hitInEvent=1 rows (default: true)
%   verbose      - print top candidates (default: true)
%
% See also: hitsAuditStrAlignment
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    strFile (1,:) char
    options.eventOffset (1,1) double = 0
    options.sampleStride (1,1) double = 1000
    options.maxRows (1,1) double = 50000
    options.minBits (1,1) double = 4
    options.maxBits (1,1) double = 24
    options.singleOnly (1,1) logical = true
    options.verbose (1,1) logical = true
end

[hits, ~] = hitsLoad(hitsFile);
[strHits, ~] = strLoad(strFile);

sel = true(height(hits), 1);
if options.singleOnly
    sel = sel & hits.parentType == 'SINGLE' & hits.hitInEvent == 1;
end
sel = sel & ~hits.isDelta;
idx = find(sel);
idx = idx(1:max(1, round(options.sampleStride)):end);
if numel(idx) > options.maxRows
    idx = idx(1:options.maxRows);
end
hitsSample = hits(idx, :);

maxEvent = max([double(hitsSample.eventIdx) + abs(options.eventOffset); double(strHits.eventIdx)]);
strRowByEvent = zeros(maxEvent + 2, 1, 'uint32');
strRowByEvent(double(strHits.eventIdx)) = uint32((1:height(strHits))');
eventQuery = double(hitsSample.eventIdx) + options.eventOffset;
inBounds = eventQuery >= 1 & eventQuery <= numel(strRowByEvent);
rows = zeros(height(hitsSample), 1, 'uint32');
rows(inBounds) = strRowByEvent(eventQuery(inBounds));
valid = rows > 0;

hitsSample = hitsSample(valid, :);
strSample = strHits(double(rows(valid)), :);

targets = buildTargets(strSample);
words = buildWords(hitsSample);

rowsOut = strings(0, 1);
endianOut = strings(0, 1);
signedOut = false(0, 1);
startBitOut = zeros(0, 1);
nBitsOut = zeros(0, 1);
targetOut = strings(0, 1);
rhoOut = zeros(0, 1);
nOut = zeros(0, 1);

for wi = 1:numel(words.names)
    word = words.values{wi};
    for startBit = 0:63
        maxLenHere = min(options.maxBits, 64 - startBit);
        for nBits = options.minBits:maxLenHere
            rawField = extractField(word, startBit, nBits);
            for signedFlag = [false true]
                if signedFlag
                    field = signedField(rawField, nBits);
                else
                    field = rawField;
                end

                for ti = 1:numel(targets.names)
                    y = targets.values{ti};
                    ok = isfinite(field) & isfinite(y);
                    if nnz(ok) < 100
                        continue;
                    end
                    rho = corr(field(ok), y(ok), 'Rows', 'complete');
                    rowsOut(end+1, 1) = words.names(wi); %#ok<AGROW>
                    endianOut(end+1, 1) = words.endian(wi); %#ok<AGROW>
                    signedOut(end+1, 1) = signedFlag; %#ok<AGROW>
                    startBitOut(end+1, 1) = startBit; %#ok<AGROW>
                    nBitsOut(end+1, 1) = nBits; %#ok<AGROW>
                    targetOut(end+1, 1) = targets.names(ti); %#ok<AGROW>
                    rhoOut(end+1, 1) = rho; %#ok<AGROW>
                    nOut(end+1, 1) = nnz(ok); %#ok<AGROW>
                end
            end
        end
    end
end

results = table(rowsOut, endianOut, signedOut, startBitOut, nBitsOut, ...
    targetOut, rhoOut, abs(rhoOut), nOut, ...
    'VariableNames', {'word', 'endian', 'isSigned', 'startBit', ...
                      'nBits', 'target', 'correlation', ...
                      'absCorrelation', 'n'});

sortScore = results.absCorrelation;
sortScore(isnan(sortScore)) = -Inf;
[~, ord] = sort(sortScore, 'descend');
results = results(ord, :);

if options.verbose
    fprintf('\n--- HITS/STR bitfield search ---\n');
    fprintf('Event offset: %+d\n', options.eventOffset);
    fprintf('Sampled aligned rows: %d\n', height(hitsSample));
    disp(results(1:min(30, height(results)), :));
end
end

function targets = buildTargets(strHits)
sx = double(strHits.detxt1) + double(strHits.detxt2);
sy = double(strHits.detyt1) + double(strHits.detyt2);
sw = double(strHits.detwt1) + double(strHits.detwt2);
targets = struct();
targets.names = ["detxRaw","detyRaw","detwRaw","sumX","sumY","sumW","tofRaw","quality"];
targets.values = {
    double(strHits.detxt1) - double(strHits.detxt2)
    double(strHits.detyt1) - double(strHits.detyt2)
    double(strHits.detwt1) - double(strHits.detwt2)
    sx
    sy
    sw
    (sx + sy + sw) / 6
    double(strHits.quality)
};
end

function words = buildWords(hits)
b = cell(1, 8);
for k = 1:8
    b{k} = uint64(hits.(sprintf('b%d', k-1)));
end

wordLE = b{1};
wordBE = b{8};
for k = 2:8
    wordLE = wordLE + bitshift(b{k}, 8*(k-1));
    wordBE = wordBE + bitshift(b{9-k}, 8*(k-1));
end

words = struct();
words.names = ["payloadLE","payloadBE"];
words.endian = ["little","big"];
words.values = {wordLE, wordBE};
end

function field = extractField(word, startBit, nBits)
mask = bitshift(uint64(1), nBits) - 1;
field = double(bitand(bitshift(word, -startBit), mask));
end

function field = signedField(rawField, nBits)
signBit = 2^(nBits - 1);
fullScale = 2^nBits;
field = rawField;
field(rawField >= signBit) = field(rawField >= signBit) - fullScale;
end
