function report = hitsAuditEposAlignment(hitsFile, eposFile, options)
% HITSAUDITEPOSALIGNMENT Audit whether a HITS file can be aligned to EPOS.
%
% report = hitsAuditEposAlignment(hitsFile, eposFile)
%
% Loads the structurally decoded HITS ion-payload stream and a paired EPOS
% table, then reports whether the EPOS file contains the raw fields needed
% for bit-packing reverse engineering (tof, detx, dety, VDC) and whether
% row counts / multiplicity metadata are compatible.
%
% This function does not decode HITS b0..b7. It is a gatekeeper before
% running bit-field searches: if counts or EPOS raw columns fail here, any
% inferred mapping is underconstrained.
%
% INPUTS:
%   hitsFile - path to .HITS file
%   eposFile - path to paired .EPOS file
%
% OPTIONS:
%   'sampleStride' - stride for quick raw-byte correlation checks
%                    (default: 100; set Inf to skip)
%   'verbose'      - print report to command window (default: true)
%
% OUTPUT:
%   report - struct with count checks, EPOS column statistics, and quick
%            raw-byte correlation diagnostics.
%
% See also: hitsLoad, posLoad
%
% (c) Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

arguments
    hitsFile (1,:) char
    eposFile (1,:) char
    options.sampleStride (1,1) double = 100
    options.verbose (1,1) logical = true
end

if ~exist(hitsFile, 'file')
    error('hitsAuditEposAlignment:hitsMissing', ...
        'HITS file does not exist: %s', hitsFile);
end

report = struct();
report.hitsFile = hitsFile;
report.eposFile = eposFile;

[hits, header] = hitsLoad(hitsFile);
report.header = header;
report.nHitsPayloads = height(hits);
report.nSingleIonPayloads = header.nSingleIonPayloads;
report.nMultiIonPayloads = header.nMultiIonPayloads;
report.nMarkerOnlySingleEvents = header.nMarkerOnlySingleEvents;
report.nEmptyMultiEvents = header.nEmptyMultiEvents;
report.nIncompleteMultiEvents = header.nIncompleteMultiEvents;
report.nMultiEventsWithTrailers = header.nMultiEventsWithTrailers;

if ~exist(eposFile, 'file')
    report.eposExists = false;
    report.status = "EPOS_MISSING";
    if options.verbose
        fprintf('\n--- HITS/EPOS alignment audit ---\n');
        fprintf('HITS payload rows: %d\n', report.nHitsPayloads);
        fprintf('EPOS missing: %s\n', eposFile);
    end
    return;
end

report.eposExists = true;
epos = posLoad(eposFile, 'quiet', true);
report.nEposRows = height(epos);
report.rowDiff = report.nEposRows - report.nHitsPayloads;
report.rowDiffPercent = 100 * report.rowDiff / max(1, report.nEposRows);

requiredRaw = ["tof", "detx", "dety", "VDC"];
eposNames = string(epos.Properties.VariableNames);
report.requiredRawColumns = requiredRaw;
report.hasRequiredRawColumns = all(ismember(requiredRaw, eposNames));

report.eposColumns = summarizeNumericColumns(epos);
report.hasNonzeroRawColumns = false;
if report.hasRequiredRawColumns
    nonzero = false(size(requiredRaw));
    for k = 1:numel(requiredRaw)
        col = epos.(char(requiredRaw(k)));
        nonzero(k) = any(isfinite(col) & col ~= 0);
    end
    report.requiredRawColumnsNonzero = table(requiredRaw', nonzero', ...
        'VariableNames', {'column', 'hasNonzeroValues'});
    report.hasNonzeroRawColumns = all(nonzero);
else
    report.requiredRawColumnsNonzero = table();
end

report.multiAgreement = struct();
if ismember("multi", eposNames) && height(epos) == height(hits)
    eposMulti = double(epos.multi);
    hitsMulti = double(hits.nHitsInEvent);
    valid = isfinite(eposMulti);
    report.multiAgreement.nCompared = sum(valid);
    report.multiAgreement.fractionEqual = mean(eposMulti(valid) == hitsMulti(valid));
    report.multiAgreement.correlation = corr(eposMulti(valid), hitsMulti(valid), ...
        'Rows', 'complete');
else
    report.multiAgreement.nCompared = 0;
    report.multiAgreement.fractionEqual = NaN;
    report.multiAgreement.correlation = NaN;
end

report.rawByteCorrelations = table();
if report.hasRequiredRawColumns && report.hasNonzeroRawColumns && ...
        isfinite(options.sampleStride)
    n = min(height(hits), height(epos));
    idx = 1:max(1, round(options.sampleStride)):n;
    report.rawByteCorrelations = quickByteCorrelations(hits(idx, :), epos(idx, :));
end

if report.hasRequiredRawColumns && report.hasNonzeroRawColumns && report.rowDiff == 0
    report.status = "READY_FOR_BIT_SEARCH";
elseif report.hasRequiredRawColumns && report.hasNonzeroRawColumns
    report.status = "RAW_COLUMNS_PRESENT_COUNT_MISMATCH";
elseif report.hasRequiredRawColumns
    report.status = "RAW_COLUMNS_ZERO_OR_EMPTY";
else
    report.status = "RAW_COLUMNS_MISSING";
end

if options.verbose
    printReport(report);
end
end

function stats = summarizeNumericColumns(tbl)
names = string(tbl.Properties.VariableNames);
numericNames = strings(0, 1);
nNonzero = zeros(0, 1);
minVal = zeros(0, 1);
maxVal = zeros(0, 1);

for k = 1:numel(names)
    col = tbl.(char(names(k)));
    if ~isnumeric(col)
        continue;
    end
    finiteCol = col(isfinite(col));
    numericNames(end+1, 1) = names(k); %#ok<AGROW>
    if isempty(finiteCol)
        nNonzero(end+1, 1) = 0; %#ok<AGROW>
        minVal(end+1, 1) = NaN; %#ok<AGROW>
        maxVal(end+1, 1) = NaN; %#ok<AGROW>
    else
        nNonzero(end+1, 1) = nnz(finiteCol ~= 0); %#ok<AGROW>
        minVal(end+1, 1) = min(finiteCol); %#ok<AGROW>
        maxVal(end+1, 1) = max(finiteCol); %#ok<AGROW>
    end
end

stats = table(numericNames, nNonzero, minVal, maxVal, ...
    'VariableNames', {'column', 'nNonzero', 'min', 'max'});
end

function C = quickByteCorrelations(hits, epos)
byteNames = "b" + string(0:7);
targetNames = ["tof", "detx", "dety", "VDC"];
rows = strings(0, 1);
targets = strings(0, 1);
rho = zeros(0, 1);

for bi = 1:numel(byteNames)
    x = double(hits.(char(byteNames(bi))));
    for ti = 1:numel(targetNames)
        y = double(epos.(char(targetNames(ti))));
        rows(end+1, 1) = byteNames(bi); %#ok<AGROW>
        targets(end+1, 1) = targetNames(ti); %#ok<AGROW>
        rho(end+1, 1) = corr(x, y, 'Rows', 'complete'); %#ok<AGROW>
    end
end

C = table(rows, targets, rho, ...
    'VariableNames', {'rawByte', 'target', 'correlation'});
end

function printReport(report)
fprintf('\n--- HITS/EPOS alignment audit ---\n');
fprintf('HITS payload rows: %d (%d SINGLE + %d MULTI)\n', ...
    report.nHitsPayloads, report.nSingleIonPayloads, report.nMultiIonPayloads);
fprintf('Marker-only SINGLE events: %d\n', report.nMarkerOnlySingleEvents);
fprintf('Empty MULTI events: %d\n', report.nEmptyMultiEvents);
fprintf('Incomplete MULTI events: %d\n', report.nIncompleteMultiEvents);
fprintf('MULTI events with trailing records: %d\n', report.nMultiEventsWithTrailers);
fprintf('EPOS rows: %d\n', report.nEposRows);
fprintf('Row diff EPOS - HITS: %+d (%+.4f%%)\n', ...
    report.rowDiff, report.rowDiffPercent);
fprintf('Status: %s\n', report.status);

fprintf('\nRequired raw EPOS columns: %s\n', ...
    strjoin(cellstr(report.requiredRawColumns), ', '));
if ~isempty(report.requiredRawColumnsNonzero)
    disp(report.requiredRawColumnsNonzero);
end

if report.multiAgreement.nCompared > 0
    fprintf('EPOS multi vs HITS nHitsInEvent: equal %.4f, corr %.4f\n', ...
        report.multiAgreement.fractionEqual, report.multiAgreement.correlation);
end

if ~isempty(report.rawByteCorrelations)
    fprintf('\nQuick raw-byte correlations, sorted by |rho|:\n');
    C = report.rawByteCorrelations;
    [~, ord] = sort(abs(C.correlation), 'descend');
    disp(C(ord(1:min(12, height(C))), :));
end
end
