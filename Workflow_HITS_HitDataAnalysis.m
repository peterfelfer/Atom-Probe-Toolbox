%% Workflow: HITS v3 — fixed 8-byte ion-payload block (d1, d2)
% Pulls d1+d2 from every structurally ion-bearing HITS event:
% recovered first-block STATUS/SINGLE/OTHER payloads and MULTI(size=3/5/7/...)
% split into one row per ion.
% The 8 bytes are analysed as a fixed block, conditioned on d1.b3
% (precision/mode byte), d2.b3, and parent event type.

addpath('RHIT_imports');

hitsFile = fullfile('exampleFiles', 'R5121_12823.HITS');
eposCandidates = {
    fullfile('exampleFiles', 'HITS_and_EPOS_links', 'R5121_12823.EPOS')
    fullfile('exampleFiles', 'R5121_12823', 'R5121_12823.EPOS')
    fullfile('exampleFiles_backup', 'external_repos', 'APT_recon', 'R5121_12823.EPOS')
};

[hits, header, events] = hitsLoad(hitsFile);

%% Pull d1+d2 from every ion-bearing event
byteVars = "b" + string(0:7);
H = table2array(hits(:, cellstr(byteVars)));
nHits = height(hits);
fprintf('\nAnalysing d1/d2 from %d ion payloads...\n', nHits);

% Also retain the marker's val16 (a4_lo16) which the doc identified as
% inter-event timing (cumsum monotonic).
markerVal16 = double(hits.markerVal16);

%% Joint distribution of d1.b3 (precision) and d2.b3 (small set)
fprintf('\n--- Joint distribution: d1.b3 (precision) x d2.b3 ---\n');
b3a = double(H(:,4));
b3b = double(H(:,8));

% d1.b3 is multiples of 4; reduce to "precision class" = b3/4
precClass = b3a / 4;
[uPrec, ~, ipc] = unique(precClass);
precCounts = accumarray(ipc, 1);
fprintf('  d1.b3 (raw): %s\n', mat2str(unique(b3a(b3a < 64))'));
fprintf('  Precision class (=d1.b3/4) counts:\n');
for k = 1:min(8, numel(uPrec))
    fprintf('    class %2d (b3=0x%02x): %d events (%.1f%%)\n', ...
        uPrec(k), uPrec(k)*4, precCounts(k), 100*precCounts(k)/nHits);
end

%% For each precision class, byte-distribution of bytes 0..2 of d1
fprintf('\n--- Per-precision-class stats on d1 bytes 0..2 ---\n');
for k = 1:min(6, numel(uPrec))
    cls = uPrec(k);
    sel = precClass == cls;
    if sum(sel) < 1000, continue; end
    fprintf('  class %d (n=%d):\n', cls, sum(sel));
    for p = 1:3
        col = double(H(sel, p));
        fprintf('    d1.b%d  mean=%6.1f std=%6.1f  range=[%d, %d]  entropy=%.2f\n', ...
            p-1, mean(col), std(col), min(col), max(col), shannon(col));
    end
end

%% Same for d2 (bytes 4..6) per class
fprintf('\n--- Per-precision-class stats on d2 bytes 0..2 ---\n');
for k = 1:min(6, numel(uPrec))
    cls = uPrec(k);
    sel = precClass == cls;
    if sum(sel) < 1000, continue; end
    fprintf('  class %d (n=%d):\n', cls, sum(sel));
    for p = 1:3
        col = double(H(sel, p+4));
        fprintf('    d2.b%d  mean=%6.1f std=%6.1f  range=[%d, %d]  entropy=%.2f\n', ...
            p-1, mean(col), std(col), min(col), max(col), shannon(col));
    end
end

%% Cumulative-sum check on inter-event a4_lo16 (should be monotonic)
fprintf('\n--- Inter-event timing (a4_lo16) sanity ---\n');
fprintf('  mean=%.1f  median=%.1f  min=%d  max=%d\n', ...
    mean(markerVal16), median(markerVal16), min(markerVal16), max(markerVal16));
cs = cumsum(markerVal16);
fprintf('  cumsum after %d events = %g  (expected monotonic)\n', ...
    nHits, cs(end));

%% EPOS alignment
eposFile = '';
for k = 1:numel(eposCandidates)
    if exist(eposCandidates{k}, 'file')
        eposFile = eposCandidates{k};
        break;
    end
end

if exist(eposFile, 'file')
    fprintf('\n--- EPOS alignment check ---\n');
    report = hitsAuditEposAlignment(hitsFile, eposFile);
else
    fprintf('\n--- No readable EPOS candidate found — skipping alignment\n');
end

%% Save the H matrix and marker timings to MAT for the brute-force step
outFile = 'hits_R5121_12823_d1d2.mat';
parentType = hits.parentType;
hitInEvent = hits.hitInEvent;
nHitsInEvent = hits.nHitsInEvent;
save(outFile, 'H', 'markerVal16', 'parentType', 'hitInEvent', ...
    'nHitsInEvent', 'header', '-v7.3');
fprintf('\nSaved hit-data block to %s\n', outFile);

% --- helper ---
function H = shannon(col)
    h = histcounts(col, 0:256);
    p = h / sum(h);
    p(p == 0) = [];
    H = -sum(p .* log2(p));
end
