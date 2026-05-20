%% Workflow: are STATUS-event payloads structurally real ions?
% The recovered-mode hitsLoad assumes that a STATUS event with a complete
% leading 8-byte block at offset+4 can carry an ion. This raised the HITS
% count from 6.24M to 7.54M on R5121, closing the EPOS gap from 16% to 1.5%.
% But before trusting it, we need to verify that bytes pulled from
% STATUS-payload events have the same byte-distribution fingerprint as those
% from SINGLE-payload events.
%
% Real-hit fingerprint (from earlier work, n=6.24M SINGLEs):
%   b0, b1, b4: ~uniform 8-bit (entropy ~8)
%   b5: range [13,242], entropy ~5.2
%   b2, b6: constrained / bimodal
%   b3: multiples of 4 with bit 0 as anomaly flag
%
% If STATUS payloads match this fingerprint, they're real ions.
% If they look totally different (e.g. b3 dominated by 0x05/0x0b/0x49 — known
% TLV tags), they're bookkeeping bytes and recovered mode is over-extracting.

addpath('RHIT_imports');

[hits, header] = hitsLoad('exampleFiles/R5121_12823.HITS');

isSingle = hits.parentType == 'SINGLE';
isStatus = hits.parentType == 'STATUS';
isOther  = hits.parentType == 'OTHER';
isMulti  = hits.parentType == 'MULTI';

groups = {'SINGLE', isSingle; 'STATUS', isStatus; ...
          'MULTI',  isMulti;  'OTHER',  isOther};

byteNames = "b" + string(0:7);

fprintf('\n--- Byte fingerprint by parentType ---\n');
fprintf('%-8s %-8s %8s %8s %8s %8s %10s\n', ...
    'group', 'byte', 'mean', 'std', 'min', 'max', 'entropy');

for gi = 1:size(groups, 1)
    name = groups{gi, 1};
    sel  = groups{gi, 2};
    n    = sum(sel);
    if n < 1000, continue; end
    fprintf('  %s (n=%d)\n', name, n);
    for p = 1:8
        col = double(hits.(char(byteNames(p)))(sel));
        h = histcounts(col, 0:256);
        prob = h / sum(h);
        prob(prob == 0) = [];
        H = -sum(prob .* log2(prob));
        fprintf('  %-8s %-8s %8.1f %8.1f %8d %8d %10.3f\n', ...
            '', char(byteNames(p)), mean(col), std(col), ...
            min(col), max(col), H);
    end
end

%% b3 top values per group (since b3 is the most diagnostic)
fprintf('\n--- d1.b3 top values per group ---\n');
for gi = 1:size(groups, 1)
    name = groups{gi, 1};
    sel  = groups{gi, 2};
    n    = sum(sel);
    if n < 1000, continue; end
    b3 = double(hits.b3(sel));
    [u, ~, ic] = unique(b3);
    cnt = accumarray(ic, 1);
    [~, ord] = sort(cnt, 'descend');
    show = min(8, numel(u));
    fprintf('  %-7s n=%d  top:', name, n);
    topMass = 0;
    for k = 1:show
        v = u(ord(k));
        c = cnt(ord(k));
        fprintf('  0x%02x:%d', v, c);
        topMass = topMass + c;
    end
    fprintf('   (top-%d mass=%.2f)\n', show, topMass / n);
end

%% Compare: are the same b3 values dominant across groups?
% If SINGLE dominantly uses 0x00, 0x04, 0x08, 0x0c, 0x10, 0x14, 0x18, 0x1c
% and STATUS dominantly uses 0x05, 0x0b, 0x49, 0x1e — STATUS payloads are
% TLV records, NOT ion data, and we are over-extracting.

fprintf('\n--- Tag-family check: how many "ion-flag" b3 values vs "TLV-tag" values per group ---\n');
ionLikeB3 = uint8([0x00 0x04 0x08 0x0c 0x10 0x14 0x18 0x1c 0x20 0x24 0x28 0x2c 0x30 0x34 0x38 0x3c]);
tlvLikeB3 = uint8([0x05 0x0b 0x09 0x49 0x1e 0x20 0x47 0x4c 0x4d 0x4e 0x4f 0x50]);

for gi = 1:size(groups, 1)
    name = groups{gi, 1};
    sel  = groups{gi, 2};
    n    = sum(sel);
    if n < 1000, continue; end
    b3 = hits.b3(sel);
    nIon = sum(ismember(b3, ionLikeB3));
    nTlv = sum(ismember(b3, tlvLikeB3));
    fprintf('  %-7s ion-like b3: %.2f%%   tlv-like b3: %.2f%%   other: %.2f%%\n', ...
        name, 100*nIon/n, 100*nTlv/n, 100*(n - nIon - nTlv)/n);
end
