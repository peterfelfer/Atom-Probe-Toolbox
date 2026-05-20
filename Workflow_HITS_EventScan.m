%% Workflow: HITS v3 event walker validation
% Runs hitsScanEvents on the paired R5121_12823 HITS file and compares
% counts against the section-3 reference in HITS_format_analysis.md.

addpath('RHIT_imports');

fileName = fullfile('exampleFiles', 'R5121_12823.HITS');

[events, header, raw] = hitsScanEvents(fileName);

%% Compare to reference counts
% From HITS_format_analysis.md, section 3:
%   STATUS = 1,225,266
%   SINGLE = 6,343,221
%   MULTI  =   352,538
%   OTHER  =     6,200
%   TOTAL  = 7,927,225
ref = struct('STATUS', 1225266, 'SINGLE', 6343221, ...
             'MULTI',   352538, 'OTHER',    6200);

obs = struct( ...
    'STATUS', sum(events.type == 'STATUS'), ...
    'SINGLE', sum(events.type == 'SINGLE'), ...
    'MULTI',  sum(events.type == 'MULTI'),  ...
    'OTHER',  sum(events.type == 'OTHER'));

fprintf('\n--- Event count comparison (R5121_12823.HITS) ---\n');
fprintf('  %-7s  %10s  %10s  %10s\n', 'Type', 'Observed', 'Reference', 'Diff');
types = fieldnames(ref);
for k = 1:numel(types)
    t = types{k};
    fprintf('  %-7s  %10d  %10d  %+10d\n', ...
        t, obs.(t), ref.(t), obs.(t) - ref.(t));
end
fprintf('  %-7s  %10d  %10d  %+10d\n', ...
    'TOTAL', numel(events.fieldIdx), ...
    ref.STATUS+ref.SINGLE+ref.MULTI+ref.OTHER, ...
    numel(events.fieldIdx) - (ref.STATUS+ref.SINGLE+ref.MULTI+ref.OTHER));

%% Header sanity check
fprintf('\n--- Header (R5121_12823.HITS) ---\n');
fprintf('  version           = %d\n', header.version);
fprintf('  detectorLabels    = %s\n', header.detectorLabels);
fprintf('  voltage           = %g V\n', header.voltage);
fprintf('  multiHitIndicator = %d\n', header.multiHitIndicator);
fprintf('  thresholds        = [%s]\n', num2str(header.thresholds));
fprintf('  walkCorrections   = [%s]\n', num2str(header.walkCorrections));
fprintf('  timingChannels    = [%s]\n', num2str(header.timingChannels));
fprintf('  nHeaderFields     = %d\n', header.nHeaderFields);

%% Per-type event-size (record gap) sanity check
fprintf('\n--- Per-type event-size distribution (records / 4 bytes) ---\n');
for t = ["STATUS","SINGLE","MULTI","OTHER"]
    sel = events.type == t;
    if ~any(sel), continue; end
    gaps = double(events.nFields(sel));
    [u, ~, ic] = unique(gaps);
    counts = accumarray(ic, 1);
    [~, ord] = sort(counts, 'descend');
    show = min(6, numel(u));
    fprintf('  %-7s n=%d   top sizes:', t, sum(sel));
    for k = 1:show
        fprintf('  %d=%d', u(ord(k)), counts(ord(k)));
    end
    fprintf('\n');
end

%% Multi-hit gap distribution vs reference (section 6)
%   gap=1 -> 171,179    gap=3 -> 163,353
%   gap=5 ->  11,540    gap=7 ->     533
fprintf('\n--- Multi-hit gap distribution ---\n');
multiGaps = double(events.nFields(events.type == 'MULTI'));
fprintf('  %-3s  %10s  %10s\n', 'Gap', 'Observed', 'Reference');
multiRef = struct('g1', 171179, 'g3', 163353, 'g5', 11540, 'g7', 533);
for g = [1 3 5 7]
    obsG = sum(multiGaps == g);
    refG = multiRef.(sprintf('g%d', g));
    fprintf('  %-3d  %10d  %10d\n', g, obsG, refG);
end
