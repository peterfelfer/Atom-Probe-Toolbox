%% Workflow: HITS v3 payload stream vs same-run STR timing audit
% Uses R5124_00368, which contains both HITS and STR files. The STR file
% carries decoded delay-line TLV channel timings despite a v3 header.

addpath('RHIT_imports');

dataDir = ['/Users/peterfelfer/Dropbox/01_research/0101_shared/' ...
           '010113_data-pool/APT_data/R5124_00368'];

hitsFile = fullfile(dataDir, 'R5124_00368.HITS');
strFile = fullfile(dataDir, 'R5124_00368.STR');

report = hitsAuditStrAlignment(hitsFile, strFile, ...
    'sampleStride', 1000, ...
    'offsetCandidates', (-2000:2000)');
