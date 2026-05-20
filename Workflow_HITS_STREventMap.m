%% Workflow: HITS v3 / STR original-event map
% Builds event-level fingerprints for R5124_00368 and reports where the
% HITS and STR streams agree/diverge before bit-level decoding.

addpath('RHIT_imports');

dataDir = ['/Users/peterfelfer/Dropbox/01_research/0101_shared/' ...
           '010113_data-pool/APT_data/R5124_00368'];

hitsFile = fullfile(dataDir, 'R5124_00368.HITS');
strFile = fullfile(dataDir, 'R5124_00368.STR');

eventMap = hitsStrEventMap(hitsFile, strFile);
