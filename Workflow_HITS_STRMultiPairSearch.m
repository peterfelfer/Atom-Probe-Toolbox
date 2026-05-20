%% Workflow: HITS v3 multi-hit pair search against same-run STR pairs
% Tests whether hit2-hit1 byte differences within HITS MULTI events map to
% nearby STR detector-channel differences in R5124_00368.

addpath('RHIT_imports');

dataDir = ['/Users/peterfelfer/Dropbox/01_research/0101_shared/' ...
           '010113_data-pool/APT_data/R5124_00368'];

hitsFile = fullfile(dataDir, 'R5124_00368.HITS');
strFile = fullfile(dataDir, 'R5124_00368.STR');

resultsMultiPairs = hitsSearchStrMultiPairs(hitsFile, strFile, ...
    'eventOffset', 0, ...
    'sampleStride', 1, ...
    'maxRows', 50000);

save('hits_str_multi_pair_search_R5124.mat', 'resultsMultiPairs');
