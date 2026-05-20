%% Workflow: HITS events that map to STR hits but have no standard payload
% Inspects R5124_00368 event indices where STR has channel data but HITS
% does not expose a standard 8-byte ion payload.

addpath('RHIT_imports');

dataDir = ['/Users/peterfelfer/Dropbox/01_research/0101_shared/' ...
           '010113_data-pool/APT_data/R5124_00368'];

hitsFile = fullfile(dataDir, 'R5124_00368.HITS');
strFile = fullfile(dataDir, 'R5124_00368.STR');

report = hitsAnalyzeStrOnlyEvents(hitsFile, strFile, ...
    'maxEventsPerGroup', 50000);
