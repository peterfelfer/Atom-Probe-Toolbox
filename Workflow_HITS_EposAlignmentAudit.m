%% Workflow: HITS v3 payload stream vs paired EPOS audit
% Checks whether the structurally decoded HITS ion-payload stream can be
% aligned to a paired EPOS file with raw tof/detx/dety/VDC columns.

addpath('RHIT_imports');

hitsFile = fullfile('exampleFiles', 'R5121_12823.HITS');
eposCandidates = {
    fullfile('exampleFiles', 'HITS_and_EPOS_links', 'R5121_12823.EPOS')
    fullfile('exampleFiles', 'R5121_12823', 'R5121_12823.EPOS')
    fullfile('exampleFiles_backup', 'external_repos', 'APT_recon', 'R5121_12823.EPOS')
};

eposFile = '';
for k = 1:numel(eposCandidates)
    if exist(eposCandidates{k}, 'file')
        eposFile = eposCandidates{k};
        break;
    end
end

if isempty(eposFile)
    fprintf('No readable EPOS candidate found. Checked:\n');
    for k = 1:numel(eposCandidates)
        fprintf('  %s\n', eposCandidates{k});
    end
else
    report = hitsAuditEposAlignment(hitsFile, eposFile);
end
