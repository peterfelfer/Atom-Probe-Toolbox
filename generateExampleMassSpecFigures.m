% generateExampleMassSpecFigures
% Creates mass spectrum figures with ions and ranges for all example
% datasets that have a matching .rrng file.
%
% The figures are saved as .fig files alongside the original data.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Setup
setupToolbox;
load('isotopeTable_naturalAbundances.mat', 'isotopeTable');
load('colorScheme.mat', 'colorScheme');

BIN_WIDTH = 0.01; % Da
SPEC_MODE = 'normalised';
MAX_CHARGE_STATES = [1 2 3];
MAX_ATOM_COUNT = [1 2 3];
USE_RRNG_COLOR = true;

%% Define datasets: {posFile, rrngFile}
baseDir = fileparts(mfilename('fullpath'));
exDir = fullfile(baseDir, 'exampleFiles', 'external_repos');

datasets = {
    fullfile(exDir, '7986279', 'aut_leoben_leitner', 'R21_08680-v02.pos'), ...
    fullfile(exDir, '7986279', 'aut_leoben_leitner', 'R21_08680.rrng');

    fullfile(exDir, '7986279', 'usa_portland_wang', 'R31_06365-v02.pos'), ...
    fullfile(exDir, '7986279', 'usa_portland_wang', 'R31_06365-v02.rrng');

    fullfile(exDir, '7986279', 'usa_denton_smith_apav_gbco', 'R5038_00333-v02.epos'), ...
    fullfile(exDir, '7986279', 'usa_denton_smith_apav_gbco', 'rng_5pj.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m13 Young OES850 --9o7 Mion R31_20551-v03.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m14 Young OES850 --19M R31_20642-v01.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m14 Young OES850 --19M R31_20642-v01.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m16 Young OES850 --8o5 Mion R31_20695-v01.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m17 Old OES850 --17M R31_20647-v01.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m17 Old OES850 --17M R31_20647-v01.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m19 old OES850 --11M R31_20650-v03.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m19 old OES850 --11M R31_20650-v03.rrng');

    fullfile(exDir, '11111719', 'Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.pos'), ...
    fullfile(exDir, '11111719', 'Enamel-2 m20 Old OES850 --8o9M R31_20694-v01.rrng');
};

%% Process each dataset
nDatasets = size(datasets, 1);
results = cell(nDatasets, 1);

for d = 1:nDatasets
    posFile = datasets{d, 1};
    rrngFile = datasets{d, 2};

    [posDir, posName, ~] = fileparts(posFile);
    figFile = fullfile(posDir, [posName '_massSpec.fig']);

    fprintf('\n=== Dataset %d/%d: %s ===\n', d, nDatasets, posName);

    try
        % Check files exist
        if ~isfile(posFile)
            warning('Pos file not found: %s', posFile);
            continue;
        end
        if ~isfile(rrngFile)
            warning('RRNG file not found: %s', rrngFile);
            continue;
        end

        % Load pos data
        fprintf('  Loading %s...\n', posFile);
        pos = posLoad(posFile);
        fprintf('  Loaded %d atoms\n', height(pos));

        % Create mass spectrum
        fprintf('  Creating mass spectrum...\n');
        spec = massSpecPlot(pos.mc, BIN_WIDTH, SPEC_MODE);

        % Read rrng file and extract elements
        rrngText = fileread(rrngFile);
        elements = elementsExtractFromText(rrngText);
        ionList = ionsCreateComplex(elements, MAX_CHARGE_STATES, isotopeTable, MAX_ATOM_COUNT);

        % Split rrng file into lines and apply ranges
        rrngLines = string(splitlines(rrngText));
        rrngLines(rrngLines == "") = [];

        cs = colorScheme; % working copy per dataset
        nRanges = 0;
        for i = 1:length(rrngLines)
            [mcBegin, mcEnd, ionName, ionVolume, ionColor] = rangeInfoFromRRNG(rrngLines(i));
            if ~isempty(mcBegin)
                [~, ~, cs] = rangeAddFromRangeInfo(spec, cs, isotopeTable, ...
                    ionList, mcBegin, mcEnd, ionName, ionVolume, ionColor, USE_RRNG_COLOR);
                nRanges = nRanges + 1;
            end
        end

        fprintf('  Applied %d ranges\n', nRanges);

        % Set figure title
        fig = ancestor(spec, 'figure');
        fig.Name = posName;

        % Save figure
        savefig(fig, figFile);
        fprintf('  Saved: %s\n', figFile);

        results{d} = struct('name', posName, 'status', 'ok', 'nAtoms', height(pos), 'nRanges', nRanges);

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results{d} = struct('name', posName, 'status', 'error', 'message', ME.message);
    end
end

%% Summary
fprintf('\n\n=== SUMMARY ===\n');
for d = 1:nDatasets
    if isempty(results{d})
        continue;
    end
    r = results{d};
    if strcmp(r.status, 'ok')
        fprintf('  OK   : %s (%d atoms, %d ranges)\n', r.name, r.nAtoms, r.nRanges);
    else
        fprintf('  FAIL : %s - %s\n', r.name, r.message);
    end
end
fprintf('Done.\n');
