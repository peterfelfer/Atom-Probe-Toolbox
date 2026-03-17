% convertHDF5Examples
% Converts all Globus HDF5 example files to .pos files and creates
% mass spectrum figures with ions and ranges.
%
% Each HDF5 file contains:
%   /mass  - mass-to-charge array (float32)
%   /xyz   - coordinates (3 x N, float32)
%   /ranges/range_N - groups with attributes: element, min_da, max_da, volume, color_rgb
%
% The compact HDF5 ion names (e.g. 'NiO2') are parsed into the spaced
% format used by the toolbox ion list (e.g. 'Ni O2').
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Setup
setupToolbox;
load('isotopeTable_naturalAbundances.mat', 'isotopeTable');
load('colorScheme.mat', 'colorScheme');

BIN_WIDTH = 0.01; % Da
SPEC_MODE = 'normalised';
CHARGE_STATES = [1 2 3];
USE_HDF5_COLOR = true;
MAX_COMPLEXITY_CAP = 4; % cap to avoid combinatorial explosion with many elements

%% Find all HDF5 files
baseDir = fileparts(mfilename('fullpath'));
globusDir = fullfile(baseDir, 'exampleFiles', 'external_repos', 'globus');
files = dir(fullfile(globusDir, '*.hdf5'));
nFiles = length(files);
fprintf('Found %d HDF5 files\n\n', nFiles);

results = cell(nFiles, 1);

%% Process each file
for f = 1:nFiles
    hdf5Path = fullfile(files(f).folder, files(f).name);
    [~, baseName, ~] = fileparts(files(f).name);
    posPath = fullfile(files(f).folder, [baseName '.pos']);
    figPath = fullfile(files(f).folder, [baseName '_massSpec.fig']);

    fprintf('\n=== [%d/%d] %s (%.0f MB) ===\n', f, nFiles, baseName, files(f).bytes/1e6);

    try
        %% 1. Read HDF5 data
        fprintf('  Reading HDF5...\n');
        mass = h5read(hdf5Path, '/mass');
        xyz = h5read(hdf5Path, '/xyz')';
        nAtoms = length(mass);
        fprintf('  %d atoms\n', nAtoms);

        %% 2. Build pos table and export
        pos = table((1:nAtoms)', xyz(:,1), xyz(:,2), xyz(:,3), mass, ...
            'VariableNames', {'ionIdx', 'x', 'y', 'z', 'mc'});
        posExport(pos, posPath);
        fprintf('  Exported .pos\n');

        %% 3. Read ranges from HDF5 and parse ion names
        rangeGrp = h5info(hdf5Path, '/ranges');
        nRanges = length(rangeGrp.Groups);

        hdf5Ranges = struct('element', {}, 'mcBegin', {}, 'mcEnd', {}, ...
            'volume', {}, 'color', {}, 'spacedName', {});

        fileElements = {};
        maxComplexity = 1;

        for r = 1:nRanges
            grp = rangeGrp.Groups(r);
            rng = struct('element', '', 'mcBegin', 0, 'mcEnd', 0, ...
                'volume', 0, 'color', [0 0 0], 'spacedName', '');
            for a = 1:length(grp.Attributes)
                att = grp.Attributes(a);
                switch att.Name
                    case 'element',   rng.element = char(att.Value);
                    case 'min_da',    rng.mcBegin = double(att.Value);
                    case 'max_da',    rng.mcEnd = double(att.Value);
                    case 'volume',    rng.volume = double(att.Value);
                    case 'color_rgb', rng.color = double(att.Value(:)');
                end
            end

            % Parse compact name (e.g. 'NiO2') into elements + atomic numbers
            tokens = regexp(rng.element, '([A-Z][a-z]?)(\d*)', 'tokens');
            atomicNums = [];
            parsedOK = true;
            for t = 1:length(tokens)
                elemSym = tokens{t}{1};
                countStr = tokens{t}{2};
                count = max(1, str2double(countStr));
                if isnan(count), count = 1; end
                try
                    atomNum = symbolConvertAtomicNumber(elemSym);
                    atomicNums = [atomicNums; repmat(atomNum, count, 1)]; %#ok<AGROW>
                    fileElements{end+1} = elemSym; %#ok<SAGROW>
                catch
                    parsedOK = false;
                    break;
                end
            end

            if parsedOK && ~isempty(atomicNums)
                rng.spacedName = ionConvertName(atomicNums, NaN, 'plain', isotopeTable);
                complexity = length(atomicNums);
                if complexity > maxComplexity
                    maxComplexity = complexity;
                end
            else
                rng.spacedName = rng.element;
            end

            hdf5Ranges(r) = rng;
        end

        fileElements = unique(fileElements);
        effectiveComplexity = min(maxComplexity, MAX_COMPLEXITY_CAP);
        if effectiveComplexity < maxComplexity
            fprintf('  %d ranges, %d elements, complexity %d (capped from %d)\n', ...
                nRanges, length(fileElements), effectiveComplexity, maxComplexity);
        else
            fprintf('  %d ranges, %d elements, max complexity %d\n', ...
                nRanges, length(fileElements), maxComplexity);
        end

        %% 4. Create per-file ion list with correct complexity
        % Close all figures first so ionsCreateComplex waitbar doesn't conflict
        close all force;
        ionList = ionsCreateComplex(fileElements, 1:effectiveComplexity, ...
            isotopeTable, CHARGE_STATES);
        % Close waitbar before creating spectrum figure
        close all force;

        %% 5. Create mass spectrum
        fprintf('  Creating mass spectrum...\n');
        spec = massSpecPlot(pos.mc, BIN_WIDTH, SPEC_MODE);

        %% 6. Apply ranges using spaced names
        cs = colorScheme;
        nApplied = 0;
        for r = 1:nRanges
            rng = hdf5Ranges(r);
            [~, ~, cs] = rangeAddFromRangeInfo(spec, cs, isotopeTable, ...
                ionList, rng.mcBegin, rng.mcEnd, rng.spacedName, ...
                rng.volume, rng.color, USE_HDF5_COLOR);
            nApplied = nApplied + 1;
        end
        fprintf('  Applied %d/%d ranges\n', nApplied, nRanges);

        %% 7. Save figure
        fig = ancestor(spec, 'figure');
        fig.Name = baseName;
        savefig(fig, figPath);
        fprintf('  Saved figure\n');

        results{f} = struct('name', baseName, 'status', 'ok', ...
            'nAtoms', nAtoms, 'nRanges', nApplied, 'nTotal', nRanges);

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results{f} = struct('name', baseName, 'status', 'error', ...
            'message', ME.message);
    end

    close all force;
end

%% Summary
fprintf('\n\n=== SUMMARY ===\n');
for f = 1:nFiles
    if isempty(results{f}), continue; end
    r = results{f};
    if strcmp(r.status, 'ok')
        fprintf('  OK   : %-55s %10d atoms  %3d/%3d ranges\n', ...
            r.name, r.nAtoms, r.nRanges, r.nTotal);
    else
        fprintf('  FAIL : %-55s %s\n', r.name, r.message);
    end
end
fprintf('Done.\n');
