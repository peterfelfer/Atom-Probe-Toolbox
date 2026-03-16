function ranges = rangesFromRNG(filename)
% RANGESFROMRNG Parse an RNG (CAMECA/IVAS text range) file into a range table.
%
% ranges = rangesFromRNG(filename)
%
% INPUT
%   filename - path to an .rng file
%
% OUTPUT
%   ranges - table with columns:
%       mcbegin      - range start (Da)
%       mcend        - range end (Da)
%       ion          - ion name string (e.g. 'Fe', 'Cr O')
%       color        - RGB color [r g b]
%
% The RNG format has two sections:
%   1. Atomic section: single-element ions with flag-based assignment
%   2. Polyatomic extension (optional): complex ion assignments for
%      ranges that have multiple element flags
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

lines = readlines(filename);

%% Parse atomic section
% Line 1: numIons numRanges
header = str2num(char(lines(1))); %#ok<ST2NM>
numIons = header(1);

% Ion definitions: pairs of lines (name, name + r g b)
ionNames = cell(numIons, 1);
ionColors = zeros(numIons, 3);
for i = 1:numIons
    lineIdx = 2 + (i-1)*2;
    ionNames{i} = strtrim(char(lines(lineIdx)));
    parts = strsplit(strtrim(char(lines(lineIdx + 1))));
    ionColors(i,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];
end

% Find separator line (starts with dashes, contains ion names)
sepIdx = 0;
for i = 1:numel(lines)
    line = strtrim(char(lines(i)));
    if startsWith(line, '-') && ~contains(line, 'polyatomic')
        sepIdx = i;
        break;
    end
end

% Parse range lines after separator
mcbegin = [];
mcend = [];
ionOut = {};
colorOut = [];

rangeStart = sepIdx + 1;
for i = rangeStart:numel(lines)
    line = strtrim(char(lines(i)));
    if isempty(line) || startsWith(line, '---')
        break;  % end of atomic section or start of polyatomic
    end
    if ~startsWith(line, '.')
        continue;
    end

    parts = strsplit(line);
    % parts{1} = '.', parts{2} = mcbegin, parts{3} = mcend, parts{4:end} = flags
    mcB = str2double(parts{2});
    mcE = str2double(parts{3});
    flags = cellfun(@str2double, parts(4:end));

    % Determine ion from flags
    if numel(flags) >= numIons
        activeIons = find(flags(1:numIons) > 0);
    else
        activeIons = find(flags > 0);
    end

    if numel(activeIons) == 1
        % Single element
        ionName = ionNames{activeIons};
        col = ionColors(activeIons, :);
    elseif numel(activeIons) > 1
        % Multi-element: build spaced name, will be overwritten by polyatomic section
        elemParts = {};
        for a = 1:numel(activeIons)
            idx = activeIons(a);
            count = flags(idx);
            if count > 1
                elemParts{end+1} = [ionNames{idx} num2str(count)];
            else
                elemParts{end+1} = ionNames{idx};
            end
        end
        ionName = strjoin(elemParts, ' ');
        col = ionColors(activeIons(1), :);
    else
        ionName = 'unknown';
        col = [0.5 0.5 0.5];
    end

    mcbegin(end+1,1) = mcB;
    mcend(end+1,1) = mcE;
    ionOut{end+1,1} = ionName;
    colorOut(end+1,:) = col;
end

%% Parse polyatomic extension (if present)
polyStart = 0;
for i = 1:numel(lines)
    if contains(char(lines(i)), 'polyatomic extension')
        polyStart = i;
        break;
    end
end

if polyStart > 0
    % Header
    polyHeader = str2num(char(lines(polyStart + 1))); %#ok<ST2NM>
    numPolyIons = polyHeader(1);

    % Ion definitions
    polyIonNames = cell(numPolyIons, 1);
    polyIonColors = zeros(numPolyIons, 3);
    for i = 1:numPolyIons
        lineIdx = polyStart + 2 + (i-1)*2;
        compactName = strtrim(char(lines(lineIdx)));
        parts = strsplit(strtrim(char(lines(lineIdx + 1))));
        polyIonColors(i,:) = [str2double(parts{2}), str2double(parts{3}), str2double(parts{4})];

        % Convert compact name (e.g. 'FeO2') to spaced name (e.g. 'Fe O2')
        tokens = regexp(compactName, '([A-Z][a-z]?)(\d*)', 'tokens');
        spacedParts = {};
        for t = 1:numel(tokens)
            elem = tokens{t}{1};
            count = tokens{t}{2};
            if isempty(count)
                spacedParts{end+1} = elem;
            else
                spacedParts{end+1} = [elem count];
            end
        end
        polyIonNames{i} = strjoin(spacedParts, ' ');
    end

    % Find separator
    polySepIdx = 0;
    for i = polyStart+2:numel(lines)
        if startsWith(strtrim(char(lines(i))), '-') && ~contains(char(lines(i)), 'polyatomic')
            polySepIdx = i;
            break;
        end
    end

    % Parse polyatomic range lines and update the main range table
    if polySepIdx > 0
        for i = (polySepIdx + 1):numel(lines)
            line = strtrim(char(lines(i)));
            if isempty(line); continue; end
            if ~startsWith(line, '.'); continue; end

            parts = strsplit(line);
            mcB = str2double(parts{2});
            mcE = str2double(parts{3});
            flags = cellfun(@str2double, parts(4:end));

            activeIons = find(flags(1:min(numPolyIons,numel(flags))) > 0);
            if isempty(activeIons); continue; end

            polyIonName = polyIonNames{activeIons(1)};
            polyCol = polyIonColors(activeIons(1), :);

            % Find and update the matching range in the main table
            for r = 1:numel(mcbegin)
                if abs(mcbegin(r) - mcB) < 0.001 && abs(mcend(r) - mcE) < 0.001
                    ionOut{r} = polyIonName;
                    colorOut(r,:) = polyCol;
                    break;
                end
            end
        end
    end
end

%% Build output table
ion = ionOut;
color = colorOut;
ranges = table(mcbegin, mcend, ion, color);

end
