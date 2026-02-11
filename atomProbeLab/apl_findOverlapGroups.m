
function [groupRanges, groupIons, rangeToGroup] = apl_findOverlapGroups(F)
% apl_findOverlapGroups finds connected overlap groups in bipartite graph

[numRanges, numIons] = size(F);

rangeToGroup = zeros(numRanges,1);
visitedRange = false(numRanges,1);
visitedIon = false(numIons,1);

groupRanges = {};
groupIons = {};

g = 0;

for r = 1:numRanges
    if visitedRange(r)
        continue;
    end
    if ~any(F(r,:) > 0)
        visitedRange(r) = true;
        continue;
    end
    % start BFS
    queueRanges = r;
    queueIons = [];
    ranges = [];
    ions = [];

    while ~isempty(queueRanges) || ~isempty(queueIons)
        if ~isempty(queueRanges)
            rr = queueRanges(1);
            queueRanges(1) = [];
            if visitedRange(rr)
                continue;
            end
            visitedRange(rr) = true;
            ranges(end+1) = rr; %#ok<AGROW>
            ionIdx = find(F(rr,:) > 0);
            queueIons = [queueIons ionIdx]; %#ok<AGROW>
        else
            ii = queueIons(1);
            queueIons(1) = [];
            if visitedIon(ii)
                continue;
            end
            visitedIon(ii) = true;
            ions(end+1) = ii; %#ok<AGROW>
            rangeIdx = find(F(:,ii) > 0)';
            queueRanges = [queueRanges rangeIdx]; %#ok<AGROW>
        end
    end

    if isempty(ranges)
        continue;
    end

    % only keep groups that include at least one overlapping range
    overlapFlag = any(sum(F(ranges,:) > 0, 2) > 1);
    if ~overlapFlag
        continue;
    end

    g = g + 1;
    groupRanges{g} = unique(ranges, 'stable'); %#ok<AGROW>
    groupIons{g} = unique(ions, 'stable'); %#ok<AGROW>
    rangeToGroup(groupRanges{g}) = g;
end
end
