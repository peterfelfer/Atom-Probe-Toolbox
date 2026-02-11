
function OP = apl_overlapProblem_getRangeCounts(OP, posData, posSubVolIndexList)
% apl_overlapProblem_getRangeCounts computes counts per range

if nargin < 3
    posSubVolIndexList = [];
end

rangeTable = OP.rangeTable;
numRanges = height(rangeTable);

if isnumeric(posData) && isvector(posData) && numel(posData) == numRanges
    OP.rangeCounts = posData(:);
    return;
end

if istable(posData)
    if ~ismember('mc', posData.Properties.VariableNames)
        error('apl_overlapProblem_getRangeCounts:missingMc', 'pos table missing mc column');
    end
    mc = posData.mc;
else
    mc = posData;
end

if ~isempty(posSubVolIndexList)
    mc = mc(posSubVolIndexList);
end

rangeCounts = zeros(numRanges,1);
for r = 1:numRanges
    rangeCounts(r) = sum(mc >= rangeTable.mcbegin(r) & mc <= rangeTable.mcend(r));
end

OP.rangeCounts = rangeCounts;
end
