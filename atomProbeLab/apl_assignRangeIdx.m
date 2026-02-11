
function rangeIdx = apl_assignRangeIdx(mc, rangeTable)
% apl_assignRangeIdx assigns each mc to a range index (0 if none)

rangeIdx = zeros(size(mc));
for r = 1:height(rangeTable)
    mask = mc >= rangeTable.mcbegin(r) & mc <= rangeTable.mcend(r);
    rangeIdx(mask) = r;
end
end
