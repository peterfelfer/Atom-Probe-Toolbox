
function [F, rangeCenters] = apl_buildFractionMatrix(rangeTable, ions, isotopeTable, minAbundance)
% apl_buildFractionMatrix builds fraction matrix F (ranges x ions)
%
% F(r,i) = fraction of ion i that falls in range r

if nargin < 4 || isempty(minAbundance)
    minAbundance = 0.21; % percent
end

numRanges = height(rangeTable);
numIons = height(ions);

mcbegin = rangeTable.mcbegin;
mcend = rangeTable.mcend;
rangeCenters = (mcbegin + mcend) / 2;

F = zeros(numRanges, numIons);

for i = 1:numIons
    cs = ions.chargeState(i);
    if isnan(cs) || cs == 0
        continue;
    end
    ionTbl = ions.ionTable{i};
    [~, abund, mass] = apl_buildIsotopeList(ionTbl, isotopeTable, minAbundance);
    if isempty(mass)
        continue;
    end
    mz = mass / abs(cs);
    for p = 1:numel(mz)
        idx = find(mz(p) >= mcbegin & mz(p) <= mcend, 1, 'first');
        if ~isempty(idx)
            F(idx, i) = F(idx, i) + abund(p);
        end
    end
    if sum(F(:,i)) > 0
        F(:,i) = F(:,i) ./ sum(F(:,i));
    end
end
end
