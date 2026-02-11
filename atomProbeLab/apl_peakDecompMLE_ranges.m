
function [ionicRangeComp, failedIons, ionAmounts] = apl_peakDecompMLE_ranges(rangeCounts, overlapTable, p0)
% apl_peakDecompMLE_ranges MLE deconvolution for overlapped ranges

C = overlapTable(:,2:end);

if nargin < 3
    p0 = [];
end

[p, ~] = apl_mleSolver(rangeCounts, C, p0);
ionAmounts = p(:)';
failedIons = ionAmounts < 0 | isnan(ionAmounts);

if any(failedIons)
    ionicRangeComp = zeros(size(C));
    return;
end

% ionicRangeComp: probability of each ion in each range
numRanges = size(C,1);
numIons = size(C,2);
ionicRangeComp = zeros(numRanges, numIons);
for r = 1:numRanges
    denom = sum(C(r,:) .* ionAmounts);
    if denom <= 0
        continue;
    end
    ionicRangeComp(r,:) = (C(r,:) .* ionAmounts) ./ denom;
end
end
