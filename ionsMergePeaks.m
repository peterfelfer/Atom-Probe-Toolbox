function [ionType, abundance, weight] = ionsMergePeaks(ionType, abundance, weight, mergeMargin)
% ionsMergePeaks takes a list of isotopic combinations and merges peaks 
% that are closer than mergeMargin
% 
% [ionType, abundance, weight] = ionsMergePeaks(ionType, abundance, weight, mergeMargin)
% 
% INPUT
% ionType:                      type of the ion for that peaks should 
%                               be merged
% 
% abundance:                    relativ abundance of the peak in 
%                               the mass spectrum
%
% weight:                       weight of the ion
%
% mergeMargin:                  peaks within this range will be merged
% 
% OUTPUT
% [ionType, abundance, weight]: cell array of merged peaks
%                               sorted by molecular weight
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

[weight, sortIdx] = sort(weight);
abundance = abundance(sortIdx);
ionType = ionType(sortIdx);

% determine which peaks are closer than a margin and group them in
% clusters
diffs = weight(2:end) - weight(1:end-1);
isClose = diffs < mergeMargin;

peakCluster = 1;
peakClusterIdx(1) = peakCluster;

numCombos = numel(ionType); % number of isotopic combinations

for i = 1:numCombos-1
    if isClose(i)
        peakClusterIdx(i+1) = peakCluster;
    else
        peakCluster = peakCluster + 1;
        peakClusterIdx(i+1) = peakCluster;
    end
end

% merge individual peaks
for i = 1:peakCluster
    weightTmp(i) = mean(weight(peakClusterIdx == i));
    abundanceTmp(i) = sum(abundance(peakClusterIdx == i));
    % ion type will be taken from the most abundant ion in the cluster
    ionTypePkClust = ionType(peakClusterIdx == i);
    ionTypeTmp{i} = ionTypePkClust{abundance(peakClusterIdx == i) == max(abundance(peakClusterIdx == i))};
end
weight = weightTmp;
abundance = abundanceTmp;
ionType = ionTypeTmp;