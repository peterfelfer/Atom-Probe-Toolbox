function sdm = spatialDistributionMap(points,binSize,radius)
%calculates a 3d spatial distribution map for the 'points', with a bin size of
% 'binSize'. After Geiser et al. DOI: 10.1017/S1431927607070948

%%unfinished

%% overlay all 3d distributions
 all = repmat(points(:,1:3),1,1,length(points));
 offset = repmat(reshape(points(:,1:3),1,3,length(points)),length(points),1,1);
 sdmPoints = all - offset;
 sdm = reshape(sdmPoints,length(points)^2,3,1);