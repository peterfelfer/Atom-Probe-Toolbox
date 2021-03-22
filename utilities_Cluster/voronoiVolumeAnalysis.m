function [numClustered, clusterCutoff, histCounts, experimentalVolumes, randomVolumes] = voronoiVolumeAnalysis(clusterPos, pos,vis,bins, Vmax)
% voronoivolumeAnalysisNEW calculates the Voronoi Volume for a given data set
%
% [numClustered clusterCutoff histCounts experimentalVolumes randomVolumes]
% = voronoiVolumeAnalysisNEW(clusterPos, pos,vis,bins,Vmax);
% [numClustered clusterCutoff histCounts experimentalVolumes randomVolumes]
% = voronoiVolumeAnalysisNEW(clusterPos, pos,vis,bins);
% [numClustered clusterCutoff histCounts experimentalVolumes randomVolumes]
% = voronoiVolumeAnalysisNEW(clusterPos, pos,vis);
%[numClustered clusterCutoff histCounts experimentalVolumes randomVolumes]
% = voronoiVolumeAnalysisNEW(clusterPos, pos);
%
% INPUT
% clusterPos:   a pos file with the atoms that are in the cluster     
% pos:          the parent pos file with all ions of the dataset
% vis:          true logical value for visualisation output
% Vmax:         maximum volume of 
% bins:         the bin size for the histogram default 50
%
% OUTPUT
% numClustered:     number of the cluster with max volume of clusterCutoff
% clusterCutoff:    max volume of the voronoi cell of the clustered atoms
% experimental:     histcounts of the experimental dataset
% random:           histcounts of the random dataset
% experimentalVolumes:  the voronoi volume of each atom of the experimental
%                       dataset
% randomVolumes:    the voronoi volume of each atom of the random dataset
% randpos:          pos file of the random dataset

%% check for a given bin number
if ~exist('bins','var')
    bins = 50;
end
%% calculating the volume of the Voronoi cell of each atom 
vol = vertexVolume(clusterPos);

%calculating the volume of the Voronoi cells of a random sample of atoms
randpos = pos(randsample(height(pos),height(clusterPos)),2:4);
randVol = vertexVolume(randpos);


%% determine maximum volume to plot to
if ~exist('Vmax','var')
    VmaxE = median(vol);
    VmaxR = median(randVol);
    
     Vmax = max(VmaxE,VmaxR) * 3;
end


volHis = histcounts(vol(vol<Vmax),bins); % blue line 
volHisRand = histcounts(randVol(randVol<Vmax),bins); % green line
emr = volHis - volHisRand;   % emr is the difference between the random histogram and the experimental, red line
cs = cumsum(emr);



%% determining the clustering parameters
numClustered = max(cs); % find the max of the histogram
mx = find(cs == numClustered,1); % calculate in which field of the histogram is the max value 

x = linspace(0,Vmax,bins); % create an x vector from 0 to the max value with the bins
clusterCutoff = x(mx);

%% plotting
if exist('vis','var')
    %plotting results
    figure
    plot(x,volHis,'LineWidth',2,'DisplayName','experimental','XDataSource','x','YDataSource','volHis');
    hold on;
    plot(x,volHisRand,'g','LineWidth',2,'DisplayName','random','XDataSource','x','YDataSource','volHisRand');
    
    %experimental - random
    
    hold on;
    plot(x,emr,'r','LineWidth',2,'DisplayName','experimental - random','XDataSource','x','YDataSource','emr');
    
    legend('experimental','random','experimental - random');
    
    
    xlabel('Voronoi volume of atom [nm3]');
    ylabel('frequency [cts]');
    %set(gca,'XLabel','Voronoi volume of atom [nm3]','YLabel','frequency [cts]');
    set(gca,'YGrid','on');
    set(gcf,'Color','w');
    %clstTxt = ['clustering level: ' num2str(numClustered/length(clusterPos)*100,3) '%'];
    %text(0,0,clstTxt);
end

%% set outputs
histCounts = table(volHis', volHisRand');
histCounts.Properties .VariableNames = {'experimental', 'random'};
experimentalVolumes = table(vol');
experimentalVolumes.Properties.VariableNames = {'expVol'};
randomVolumes = table(randVol', randpos.x, randpos.y, randpos.z);
randomVolumes.Properties.VariableNames = {'randVol', 'x', 'y', 'z'};
end

