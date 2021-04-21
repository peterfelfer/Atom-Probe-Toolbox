function Nmin = clusterSizeAnalyse(expClusterSizes,ranClusterSizes)
% clusterSizeAnalyse plotts the cluster sizes as a histogram
% the y axis is the 'concentration' the occurences. This is the invertse of
% the 1d voronoi tessellation of the occurences, or N for N occurences
% of the same value.
%
%
%
% INPUT
% expClusterSizes: experimental cluster sizes
% ranClustersizes: random cluster sizes
%
% OUTPUT
% Nmin 

%% experimental cluster sizes
expClusterSizes = sort(expClusterSizes);
uniqueSizes = unique(expClusterSizes);

for idx = 1:length(uniqueSizes)
    multiplicity(idx) = sum(expClusterSizes == uniqueSizes(idx));
end


dist = uniqueSizes(2:end) - uniqueSizes(1:end-1);
domainSize = [dist(1) (dist(1:end-1)+dist(2:end))/2 dist(end)];

y = multiplicity./domainSize;


%% random cluster sizes
ranClusterSizes = sort(ranClusterSizes);
r_uniqueSizes = unique(ranClusterSizes);

for idx = 1:length(r_uniqueSizes)
    r_multiplicity(idx) = sum(ranClusterSizes == r_uniqueSizes(idx));
end

r_dist = r_uniqueSizes(2:end) - r_uniqueSizes(1:end-1);
r_domainSize = [r_dist(1) (r_dist(1:end-1)+r_dist(2:end))/2 r_dist(end)];

r_y = r_multiplicity./r_domainSize;



%% plot levels at which 50% of all clusters are non-random.

plot_limit = 0.5;

% the sum of both curves, starting at N = 2:
for N = 2: max(max(uniqueSizes),max(r_uniqueSizes))
    
    
    ex = y(uniqueSizes == N);
    if isempty(ex)
        ex = 0;
    end
    
    rand = r_y(r_uniqueSizes == N);
    if isempty(rand)
        rand =0;
    end
    
    
    ex_curve(N) = ex * N;
    rand_curve(N) = rand * N;
    
end

rand_cumulative = cumsum(rand_curve);
ex_cumulative = cumsum(ex_curve);

ratio = ex_curve./rand_curve;
ratio(1) = 1;


%% definfing Nmin as 50% threshold
pct = 1/(1 - plot_limit);

t_lim = find(ratio > pct);

if ~isempty(t_lim)
    plot_limits_N = min(t_lim);
else
    plot_limits_N = 1;
end




max_y = max(max(y(2:end)),max(r_y(2:end)));
%since we were starting at N = 2:
plot_limits_N = plot_limits_N + 1;


Nmin = plot_limits_N;


%% classification errors
% classification error for clustered atoms:
%classErr_clust = (ex_cumulative(Nmin)-rand_cumulative(Nmin))/ex_cumulative(end);
% classification error for unclustered atoms
%classErr_rand = (rand_cumulative(end)-rand_cumulative(Nmin))/ex_cumulative(end);



%% plotting

plot(uniqueSizes(2:end),y(2:end).*uniqueSizes(2:end),'-k','LineWidth',2,'Marker','o');
hold on
plot(r_uniqueSizes(2:end),r_y(2:end).*r_uniqueSizes(2:end),':k','LineWidth',2);

set(gca,'YScale','log');
set(gcf,'Color','w');


xlabel('cluster size [atoms]');
ylabel('frequency [cts]');

leg = "Nmin = " + string(Nmin);
stem(plot_limits_N,max_y,'-r','LineWidth',2,'Marker','none');
legend('cluster size distribution','cluster size distribution (randomized)', leg);










