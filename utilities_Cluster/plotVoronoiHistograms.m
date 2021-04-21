function handle = plotVoronoiHistograms(expVols, randVols, class, bins)

% vols: voronoi volumes of atoms
% class: cluster class of atoms 0: unclustered, 1..n cluster class

if ~exist('bins','var')
    bins = 50;
end


vol = {};
Vmax = 0;

%% creating the individual volume variables
for cl = 0:max(class)
    
    vol{cl+1} = expVols(class == cl);
    
    if median(vol{cl+1}) > Vmax
        Vmax = median(vol{cl+1});
    end
    
end

Vmax = Vmax *3;
x = linspace(0,Vmax,bins);




%% calculating and plotting the histograms
f = figure();


for cl = 1:max(class)+1
    
    maxVol = max(vol{cl});
    maxVol = min(maxVol,Vmax);
    
    vol{cl} = vol{cl}(vol{cl}<Vmax);
    his{cl} = hist(vol{cl},bins)*Vmax/maxVol *fac(cl);
    plot(linspace(0,max(vol{cl}),bins),his{cl},'LineWidth',2);
    hold on;
    
end




randVols = randVols(class == 0);


maxVol = max(randVols);
maxVol = min(maxVol,Vmax);

randVols = randVols(randVols<Vmax);
randHis = hist(randVols,bins(1))*Vmax/maxVol;

plot(linspace(0,max(randVols),bins),randHis,'-r','LineWidth',2);

end