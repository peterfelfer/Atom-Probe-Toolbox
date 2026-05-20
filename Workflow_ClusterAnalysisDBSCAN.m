%[text] # **Cluster Analysis with DBSCAN**
%[text] **Density-based clustering with post-cluster composition analysis**
%[text] This workflow demonstrates how to use the DBSCAN algorithm to identify clusters (e.g., precipitates, solute-enriched regions) in atom probe data. DBSCAN is a density-based approach that does not require a predefined number of clusters. The toolbox implementation supports GPU acceleration for large datasets.
%[text] For the Voronoi-based clustering approach, see the live script ***ClusterDetermination***.
%[text]
%[text] **Prerequisites:** A decomposed pos variable must be in the workspace (see ***FirstSteps*** or ***EndToEndAnalysis***). The *colorScheme* and *isotopeTable* should also be loaded.
%[text]
%[text] **NOTE on data types:** The DBSCAN function accepts both pos tables and Nx3 double arrays. If you encounter type errors with table input, extract coordinates as *double(\[pos.x, pos.y, pos.z\])*.
%%
%[text] ## 1. Select solute atoms
%[text] Choose the solute species whose clustering behaviour you want to analyse. Only atoms of this species are used for the density-based search.
load colorScheme.mat
solute = 'Cu'; % solute species to cluster %[control:editfield:0401]{"position":[11,15]}
clusterPos = pos(pos.ion == solute, :);
fprintf('Selected %d %s atoms for clustering (%.1f%% of ranged atoms)\n', ...
    height(clusterPos), solute, height(clusterPos)/height(pos(pos.ion ~= 'unranged',:))*100); %[output:04020000]
  %[control:button:0402]{"position":[1,2]}
%%
%[text] ## 2. Run DBSCAN
%[text] The two key parameters are:
%[text] - *epsilon*: the neighbourhood radius in nm. Atoms within this distance of each other are considered neighbours.
%[text] - *minPts*: the minimum number of neighbouring atoms required to form a cluster core.
%[text] Typical starting values: epsilon = 0.5-1.0 nm, minPts = 5-20. Lower epsilon finds tighter clusters; higher minPts rejects smaller features.
%[text] The function uses GPU acceleration by default if available.
epsilon = 0.7; % neighbourhood radius in nm %[control:editfield:0403]{"position":[11,14]}
minPts = 10; % minimum points per cluster %[control:editfield:0404]{"position":[10,12]}
[clusterIdx, clusterInfo] = clusterDBSCAN(clusterPos, epsilon, minPts); %[output:04050000]
fprintf('Found %d clusters\n', max(clusterIdx));
fprintf('Clustered atoms: %d (%.1f%% of solute)\n', sum(clusterIdx > 0), sum(clusterIdx > 0)/height(clusterPos)*100);
  %[control:button:0405]{"position":[1,2]}
%%
%[text] ## 3. Visualise clusters
%[text] Plot the clustered atoms in 3D. Noise points (clusterIdx == 0) are excluded.
clusteredAtoms = clusterPos(clusterIdx > 0, :);
clusteredIdx = clusterIdx(clusterIdx > 0);

figure;
scatter3(clusteredAtoms.x, clusteredAtoms.y, clusteredAtoms.z, 3, clusteredIdx, 'filled');
axis equal; rotate3d on;
axisSpatialAptify;
colorbar;
title(sprintf('DBSCAN clusters (\\epsilon = %.2f nm, minPts = %d)', epsilon, minPts)); %[output:04060000]
  %[control:button:0406]{"position":[1,2]}
%%
%[text] ## 4. Cluster size distribution
%[text] Analyse the size distribution of the identified clusters. The Guinier radius provides an effective cluster size measure.
numClusters = max(clusterIdx);
clusterSizes = histcounts(clusterIdx(clusterIdx > 0), numClusters);

figure;
histogram(clusterSizes, 'BinMethod', 'integers');
xlabel('cluster size [atoms]');
ylabel('count');
title('Cluster size distribution');

fprintf('Number of clusters: %d\n', numClusters);
fprintf('Mean cluster size: %.1f atoms\n', mean(clusterSizes));
fprintf('Median cluster size: %.0f atoms\n', median(clusterSizes));
fprintf('Largest cluster: %d atoms\n', max(clusterSizes)); %[output:04070000]
  %[control:button:0407]{"position":[1,2]}
%%
%[text] ## 5. Composition of clustered vs. matrix atoms
%[text] Compare the composition inside and outside the clusters. The detection efficiency must match the instrument used.
%[text] **Instrument detection efficiencies:** 0.37 (LEAP 4000X HR), 0.52 (LEAP 5000 XR), 0.80 (LEAP 5000 XS), 0.50 (EIKOS).
detEff = 0.52; % detection efficiency %[control:editfield:0408]{"position":[10,14]}
excludeList = {'unranged'};

% Identify all atoms in clusters (not just solute, but all species within cluster volumes)
% For a simple approach, use the solute-only clustered atoms:
concClusters = posCalculateConcentrationSimple(clusteredAtoms, detEff, excludeList, 'clusters', 'mode', 'atomic');

% Matrix composition (non-clustered solute atoms for comparison)
matrixAtoms = clusterPos(clusterIdx == 0, :);
concMatrix = posCalculateConcentrationSimple(matrixAtoms, detEff, excludeList, 'matrix', 'mode', 'atomic');

fprintf('\n--- Cluster composition ---\n');
disp(concClusters);
fprintf('\n--- Matrix composition (solute only) ---\n');
disp(concMatrix); %[output:04090000]
  %[control:button:0409]{"position":[1,2]}
%%
%[text] ## 6. Composition of individual clusters
%[text] Loop through each cluster and calculate its composition separately. Results are stored in a table.
clusterCompositions = table();

for i = 1:numClusters
    atomsInCluster = clusteredAtoms(clusteredIdx == i, :);
    if height(atomsInCluster) >= minPts % only analyse clusters above the threshold
        concI = posCalculateConcentrationSimple(atomsInCluster, detEff, excludeList, ...
            sprintf('cluster_%d', i), 'mode', 'atomic');
        concI = concI(concI.format == 'concentration', :);
        concI.clusterID = repmat(i, height(concI), 1);
        concI.clusterSize = repmat(height(atomsInCluster), height(concI), 1);
        clusterCompositions = [clusterCompositions; concI];
    end
end

disp(clusterCompositions); %[output:040a0000]
  %[control:button:040a]{"position":[1,2]}
%%
%[text] ## 7. Alpha hull of clusters (optional)
%[text] Generate an alpha hull around each cluster for visualisation or volume estimation. The *patchCreateSampledAlphaHull* function creates a surface mesh.
%[text] This step can be slow for many clusters. Adjust the cluster index to visualise individual clusters.
clusterToShow = 1; % index of the cluster to visualise %[control:editfield:040b]{"position":[17,18]}
atomsInCluster = clusteredAtoms(clusteredIdx == clusterToShow, :);

if height(atomsInCluster) >= 20 % need enough atoms for a hull
    alpha = 5; % concavity radius
    hull = patchCreateSampledAlphaHull(atomsInCluster, alpha);
    figure;
    patch(hull, 'FaceColor', [0.8 0.2 0.2], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    hold on;
    scatter3(atomsInCluster.x, atomsInCluster.y, atomsInCluster.z, 5, 'filled');
    axis equal; rotate3d on; light;
    title(sprintf('Cluster %d — %d atoms', clusterToShow, height(atomsInCluster))); %[output:040c0000]
else
    fprintf('Cluster %d has only %d atoms — too few for an alpha hull.\n', clusterToShow, height(atomsInCluster));
end
  %[control:button:040c]{"position":[1,2]}
%%
%[text] ##

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":32}
%---
%[control:editfield:0401]
%   data: {"defaultValue":"'Cu'","label":"solute","run":"Nothing","valueType":"Char"}
%---
%[control:button:0402]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0403]
%   data: {"defaultValue":0.7,"label":"epsilon","run":"Nothing","valueType":"Double"}
%---
%[control:editfield:0404]
%   data: {"defaultValue":10,"label":"minPts","run":"Nothing","valueType":"Double"}
%---
%[control:button:0405]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0406]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0407]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0408]
%   data: {"defaultValue":0.52,"label":"detEff","run":"Nothing","valueType":"Double"}
%---
%[control:button:0409]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:040a]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:040b]
%   data: {"defaultValue":1,"label":"clusterToShow","run":"Nothing","valueType":"Double"}
%---
%[control:button:040c]
%   data: {"label":"Run","run":"Section"}
%---
