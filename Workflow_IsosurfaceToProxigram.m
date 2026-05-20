%[text] # **Isosurface to Proxigram**
%[text] **Unified pipeline from voxelisation through proxigram**
%[text] This workflow combines the isosurface creation and proxigram calculation into a single guided pipeline. It is the most common analysis sequence for studying interfaces and precipitates in APT data.
%[text] **Prerequisites:** A decomposed pos variable must be present in the workspace (see ***FirstSteps*** or ***EndToEndAnalysis***). The *colorScheme* must also be loaded.
%[text]
%[text] **NOTE on data types:** If pos table columns are stored as single precision, the proxigram functions may require double precision. Wrap coordinates in *double()* if you encounter type errors.
%%
%[text] ## 1. Voxelisation
%[text] Create a 3D grid and voxelise the data. The voxel size controls the spatial resolution of the concentration map.
load colorScheme.mat
dist = [pos.x pos.y pos.z];
voxBin = [2 2 2]; % voxel size in nm (typical: 1=fine, 2=general, 3=fast) %[control:editfield:0301]{"position":[10,18]}
mode = 'distance';
[binCenters, binEdges] = binVectorsFromDistance(dist, voxBin, mode);
gridVec = binCenters;
vox = posToVoxel(pos, gridVec); % voxelise all atoms
  %[control:button:0302]{"position":[1,2]}
%%
%[text] ## 2. Concentration map
%[text] Calculate the concentration of the species of interest. The result is a 3D array where each voxel contains the atomic fraction of the selected species.
species = {'Al'}; % species for the isosurface %[control:editfield:0303]{"position":[12,18]}
voxIon = posToVoxel(pos, gridVec, species); % voxelise selected species only
concMap = voxIon ./ vox; % concentration in each voxel (as fraction, not %)
  %[control:button:0304]{"position":[1,2]}
%%
%[text] ## 3. Create isosurface
%[text] Generate an isosurface at a chosen concentration threshold. The isovalue is in at.% (e.g., 8 = 8 at.%). It is common practice to try several values to find the one that best captures the feature of interest.
%[text] **IMPORTANT: Axis ordering** — MATLAB's *isosurface* treats the first dimension as Y and the second as X. Grid vectors must be passed as *gridVec\{2\}, gridVec\{1\}, gridVec\{3\}* (swapped). If your isosurface appears mirrored, check this ordering.
isovalue = "8"; % isovalue in at.% (typical exploration range: 1-40) %[control:editfield:0305]{"position":[12,15]}
fv = isosurface(gridVec{2}, gridVec{1}, gridVec{3}, concMap, isovalue);
figure;
p = patch(fv, 'FaceColor', [0.8 0.8 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
axis equal; rotate3d on; light;
axisSpatialAptify;
title(sprintf('Isosurface at %s at.%%', isovalue)); %[output:03060000]
  %[control:button:0306]{"position":[1,2]}
%%
%[text] ## 4. Single-ion proxigram
%[text] Calculate the proxigram for one species relative to the isosurface. The bin width controls the distance resolution.
proxiSpecies = pos(pos.ion == 'Al', :); % select species %[control:editfield:0307]{"position":[30,34]}
proxiBin = 0.5; % bin width in nm (typical: 0.25 near interfaces, 0.5 default) %[control:editfield:0308]{"position":[12,15]}
[proxi, binVector] = patchCreateProxigram(proxiSpecies, pos, fv, proxiBin); %[output:03090000]
figure;
plot(binVector, proxi * 100, 'LineWidth', 2);
xlabel('distance from interface [nm]');
ylabel('concentration [%]');
title('Single-ion proxigram');
grid on;
  %[control:button:0309]{"position":[1,2]}
%%
%[text] ## 5. Multi-ion proxigram
%[text] Calculate the proxigram for multiple species at once. All concentrations are stored in a single table for easy comparison and plotting.
ionList = {'Al', 'Ni', 'Cr', 'Co'}; % ions to include in the proxigram %[control:editfield:030a]{"position":[12,38]}
proxiBin = 0.5; %[control:editfield:030b]{"position":[12,15]}
proxiAll = table();

for i = 1:length(ionList) %[output:group:030c0000]
    posSpecies = pos(pos.ion == ionList{i}, :);
    [proxi, binVector] = patchCreateProxigram(posSpecies, pos, fv, proxiBin);
    proxi = proxi * 100; % convert to %
    if i == 1
        proxiAll = addvars(proxiAll, binVector', proxi', 'NewVariableNames', {'binVector', ionList{i}});
    else
        proxiAll = addvars(proxiAll, proxi', 'NewVariableNames', ionList(i));
    end
end %[output:group:030c0000]
  %[control:button:030c]{"position":[1,2]}
%%
%[text] ## 6. Plot multi-ion proxigram
%[text] Plot all species in one figure with colors from the color scheme.
figure;
lineWidth = 2;

for k = 1:length(ionList)
    idxElement = find(string(proxiAll.Properties.VariableNames) == ionList{k});
    plot(proxiAll.binVector, table2array(proxiAll(:, idxElement)), ...
        'DisplayName', ionList{k}, 'LineWidth', lineWidth, ...
        'Color', colorScheme.color(colorScheme.ion == ionList{k}, :));
    hold on
end

hold off
legend('Location', 'best');
xlabel('distance from interface [nm]');
ylabel('concentration [%]');
title('Multi-ion proxigram');
grid on; %[output:030d0000]
  %[control:button:030d]{"position":[1,2]}
%%
%[text] ##

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":32}
%---
%[control:editfield:0301]
%   data: {"defaultValue":"[2 2 2]","label":"voxBin","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:button:0302]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0303]
%   data: {"defaultValue":"{'Al'}","label":"species","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:button:0304]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0305]
%   data: {"defaultValue":"\"8\"","label":"isovalue","run":"Nothing","valueType":"String"}
%---
%[control:button:0306]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0307]
%   data: {"defaultValue":"'Al'","label":"proxiSpecies","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:0308]
%   data: {"defaultValue":0.5,"label":"proxiBin","run":"Nothing","valueType":"Double"}
%---
%[control:button:0309]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:030a]
%   data: {"defaultValue":"{'Al', 'Ni', 'Cr', 'Co'}","label":"ionList","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:030b]
%   data: {"defaultValue":0.5,"label":"proxiBin","run":"Nothing","valueType":"Double"}
%---
%[control:button:030c]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:030d]
%   data: {"label":"Run","run":"Section"}
%---
