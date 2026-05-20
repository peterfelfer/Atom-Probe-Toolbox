%[text] # **End-to-End Atom Probe Analysis**
%[text] **Complete guided pipeline from data loading to proxigram**
%[text] This workflow walks through the full atom probe analysis pipeline in a single document. Each section is self-contained and can be run independently (provided previous sections have been executed). For more detail on any step, consult the dedicated workflow referenced in each section.
%[text]
%[text] ### Parameter reference (empirical defaults from common usage)
%[text] | Parameter | Default | Range | Notes |
%[text] | Mass spec bin | 0.01 Da | 0.001-0.1 | 0.04 Da for quick overview |
%[text] | Voxel size | 2 nm | 1-3 nm | 1 = fine detail, 3 = fast |
%[text] | 1D profile bin | 0.5 nm | 0.25-1.0 | 0.25 nm near interfaces |
%[text] | Scatter sample | 100,000 pts | 10k-all | Performance cap |
%[text] | Det. efficiency | per instrument | 0.37-0.80 | Critical for composition |
%[text]
%%
%[text] ## 1. Setup and ground truths
%[text] Run *setupToolbox* at the start of every MATLAB session to add toolbox paths. Then load the isotope table and color scheme that are used throughout the analysis.
setupToolbox;
load isotopeTable_naturalAbundances.mat  % provides: isotopeTable
load colorScheme.mat                      % provides: colorScheme
  %[control:button:0101]{"position":[1,2]}
%%
%[text] ## 2. Load data
%[text] Use *posLoad* to load a .pos, .epos, .apt, or .h5 file. A dialog will open to select the file. For .epos files, detector coordinates (detx, dety) and time-of-flight (tof) are also loaded.
%[text] See the live script ***FirstSteps*** for more details on data loading.
posIn = posLoad;
disp(head(posIn)); % quick look at the first few rows
  %[control:button:0102]{"position":[1,2]}
%%
%[text] ## 3. Mass spectrum
%[text] Create a mass spectrum with a bin width of 0.01 Da (normalised mode). This is the foundation for all subsequent ion identification and ranging.
%[text] See the live script ***FirstSteps*** for more details.
bin = 0.01; % bin width in Da (typical: 0.01 for detailed analysis, 0.04 for quick overview) %[control:slider:0103]{"position":[7,11]}
mode = 'normalised'; %[control:dropdown:0104]{"position":[8,20]}
spec = massSpecPlot(posIn, bin, mode); %[output:01050000]
  %[control:button:0105]{"position":[1,2]}
%%
%[text] ## 4. Ion identification
%[text] Add ion markers to the mass spectrum using *ionAdd*. For each element expected in the sample, specify the element symbol and charge state(s). Common charge states are 1+ and 2+ for most elements, 3+ for some transition metals.
%[text] If you do not know which ions are present, use *ionFind* with an ionList created by *ionsCreateComplex* (see ***FirstSteps***).
%[text] Example for a Ni-based superalloy:
ionAdd(spec, 'Ni', [1 2], isotopeTable, colorScheme);
ionAdd(spec, 'Al', [1 2], isotopeTable, colorScheme);
ionAdd(spec, 'Cr', [1 2], isotopeTable, colorScheme);
ionAdd(spec, 'Co', [1 2], isotopeTable, colorScheme);
  %[control:button:0106]{"position":[1,2]}
%%
%[text] ## 5. Range definition
%[text] Define ranges around the identified ions. The function *rangeAddAll* guides through each unranged ion peak. Set *useMin* to true to interactively define a minimum peak height threshold. Press Enter to skip a peak.
%[text] The *rangeMargin* controls how many Da before and after each peak are shown (typically 3 Da).
rangeMargin = 3; % Da displayed before and after the peak %[control:slider:0107]{"position":[15,16]}
useMin = true; %[control:dropdown:0108]{"position":[10,14]}
rangeAddAll(spec, colorScheme, rangeMargin, useMin);
  %[control:button:0109]{"position":[1,2]}
%%
%[text] ## 6. Range extraction and ion allocation
%[text] Extract the defined ranges from the mass spectrum into a table, then allocate each ion hit to its range. The *'decomposed'* option breaks molecular ions into constituent atoms; *'raw'* keeps them as-is.
rangeTable = rangesExtractFromMassSpec(spec);
option = 'decomposed'; %[control:dropdown:010a]{"position":[10,22]}
pos = posAllocateRange(posIn, rangeTable, option);
  %[control:button:010b]{"position":[1,2]}
%%
%[text] ## 7. Composition
%[text] Calculate the bulk composition of the dataset. The detection efficiency *must match the instrument* used for data acquisition.
%[text] **Instrument detection efficiencies:** 0.37 (LEAP 4000X HR), 0.52 (LEAP 5000 XR), 0.80 (LEAP 5000 XS), 0.50 (EIKOS).
%[text] See the live script ***Concentration*** for advanced composition analysis (background correction, deconvolution).
detEff = 0.52; % detection efficiency — must match the instrument %[control:editfield:010c]{"position":[10,14]}
excludeList = {'unranged'}; % ions to exclude from composition %[control:editfield:010d]{"position":[16,28]}
conc = posCalculateConcentrationSimple(pos, detEff, excludeList, '', 'mode', 'atomic');
disp(conc); %[output:010e0000]
  %[control:button:010e]{"position":[1,2]}
%%
%[text] ## 8. 3D visualisation
%[text] Create a 3D scatter plot of the data. Specify which species to display and the sampling rate (fraction <1 or fixed count >1). A sample of 100,000 points is typically sufficient for an overview.
%[text] See the live script ***3D_Visualisation_of_APT_data*** for turntable animations and advanced plotting.
species = {'Ni', 'Al', 'Cr'}; % species to display %[control:editfield:010f]{"position":[12,30]}
sample = 100000; % number of points to display %[control:editfield:0110]{"position":[10,16]}
[p, ax] = scatterPlotPosData(pos, species, sample, colorScheme); %[output:01110000]
  %[control:button:0111]{"position":[1,2]}
%%
%[text] ## 9. Voxelisation and isosurface
%[text] Create a 3D concentration map by voxelising the data and calculating the local concentration. Then generate an isosurface at a chosen isovalue.
%[text] **IMPORTANT: Axis ordering** — MATLAB's *isosurface* function treats the first dimension as Y and the second as X. Grid vectors must be passed as *gridVec\{2\}, gridVec\{1\}, gridVec\{3\}* (X and Y swapped).
%[text] See the live scripts ***Isosurface*** and ***IsosurfaceToProxigram*** for more detail.
dist = [pos.x pos.y pos.z];
voxBin = [2 2 2]; % voxel size in nm (typical: 1=fine, 2=general, 3=fast) %[control:editfield:0112]{"position":[10,18]}
[binCenters, binEdges] = binVectorsFromDistance(dist, voxBin, 'distance');
gridVec = binCenters;
vox = posToVoxel(pos, gridVec);
isoSpecies = {'Al'}; % species for the isosurface %[control:editfield:0113]{"position":[15,21]}
voxIon = posToVoxel(pos, gridVec, isoSpecies);
concMap = voxIon ./ vox;
isovalue = "8"; % isovalue in at.% (typical exploration range: 1-40) %[control:editfield:0114]{"position":[12,15]}
fv = isosurface(gridVec{2}, gridVec{1}, gridVec{3}, concMap, isovalue);
p = patch(fv, 'FaceColor', [1 1 0]); axis equal; rotate3d on; %[output:01150000]
axisSpatialAptify;
  %[control:button:0115]{"position":[1,2]}
%%
%[text] ## 10. Proxigram
%[text] Calculate a proximity histogram from the isosurface. This shows the concentration profile as a function of distance from the interface. The bin width controls the resolution (typical: 0.25-1.0 nm).
%[text] **NOTE on data types:** If you encounter type errors, ensure pos table columns are double precision by wrapping in *double()*.
%[text] See the live scripts ***Proxigram*** and ***IsosurfaceToProxigram*** for multi-ion proxigrams and plotting.
proxiSpecies = pos(pos.ion == 'Al', :); % select species for proxigram %[control:editfield:0116]{"position":[30,34]}
proxiBin = 0.5; % bin width in nm (typical: 0.25 near interfaces, 0.5 default, 1.0 overview) %[control:editfield:0117]{"position":[12,15]}
[proxi, binVector] = patchCreateProxigram(proxiSpecies, pos, fv, proxiBin); %[output:01180000]
figure;
plot(binVector, proxi * 100, 'LineWidth', 2);
xlabel('distance [nm]');
ylabel('concentration [%]');
  %[control:button:0118]{"position":[1,2]}
%%
%[text] ## 11. Optional: ROI and 1D concentration profile
%[text] Define a region of interest and extract a 1D concentration profile along a direction. See the live scripts ***Defining_ROI_and_cropping_data*** and ***1DConcentrationProfile*** for full details.
%[text] The ROI creation functions (*roiCreateBox*, *roiCreateCylinder*, *roiCreateSphere*) require interactive manipulation in the MATLAB figure window.
% roiCreateBox;
% roiCreateCylinder;
%%
%[text] ## 12. Save results
%[text] Export the ranged mass spectrum figure and the allocated pos table for later use.
% savefig(gcf, 'massSpectrum_ranged.fig');
% posExport(pos, 'allocated_pos.pos');
% writetable(conc, 'composition.csv');
%%
%[text] ##

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":32}
%---
%[control:button:0101]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0102]
%   data: {"label":"Run","run":"Section"}
%---
%[control:slider:0103]
%   data: {"defaultValue":0.01,"label":"bin","max":0.1,"min":0.001,"run":"Nothing","step":0.01}
%---
%[control:dropdown:0104]
%   data: {"defaultValue":"'normalised'","itemLabels":["'count'","'normalised'"],"items":["'count'","'normalised'"],"label":"mode","run":"Nothing"}
%---
%[control:button:0105]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:0106]
%   data: {"label":"Run","run":"Section"}
%---
%[control:slider:0107]
%   data: {"defaultValue":3,"label":"rangeMargin","max":10,"min":1,"run":"Nothing","step":1}
%---
%[control:dropdown:0108]
%   data: {"defaultValue":"true","itemLabels":["true","false"],"items":["true","false"],"label":"useMin","run":"Nothing"}
%---
%[control:button:0109]
%   data: {"label":"Run","run":"Section"}
%---
%[control:dropdown:010a]
%   data: {"defaultValue":"'decomposed'","itemLabels":["'decomposed'","'raw'"],"items":["'decomposed'","'raw'"],"label":"option","run":"Nothing"}
%---
%[control:button:010b]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:010c]
%   data: {"defaultValue":0.52,"label":"detEff","run":"Nothing","valueType":"Double"}
%---
%[control:editfield:010d]
%   data: {"defaultValue":"{'unranged'}","label":"excludeList","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:button:010e]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:010f]
%   data: {"defaultValue":"{'Ni', 'Al', 'Cr'}","label":"species","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:0110]
%   data: {"defaultValue":100000,"label":"sample","run":"Nothing","valueType":"Double"}
%---
%[control:button:0111]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0112]
%   data: {"defaultValue":"[2 2 2]","label":"voxBin","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:0113]
%   data: {"defaultValue":"{'Al'}","label":"isoSpecies","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:editfield:0114]
%   data: {"defaultValue":"\"8\"","label":"isovalue","run":"Nothing","valueType":"String"}
%---
%[control:button:0115]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:0116]
%   data: {"defaultValue":"'Al'","label":"proxiSpecies","run":"Nothing","valueType":"Char"}
%---
%[control:editfield:0117]
%   data: {"defaultValue":0.5,"label":"proxiBin","run":"Nothing","valueType":"Double"}
%---
%[control:button:0118]
%   data: {"label":"Run","run":"Section"}
%---
