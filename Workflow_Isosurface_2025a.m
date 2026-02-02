%[text] # **Visualisation of Isosurfaces within a reconstructed APT tip**
%[text] **Voxelisation of atom probe data in 3D, calculating concentrations, creating isosurfaces**
%[text] In order to visualise isosurfaces with an APT data set, it is necessary to voxelise the data and calculate the corresponding concentrations. Also, it is essential, that the decomposed pos variable (*pos*) is present in the workspace. These steps are extensively described in the live script ***FirstSteps***.
%[text] In this version, compatible with Matlab 2025a and later, new functions such as linking variables to plots and gif generation are used to create interactivity and visual feedback and documentation without leaving the command line. This enhances workflow documentation.
%[text] ### Create grid vectors
%[text] In the first step, a set of grid vectors in 3D binning is created by the function *binVectorsFromDistance*. As inputs the distance variable and the bin width need to be specified. For the 3D binning the x, y, and z coordinates of the pos file are predestined as the *dist* input. The bin width can be either isotropic (e.g., \[1 1 1\] or anisotropic (e.g., \[1 2 2\]), whereas the positions in *dist* correspond with the position within *bin*. Here, as mode only *'distance'* is possible since the binning needs to be executed in three dimensions.
dist = [pos.x pos.y pos.z]; % distance variable for the binning %[control:editfield:0676]{"position":[8,27]} %[output:48415c1e]
bin = [1 1 1]; % bin width of each voxel in the direction of the defined distances, in nm %[control:editfield:6d64]{"position":[7,14]}
mode = 'distance'; %[control:dropdown:9870]{"position":[8,18]}
[binCenters, binEdges] = binVectorsFromDistance(dist,bin,mode); % creates the bin centers and bin edges of a grid
   %[control:button:64f1]{"position":[2,3]}
%%
%[text] ## Voxelisation of dataset
%[text] In the second step, the dataset is voxelised. In this newer version, this creates a cell variable, where each cell contains a sub pos variable with the same information as the parent pos file. The upside is that we can apply functions to each voxel to extract information. In many cases, this will be a concentration calculation, but it can be all kinds of operations.
vox = posNdBin(pos, [pos.x pos.y pos.z], binEdges);
%[text] This function forms the basis of any subdivision of pos files based on a property. This can be a voxelisation, a 1D concentration profile, an FDM or a proxigram. In case of a voxelisation, the distance to be binned is x y and z.
%[text] 
%%
%[text] ## Concentration calculation
%[text] Now we can calculate the atomic or ionic concentration in each voxel using a 'concentration kernel'. A concentration kernel is a small program that calculates a concentration from a pos variable. We can apply it to the entire dataset, or a subset. You can even write your own, if you want e.g. isotopic information, background subtraction, deconvolution etc. The results will remain usable in the toolbox so long as you follow the format in the supplied concentration kernels. 
%[text] The concentration kernel will be defined as an anonymous function, i.e. a function that is represented by a variable. Some inputs are set at the creation, whereas other inputs are taken each time it is called. For the function details, consult the function reference. Here, on each call, only a pos variable is parsed. The detection efficiency of the instrument and the elements/ions to exclude remain the same.
detectionEfficiency = 0.8;
excludeList = {'unranged'}
ck = @(pos) posCalculateConcentrationSimple(pos, detectionEfficiency, excludeList);
%[text] We can then use the function that applies this concentration kernel to all subvolumes ('bins'). The associated coordinates are the binCenters and 'vol' is the name of the volumes. The latter is not that important for isosurfaces that span entire datasets, but useful if a dataset is e.g. split in ROIs.
conc = binApplyConcentrationKernel(vox,binCenters,ck,'vol');
%[text] Now, the variable conc contains concentration calculations for each of the subvolumes, including atom/ion counts, fractions and variances. To calculate an isosurface, we need to choose a specific value out of these. Typically, this will be a concentration. Note: concentrations are in fractions.
%[text] This will be illustrated on the example of Ti
concTi = cellfun(@(x) x.Ti(2,1),conc);
%[text] This now gives a volume variable with the Ti concentration, that can be used as the basis for an isosurface concentration, or for direct visualisation.
%%
%[text] ## Isosurface generation
%[text] Isosurface generation is built-in to Matlab. We will therefore use the built in function. This function however has one inconsistency: in the generated mesh, the x and y coordinates are swapped relative to the input data. This will be rectified in our script. 
%[text] We will also make use of a new function in Matlab: plot linking. This allows us to link the content of a plot (the isosurface) to be updated any time the content of the variable changes. As a result, we can programatically change the isovalue and observe the effect on the isosurface. We can also run a whole range of isovalues and animate the resulting change in isovalue in e.g. a gif.
% single isovalue
isovalue = 0.05;
fv = isosurface(binCenters{2},binCenters{1},binCenters{3},concTi,isovalue);
fv.vertices = [fv.vertices(:,2), fv.vertices(:,1), fv.vertices(:,3)];
linkdata on;
% Visualize the isosurface
figure;
patch(fv, 'FaceColor', 'red','FaceAlpha',0.6, 'EdgeColor', 'black');
camlight; 
lighting gouraud;
axis equal;
axisSpatialAptify;
rotate3d on;
%%
% animation of isovalue
isoMin = 0.01;
isoMax = 0.1;
isoSteps = 20;
isovalueRange = linspace(isoMin, isoMax, isoSteps); % Define a range of isovalues for animation
for i = 1:length(isovalueRange)
    fv = isosurface(binCenters{2}, binCenters{1}, binCenters{3}, concTi, isovalueRange(i));
    fv.vertices = [fv.vertices(:,2), fv.vertices(:,1), fv.vertices(:,3)];
    pause(0.1); % Pause for a brief moment to create an animation effect
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":7.3}
%---
%[control:editfield:0676]
%   data: {"defaultValue":"[pos.x pos.y pos.z]","label":"Edit field","run":"Section","valueType":"MATLAB code"}
%---
%[control:editfield:6d64]
%   data: {"defaultValue":"[1 1 1]","label":"Edit field","run":"Section","valueType":"MATLAB code"}
%---
%[control:dropdown:9870]
%   data: {"defaultValue":"'distance'","itemLabels":["distance","count"],"items":["'distance'","'count'"],"label":"mode","run":"Section"}
%---
%[control:button:64f1]
%   data: {"label":"Run","run":"Section"}
%---
%[output:48415c1e]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unable to resolve the name 'pos.x'."}}
%---
