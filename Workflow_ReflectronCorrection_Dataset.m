%[text] # **Reflectron correction of a dataset**
%[text] **This workflow shows how to correct the detx and dety coordinates with the correction mesh grid.** 
%[text] Prerequisite
%[text] - intersection Table
%[text] - .epos file - to be corrected \
%%
%[text] ## Correction for 3000 XHR LEAP
% reference data loading, creates a variable 'intersections' with the
% imaging distortions. Only needed once for multiple datasets
load('3000XHR_Leoben_21_14093_intersections.mat')

% atom probe data loading as *.epos
[transPosName, transPosPath] = uigetfile({'*.epos'},'Select an epos file');
transPos = posToTable([transPosPath transPosName]);

% actual data transformation
posIn = posReflectronCorrection(transPos, intersections);

% atom probe data export as *.epos, optional
fileName = [transPosName(1:end-5) '_corrected.epos']
posExport(posIn, fileName);
  %[control:button:5c5e]{"position":[1,2]}
%%
%[text] ## Correction for 4000 XHR LEAP
% reference data loading, creates a variable 'intersections' with the
% imaging distortions. Only needed once for multiple datasets
load('4000XHR_Erlangen_56_4833_intersections.mat')

% atom probe data loading as *.epos
[transPosName, transPosPath] = uigetfile({'*.epos'},'Select an epos file');
transPos = posToTable([transPosPath transPosName]);

% actual data transformation
posIn = posReflectronCorrection(transPos, intersections);

% atom probe data export as *.epos, optional
fileName = [transPosName(1:end-5) '_corrected.epos'] %[output:3c6770a2]
posExport(posIn, fileName);
  %[control:button:792e]{"position":[1,2]}
%%
%[text] ## Correction for 5000 XHR LEAP - Leoben
% reference data loading, creates a variable 'intersections' with the
% imaging distortions. Only needed once for multiple datasets
load('5000XR_Leoben_5124_368_intersections.mat')

% atom probe data loading as *.epos
[transPosName, transPosPath] = uigetfile({'*.epos'},'Select an epos file');
transPos = posToTable([transPosPath transPosName]);

% actual data transformation
posIn = posReflectronCorrection(transPos, intersections);

% atom probe data export as *.epos, optional
fileName = [transPosName(1:end-5) '_corrected.epos'] %[output:9cf3dff3]
posExport(posIn, fileName);
  %[control:button:397e]{"position":[1,2]}
%%
%[text] ## Correction for 5000 XHR LEAP - Oxford
% reference data loading, creates a variable 'intersections' with the
% imaging distortions. Only needed once for multiple datasets
load('5000XR_Oxford_5083_23091_intersections.mat')

% atom probe data loading as *.epos
[transPosName, transPosPath] = uigetfile({'*.epos'},'Select an epos file');
transPos = posToTable([transPosPath transPosName]);

% actual data transformation
posIn = posReflectronCorrection(transPos, intersections);

% atom probe data export as *.epos, optional
fileName = [transPosName(1:end-5) '_corrected.epos'] %[output:26a156f0]
posExport(posIn, fileName);
  %[control:button:7061]{"position":[1,2]}
%%
%[text] ## Correction for 6000 XHR LEAP
% reference data loading, creates a variable 'intersections' with the
% imaging distortions. Only needed once for multiple datasets
load('6000XR_Chalmers_6002_770_intersections.mat')

% atom probe data loading as *.epos
[transPosName, transPosPath] = uigetfile({'*.epos'},'Select an epos file');
transPos = posToTable([transPosPath transPosName]);

% actual data transformation
posIn = posReflectronCorrection(transPos, intersections);

% atom probe data export as *.epos, optional
fileName = [transPosName(1:end-5) '_corrected.epos'] %[output:75c5ea3e]
posExport(posIn, fileName);
  %[control:button:1a98]{"position":[1,2]}
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":26}
%---
%[control:button:5c5e]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:792e]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:397e]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:7061]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:1a98]
%   data: {"label":"Run","run":"Section"}
%---
%[output:3c6770a2]
%   data: {"dataType":"textualVariable","outputData":{"name":"fileName","value":"'R56_06712-v01_corrected.epos'"}}
%---
%[output:9cf3dff3]
%   data: {"dataType":"textualVariable","outputData":{"name":"fileName","value":"'R5117_36845_corrected.epos'"}}
%---
%[output:26a156f0]
%   data: {"dataType":"textualVariable","outputData":{"name":"fileName","value":"'fertig_corrected.epos'"}}
%---
%[output:75c5ea3e]
%   data: {"dataType":"textualVariable","outputData":{"name":"fileName","value":"'fertig_corrected.epos'"}}
%---
