%[text] # Workflow Peak Background correction
%[text] This workflow shows how to calculate a concentration of a peak with background correction. The background correction is based on a linear fit of the background of a peak. 
%[text] 
%[text] #### One single peak
%[text] Go into your massSpec and zoom to your peak of interest
%[text] The current figure should be your massSpec and the pos file should be in RAW format
%[text] copy and past the code line - it wont work in the live script
[peakData] = peakBGCorrectedCount(pos,rangeTavle,ionName);
%[text] run the function. 
%[text] You have to select 4 points in the massSpec - zwo before the peak and two after the peak. The points mark the background range that is used for the background correction. 
%[text] As output, the figure shows the peak with the background correction and the valuable information. The information is also stored in the struct variable peak data. 
%%
%[text] #### Save your data and calculate the ratio between two elements
%[text] First time you run the code
rangeTable = rangesExtractFromMassSpec(spec);
peakName = rangeTable.rangeName(rangeTable.mcbegin<peakData.loc & rangeTable.mcend>peakData.loc);
peakChargeState = rangeTable.chargeState(rangeTable.mcbegin<peakData.loc & rangeTable.mcend>peakData.loc);
peakConc = table (peakName, peakChargeState, peakData.loc, peakData.counts, peakData.pct, 'VariableNames', ...
    {'peakName', 'peakChargeState', 'peakLoc', 'peakCounts', 'peakPct'});
clearvars peakName peakChargeState
  %[control:button:1039]{"position":[1,2]}

%%
%[text] For each other time 
%[text] run the function on the next peak and then hit the run button at the livescript and the new peak is added to the table
rangeTable = rangesExtractFromMassSpec(spec);
peakName = rangeTable.rangeName(rangeTable.mcbegin<peakData.loc & rangeTable.mcend>peakData.loc);
peakChargeState = rangeTable.chargeState(rangeTable.mcbegin<peakData.loc & rangeTable.mcend>peakData.loc);
peakNew = table (peakName, peakChargeState, peakData.loc, peakData.counts, peakData.pct, 'VariableNames', ...
    {'peakName', 'peakChargeState', 'peakLoc', 'peakCounts', 'peakPct'});
peakConc = [peakConc; peakNew]; %[output:55379602]
clearvars peakName peakChargeState
  %[control:button:68f7]{"position":[1,2]}

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[control:button:1039]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:68f7]
%   data: {"label":"Run","run":"Section"}
%---
%[output:55379602]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'peakConc'."}}
%---
