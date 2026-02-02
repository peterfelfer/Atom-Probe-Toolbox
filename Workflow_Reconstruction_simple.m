%[text] # Workflow Reconstruction simple
%[text] 
%[text] Basic reconstruction algorithm to reconstruct your atom probe data from a staight flight path atom probe after: Gault et al., Ultramicroscopy 111 (2011) 448 - 457 
%[text] The Workflow needs a pos file with the detector hit coordinates in nm and the voltage curve as well s some reconstruction Parameters. 
%[text] In the following section these parameters are adjusted by the user and are stored in a reconPara struct array. All of the calculated reconstruction parameters are stored in the table recon Data. The final new calculated x, y, and z coordinates are automatically saved in the pos file. 
%%
%[text] #### Prerequisites
% pos = posToTable(); % load an epos file
reconPara.detEff = 0.36; % detector efficiency %[control:editfield:9de6]{"position":[20,24]}
reconPara.flightLength = 110; % in nm %[control:editfield:0b8f]{"position":[26,29]}
reconPara.evapF = 32.49; % evaporation field in V/nm %[control:editfield:6251]{"position":[19,24]}
reconPara.kf = 3.3; % Field factor %[control:slider:97e5]{"position":[16,19]}
reconPara.ICF = 1.4; % Image compression factor %[control:slider:2c59]{"position":[17,20]}
reconData = table;
  %[control:button:836a]{"position":[1,2]}
%%
%[text] if you don't know average atomic density - you can choose the element that is the main part of your sample
% reconPara.avgDens = 85; % atomic density in atoms/nm³ %[control:editfield:36bb]{"position":[23,25]}
load isotopeTable_naturalAbundances.mat;
element = "Fe"; %[control:editfield:842f]{"position":[11,15]}
reconPara.avgDens = mean(isotopeTable.atomDensity(isotopeTable.element == element));
clearvars element isotopeTable
  %[control:button:706e]{"position":[1,2]}
%%
%[text] #### Constants and variable setup
%[text] Calculation of the detector coordinates in polar form
[reconData.ang, reconData.rad] = cart2pol(pos.detx, pos.dety);
%[text] Calcualting effective detector area:
reconPara.areaDet = (max(reconData.rad))^2 * pi();
%[text] Radius evolution from voltage curve (in nm)
reconData.radEvol = pos.VDC/(reconPara.kf * reconPara.evapF);
  %[control:button:773c]{"position":[1,2]}
%%
%[text] #### Calculate x and y coordinates
%[text] launch angle relative to specimen axis
reconData.thetaP = atan(reconData.rad / reconPara.flightLength); % mm/mm
reconData.theta = reconData.thetaP + asin((reconPara.ICF - 1) * sin(reconData.thetaP)); % ??
%[text] distance from axis and z shift of each hit
[reconData.zP, reconData.distance] = pol2cart(reconData.theta, reconData.radEvol); % nm
%[text] x and y coordinates from the angle on the detector and the distance to
%[text] the specimen axis.
[pos.x, pos.y] = pol2cart(reconData.ang, reconData.distance); % nm
  %[control:button:499f]{"position":[1,2]}
%%
%[text] #### Calculate z coordinate
%[text] the z shift with respect to the top of the cap is Rspec - zP
reconData.zP = reconData.radEvol - reconData.zP;
%[text] accumulative part of z
reconPara.omega = 1 ./ reconPara.avgDens; 
reconData.dz = reconPara.omega * reconPara.flightLength^2 * reconPara.kf^2 * reconPara.evapF^2 / (reconPara.detEff * reconPara.areaDet * reconPara.ICF^2) * pos.VDC.^-2; 
%[text]  in nm^3 \* mm^2 \* V^2/nm^2 / (mm^2 \* V^2)
%[text] wide angle correction
reconData.cumZ = cumsum(double(reconData.dz));
pos.z = reconData.cumZ + reconData.zP;
  %[control:button:5076]{"position":[1,2]}
%[text] #### 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":17.3}
%---
%[control:editfield:9de6]
%   data: {"defaultValue":0.36,"label":"detEff","run":"Section","valueType":"Double"}
%---
%[control:editfield:0b8f]
%   data: {"defaultValue":110,"label":"detEff","run":"Section","valueType":"Double"}
%---
%[control:editfield:6251]
%   data: {"defaultValue":32.49,"label":"detEff","run":"Section","valueType":"Double"}
%---
%[control:slider:97e5]
%   data: {"defaultValue":3.3,"label":"kf","max":7,"min":2,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:slider:2c59]
%   data: {"defaultValue":1.4,"label":"ICF","max":3,"min":1,"run":"Section","runOn":"ValueChanged","step":0.1}
%---
%[control:button:836a]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:36bb]
%   data: {"defaultValue":85,"label":"detEff","run":"Section","valueType":"Double"}
%---
%[control:editfield:842f]
%   data: {"defaultValue":"\"Fe\"","label":"element","run":"Section","valueType":"String"}
%---
%[control:button:706e]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:773c]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:499f]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:5076]
%   data: {"label":"Run","run":"Section"}
%---
