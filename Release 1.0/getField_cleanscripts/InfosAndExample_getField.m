
% These scripts can calculate Kingham curves for almost all Elements.
% Originaly, these scripts were written by Yao and Li, and (to MM's
% knowledge) used to calulate the Kingham curves in the book "Atom Probe
% Microscopy" (Gault, Moody, Cairney, Ringer).
% The subfolder "Kingham" contains the files from Yao and Li - Baptiste
% Gault has kindly provided the scripts to MM (martin meier).


% In order to calculate run the original scripts to automatically calulate
% (and save to image files) Kingham curves, add the "Kingham" folder to
% your matlab path and run the FinalMain2WithComments.m:

addpath("Kingham\");
FinalMain2WithComments();
% this will save Kingham cuves as tiff files into your current matlab
% working directory.





% Based on these scripts, MM wrote thje Scripts in get_Field, which can
% compute the electric Field when provided with an experimental CSR. The
% computation is based on interpolation between the values as calculated by
% Yao&Li's script and is therfore not as precise as it could possibly, be,
% but hopefully "good enough". User need to be cautious when CSRs are
% partiocularly high or low (say, above 100 or below 0.01). 
% MM really hopes that these scripts calculate sensible results, but cannot
% promise that this is the case. It is a good idea to sanity-check the
% results with another source, such as the Curves in "Atom Probe
% Microscopy"
% To use the script, add the get_Field Folder to the matlab path, and call
% the function CSR2Field with the experimentially mesured CSR and the CSR
% type as parameters. Example:

addpath("get_Field\");
fieldstrength = CSR2Field(1, "Fe++/Fe+")
% this will calculate the electric field for a APT experiment where the
% ratio of Fe++ to Fe+ is 1 (ie, there are as many Fe++ counts as Fe+
% counts).

%More examples:
fieldstrength = CSR2Field(5, "W+++/W++") %Field in a experiment where W+++/W++ = 5
fieldstrength = CSR2Field(0.5, "Pd++/Pd+") %Field in experiment where Pd++/Pd+ = 0.5








