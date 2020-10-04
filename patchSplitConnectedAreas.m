function fvs = patchSplitConnectedAreas(fv)
% patchSplitConnectedAreas splits a patch into individually connected
% patches, e.g. to isolate objects
%
% fvs = patchSplitConnectedAreas(fv);
%
% INPUTS:
% fv        patch consisting of fv.vertices and fv.faces that is to be
%           split
%
% OUTPUTS:
% fvs       vector of patches with fv.faces and fv.vertices with each patch
%           individually fully connected

fvs = splitFV(fv.faces,fv.vertices);