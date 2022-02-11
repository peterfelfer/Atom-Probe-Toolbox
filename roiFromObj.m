function roi_Obj = roiFromObj(obj)
% ionVolumesFromAtomVolumes outputs a list of ion volumes based on the ions
% that are contained in the variable ions. 
% As input parameters table, categorical or an array of elements is
% possible.
%
% [ionVolumeList, ionsOut] = ionVolumesFromAtomVolumes(ions,isotopeTable);
%
% INPUT
%
% obj = obj that was exportet from e.g. Blender
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
%% loads obj
obj = objToPatch();
%% makes patch from obj
roi_Obj = patch (gca,'Vertices',obj.vertices,'Faces',obj.objects.faces)
%% changes color and name
roi_Obj.FaceColor = [.5 , .5 , .5];
roi_Obj.FaceAlpha = 0.5;
roi_Obj.DisplayName = 'ROI from Obj';
%% coloring of vertices
cols = [...
    1, 0, 0;... x axis coloring
    0, 0, 0;... 
    0, 0, 0;...
    0, 1, 0;... y axis coloring
    0, 0, 1;... z axis coloring
    0, 0, 0;...
    0, 0, 0;...
    0, 0, 0;];
roi_Obj.FaceVertexCData = cols;
roi_Obj.EdgeColor = 'flat';
roi_Obj.LineWidth = 2;