function volumes = ionVolumesFromAtomVolumes(ions,isotopeTable)
% ionVolumesFromAtomVolumes outputs a list of ion volumes based on the ions
% that are contained in the variable ions. This could be a categorical
% array of ions, or a whole pos variable. The ion volumes are simple adds
% of the individual atomic volumes as stored in the isotopeTable. These
% volumes are based on ????