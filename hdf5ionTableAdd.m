function hdf5ionTableAdd(fileName,ionTable)
% hdf5ionTableAdd can be used to write an ion table into an HDF5 file.
% Written are only the ion types, not the individual isotopic combinations.
% The ions are written in string format in the identifiedIon group:
% /atomProbeTomography/identifiedIon/ion [string] = NULL [] % e.g. Fe2 O3 ++
%
%