function hdf5rangeTableAdd(fileName,rangeTable)
% hdf5rangeTableAdd can be used to add range information to an HDF5 file.
% For each range, written are:
%
%/atomProbeTomography/massToChargeRange/name [string] = NULL [] % e.g. Fe2 O3
%/atomProbeTomography/massToChargeRange/ion [string] = NULL [] % e.g. 56Fe2 16O3 ++
%/atomProbeTomography/massToChargeRange/begin [float32] = NULL [Da] % e.g. 32.34
%/atomProbeTomography/massToChargeRange/end [float32] = NULL [Da] % e.g. 33.34
%/atomProbeTomography/massToChargeRange/reconstructionVolume [float32] = NULL [at/nm3] % e.g. 70