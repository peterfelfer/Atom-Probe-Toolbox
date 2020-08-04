function pos = posTableFromHDF5(fileName)
% posTableFromHDF5 extracts all per - hit/ion/atom information in an hdf5
% file. 

% check is file is decomposed or not


% check which variables are contained


/atomProbeTomography/detectorEvent/laserPulseCount [uint64] = NULL [1]
/atomProbeTomography/detectorEvent/voltagePulseCount [uint64] = NULL [1]
/atomProbeTomography/detectorEvent/x [float32] = NULL [mm]
/atomProbeTomography/detectorEvent/y [float32] = NULL [mm]
/atomProbeTomography/detectorEvent/timeOfFlight [float32] = NULL [nm]
/atomProbeTomography/detectorEvent/standingVoltage [float32] = NULL [V]
/atomProbeTomography/detectorEvent/voltagePulseAmplitude [float32] = NULL [V]
/atomProbeTomography/detectorEvent/fieldDesorptionMap [uint32] = NULL, NULL [1]

/atomProbeTomography/reconstruction/isDecomposed [bool] = true []
/atomProbeTomography/reconstruction/ion/x [float32] = NULL [nm]
/atomProbeTomography/reconstruction/ion/y [float32] = NULL [nm]
/atomProbeTomography/reconstruction/ion/z [float32] = NULL [nm]
/atomProbeTomography/reconstruction/ion/massToChargeState [float32] = NULL [nm]
/atomProbeTomography/reconstruction/ion/chargeState [uint8] = NULL []
/atomProbeTomography/reconstruction/ion/reconstructionVolume [float32] = NULL []
/atomProbeTomography/reconstruction/ion/fieldEvaporationSequenceIndex [int32] = NULL [nm]
/atomProbeTomography/reconstruction/ion/ion [string] = NULL []

/atomProbeTomography/reconstruction/atom/element [string] = NULL []
/atomProbeTomography/reconstruction/atom/isotope [uint8] = NULL []
/atomProbeTomography/reconstruction/atom/x [float32] = NULL [nm]
/atomProbeTomography/reconstruction/atom/y [float32] = NULL [nm]
/atomProbeTomography/reconstruction/atom/z [float32] = NULL [nm]