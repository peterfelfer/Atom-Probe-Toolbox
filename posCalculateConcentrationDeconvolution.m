function conc = posCalculateConcentrationDeconvolution(pos,detEff,excludeList,ionTable,rangeTable,isotopeTable,volumeName)
% posCalculateConcentrationDeconvolution calvulates the concentration in a
% pos table variable based on the defined ions and ranges. It uses relative
% abundances as stated in the isotopeTable
% WARNING: this function is optimised for readability not execution speed