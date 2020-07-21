function pos = posUnDecompose(pos)
% posUnDecompose takes a decomposed pos table variable and turns it into a 
% state with only ranges allocated
%
% INPUT 
% pos: decomposed pos file or not decomposed file with ranges allocated
%
% OUTPUT 
% pos: undecomposed file with only the ranges allocated to the ions
%
% NOTE: the output is the same as the result of the function
% posAllocateRange when using the initial pos file (posIn) and the option
% 'raw' as inputs

[~, uniqueIonIdx] = unique(pos.ionIdx);

pos = table(pos.ionIdx(uniqueIonIdx),...
    pos.x(uniqueIonIdx),pos.y(uniqueIonIdx),pos.z(uniqueIonIdx),...
    pos.mc(uniqueIonIdx),pos.ion(uniqueIonIdx),pos.chargeState(uniqueIonIdx));

pos.Properties.VariableNames = {'ionIdx','x','y','z','mc','ion','chargeState'};