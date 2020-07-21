function pos = posUnDecompose(pos)
% takes a decomposed pos table variable and turns it into a state with
% only ranges allocated
%
% INPUT: 
% pos: decomposed pos file or not decomposed file with ranges allocated
%
% OUTPUT: 
% pos: undecomposed file with only the ranges allocated to the ions

[~, uniqueIonIdx] = unique(pos.ionIdx);

pos = table(pos.ionIdx(uniqueIonIdx),...
    pos.x(uniqueIonIdx),pos.y(uniqueIonIdx),pos.z(uniqueIonIdx),...
    pos.mc(uniqueIonIdx),pos.ion(uniqueIonIdx),pos.chargeState(uniqueIonIdx));

pos.Properties.VariableNames = {'ionIdx','x','y','z','mc','ion','chargeState'};