function isoList = nucleideList
% nucleideList provides a list of all isotopes of all elements with all
% weights and abundances
%
% isoList = nucleideList
%
% OUTPUT
% isoList:  List with all isotopes of all elements with all weights and
%           abundances
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

isoList = [];

for Z = 1:92
    isoWeightsZ = isotopicWeightAll(Z);
    
    isoWeightsZ = [repmat([Z],length(isoWeightsZ(:,1)),1) isoWeightsZ];
    
    isoList = [isoList; isoWeightsZ];
end

isoList = isoList(~~isoList(:,2),:);


end