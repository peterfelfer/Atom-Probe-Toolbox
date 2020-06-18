function isoList = nucleideList

isoList = [];

for Z = 1:92
    isoWeightsZ = isotopicWeightAll(Z);
    
    isoWeightsZ = [repmat([Z],length(isoWeightsZ(:,1)),1) isoWeightsZ];
    
    isoList = [isoList; isoWeightsZ];
end

isoList = isoList(~~isoList(:,2),:);


end