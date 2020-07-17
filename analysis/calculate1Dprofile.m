function [profile, binVec] = calculate1Dprofile(axis,bin,solutes,allAtoms)
%calculates a 1D concentration profile along the axis given. allAtoms is
%optional, if not given, profile will be counts. 
solutes = solutes(:,1:3);
numSol = length(solutes(:,1));

origin = axis(1,:);
axis = axis(2,:) - axis(1,:);
solutes = solutes - repmat(origin,[numSol,1]);

axis = axis/norm(axis); %make axis unit length

distSol = sum(solutes.*repmat(axis,[numSol,1]),2);% dot product with scale axis for projected length
minDist = min(distSol);
maxDist = max(distSol);

if exist('allAtoms','var')
    allAtoms = allAtoms(:,1:3);
    numAtoms = length(allAtoms(:,1));
    allAtoms = allAtoms - repmat(origin,[numAtoms,1]);
    distAll = sum(allAtoms.*repmat(axis,[numAtoms,1]),2);% dot product with scale axis for projected length
    minDist = min(distAll);
    maxDist = max(distAll);
end

% calculating bin centers
binVec = linspace(0,10000*bin,10001);
binVec = [fliplr(uminus(binVec(2:end))) binVec];
binVec(binVec<minDist | binVec>maxDist) = [];

histSol = hist(distSol,binVec);

if exist('allAtoms','var')
    histAll = hist(distAll,binVec);
    profile = histSol./histAll;
else
    profile = histSol;
end


