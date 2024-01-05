function pos = reflectronImageTransform(pos, intersectionCoordinates)
% this function transforms the real hit locations of an instrument with a
% reflectron into hit locations on a virtual detector, located at the back plane of the
% local electrode. These this plane is xx mm from the aperture of the local
% electrode. They can be used to carry out a 3D reconstruction of the data.

%INPUT
% pos:          table variable with the entries detx and dety as the
%               physical detector hit coordinates in double precision float
% intersectionCoordinates: intersection coordinates that represent known
%               coordinates in the virtual detector plane. Provided for
%               each type of instrument.
%               currently avialable:
%               LEAP 4000X HR / LEAP 4000 HR
%               LEAP 5000XR / LEAP 5000R (to come soon)

%OUTPUT
% pos:          table variable with all entries identical to the input pos
%               variable except for the detx and dety, which have been
%               corrected for the reflectron image distortions.
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

ALPHA = 10; % alpha value for alpha shape the triangulation is based on.

% triangulate the grid in physical detector coordinates
gridTris = alphaTriangulation(...
    alphaShape(intersectionCoordinates.detectorX,intersectionCoordinates.detectorY,ALPHA));

%reminder: X lines are the one closer to parallel to the x axis. Their
%position represents y coordinates!

%find triangle for all hit positions inside triangulation on the detector
triangleNumber = tsearchn([intersectionCoordinates.detectorY, intersectionCoordinates.detectorX],gridTris,...
    [pos.detx, pos.dety]);

% for stray hit coordinates outside detector: find closest triangle of triangulation
isOutside = find(isnan(triangleNumber));
closestVertex = dsearchn([intersectionCoordinates.detectorY, intersectionCoordinates.detectorX],gridTris,...
    [pos.detx, pos.dety]);
for idx = 1:length(isOutside) % MATLAB is finally JIT compiled! :D
    tetIdx = find(sum(gridTris == closestVertex( isOutside(idx) ) ,2));
    triangleNumber(isOutside(idx)) = tetIdx(1);
end

%% performing affine transformation
% defining triangulation objects
detectorTriangulation = triangulation(gridTris,[intersectionCoordinates.detectorY, intersectionCoordinates.detectorX]);
gridTriangulation = triangulation(gridTris,[intersectionCoordinates.gridX, intersectionCoordinates.gridY]);

% carrying out transformations
barycentricCoordinatesPos = cartesianToBarycentric(detectorTriangulation,triangleNumber,[pos.detx, pos.dety]);
coords = barycentricToCartesian(gridTriangulation,triangleNumber,barycentricCoordinatesPos);

% creating corrected pos table
pos.detx = uminus(coords(:,2));
pos.dety = uminus(coords(:,1));