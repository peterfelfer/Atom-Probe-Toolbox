function posCorrected = posReflectronCorrection(transPos, intersections)
% posReflectronCorrected - corrects the reflectron image distortion in a
% dataset, the function corrects the detx and dety coordinates of a 
%
% INPUT:
% transPos =    pos table that was measured with a reflectron instrument
% intersections = variable for the specific instrument - measured with a
%                   grid electrode
% 
% OUTPUT:
% posCorrected = corrected pos table with the corrected detx and dety
% coordinates - an adjacent reconstruction of the dataset is needed!
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg
%% triangulation of the dataset
alpha = 10; % [mm] alpha value for triangulation
gridTris = alphaTriangulation(alphaShape(intersections.detectorX,intersections.detectorY, alpha));
% alphaTriangulation(shp) returns a triangulation that defines the domain of the alpha shape. 
% Each row in tri specifies a triangle defined by vertex IDs 

%% Create translation mesh allocation of detector hit points

triangleNumber = tsearchn([intersections.detectorY, intersections.detectorX],gridTris,...
    [transPos.detx, transPos.dety]);
isOutside = find(isnan(triangleNumber)); 
closestVertex = dsearchn([intersections.detectorY, intersections.detectorX],gridTris,...
    [transPos.detx, transPos.dety]); 
for idx = 1:length(isOutside) % MATLAB is finally JIT compiled! :D 
    tetIdx = find(sum(gridTris == closestVertex( isOutside(idx) ) ,2)); 
    triangleNumber(isOutside(idx)) = tetIdx(1);
end

%% Performing affine transformation through barycentric coordinate approximation 

% defining triangulation objects
detectorTriangulation = triangulation(gridTris,[intersections.detectorY, intersections.detectorX]); 
gridTriangulation = triangulation(gridTris,[intersections.gridX, intersections.gridY]); 

% carrying out transformations
barycentricCoordinatesPos = cartesianToBarycentric(detectorTriangulation,triangleNumber,[transPos.detx, transPos.dety]); 
coords = barycentricToCartesian(gridTriangulation,triangleNumber,barycentricCoordinatesPos); 

% creating corrected pos table
posCorrected = transPos;
posCorrected.dety = coords(:,2);
posCorrected.detx = coords(:,1);

end

