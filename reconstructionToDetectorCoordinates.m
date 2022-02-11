function detectorCoordinates = reconstructionToDetectorCoordinates(reconstructionCoordinates,pos,detx,dety,ionIdx)
%translates the reconstructionCoordinates into detectorCoordinates 
% for co - reconstruction. 
%
% detectorCoordinates = reconstructionToDetectorCoordinates(reconstructionCoordinates,pos);
% detectorCoordinates = reconstructionToDetectorCoordinates(reconstructionCoordinates,pos,detx,dety,ionIdx);
%
% INPUT
% reconstructionCoordinates:  coordinates of the points that are to be
%               translated into detector coordinates as an 3xN array
%
% pos:          pos variable containing the reconstructed corrdinates of
%               the atom probe dataset. if a pos table containing both
%               detector data (pos.detx, pos.dety) and reconstruction data 
%               (pos.x, pos.y, pos.z) is provided, no
%               further input is required. If pos is provided as a simple
%               coordinate list, the other coordinates can be provided
%               sperately
%
% detx/dety:    detector x/y coordinates
%
% ionIdx:       hit index of the ions in the pos/det variable if
%               downsampling is used
%
% OUTPUT
% detectorCoordinates: coordinates of the points in
%               reconstructionCoordinates in detector space.
%
% IMPORTANT: downsampling should be used for the pos and detector
% coordinated, as large datasets lead to excessive memory use and compute
% times. Downsampled dataset sizes of 10^5 - 10^6 are reasonable on a
% machine with 16GB RAM. Use e.g.: 
%
% numSampledAtoms = 1E5;
% sample = sort(randsample(height(pos),numSampledAtoms));
% fvDetector.vertices = object2detectorCoords(fv.vertices,pos(sample,:));
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


% convert everything to double precision for Qhull
if istable(pos)
    detx = pos.detx;
    dety = pos.dety;
    ionIdx = pos.ionIdx;
    pos = [double(pos.x), double(pos.y), double(pos.z)];

    % if epos as matrix is used, can get detector hit coordinates from there
elseif length(pos(1,:) == 9) && ~exist('detx','var')
    detx = pos(:,8);
    dety = pos(:,9);
    pos = double(pos(:,1:3));
    if not(exist(ionIdx))
            ionIdx = (1:length(detx))';
    end
end

reconstructionCoordinates = double(reconstructionCoordinates);
detx = double(detx);
dety = double(dety);


%% find delaunay simplex each point is in
disp('calculating delaunay triangulation of posfile'); tic;
dt = delaunayn(pos);
toc;

disp('finding simplices'); tic;
ti = tsearchn(pos,dt,reconstructionCoordinates);

% finding closest simplex for all points outside convex hull
outside = find(isnan(ti));
closest = dsearchn(pos,dt,reconstructionCoordinates(outside,:));

% finding a simplices this point is part of
for idx = 1:length(outside)
    isIn = find(sum(dt == closest(idx),2));
    ti(outside(idx)) = isIn(1);
end

toc;

% triangulations of reconstruciton and detector space
reconDel = triangulation(dt,pos);

detectorDel = triangulation(dt,detx, dety, ionIdx);

%% barycentric conversions
disp('barycentric conversion'); tic;
% calculate barycentric coordinates of each vertex in recon space
B = cartesianToBarycentric(reconDel,ti,reconstructionCoordinates);

% convert barycentric coordinates into detector coordinates
detectorCoordinates = barycentricToCartesian(detectorDel,ti,B);

toc;

end



% alternative function approximating the reconstruction coordinates with
% the closest point in the dataset instead of a barycentric transform
%function [detx, dety, N] = object2detectorCoords(fv,pos,detxP, detyP)
%translates the coordinates of the object fv (struct with fv.vertices and
%fv.faces) into detector coordinates for co - reconstruction. pos is the
%pos-file that was reconstructed and detxP and detyP are the detector hit
%coordinates of that posfile.

%pos = pos(:,1:3);


% find closest point in pos for each vertex of fv.
%[N dist] = dsearchn(pos, fv.vertices);


% take detector hit coordinate and sequence number of that point.
%detx = detxP(N);
%dety = detyP(N);


%end