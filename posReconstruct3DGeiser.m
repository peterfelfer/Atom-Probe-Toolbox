function [pos, objects] = posReconstruct3DGeiser(pos,flightPathLength,...
    detectorEfficiency,effectiveDetectorArea,ionVolume,radiusEvolution,...
    ICF,objects)
%atom probe reconstruction after: Gault et al., Ultramicroscopy 111 (2011) 448 - 457
%detx, dety are the detector hit coordinates in mm
%kf is the field factor and ICF is the image compression factor
% for reflectron data an image correction is required!

%% fundamental data reconstruction
% detector coordinates in polar form
[ang, rad] = cart2pol(pos.detx, pos.dety);

% launch angle relaive to specimen axis
thetaP = atan(rad / flightPathLength); % mm/mm

% image compression correction
theta = thetaP + asin((ICF - 1) * sin(thetaP));

% distance from axis and z shift of each hit
[zP, d] = pol2cart(theta, radiusEvolution); % nm

% x and y coordinates from the angle on the detector and the distance to
% the specimen axis.
[pos.x, pos.y] = pol2cart(ang, d); % nm

%% calculate z coordinate
% the z shift with respect to the top of the cap is Rspec - zP
zP = radiusEvolution - zP;

% accumulative part of z
omega = 1 ./ ionVolume; % atmic volume in nm^3

% magnification M at ion index
M = flightPathLength./(ICF * radiusEvolution);
% currently evaporating area of the specimen
specArea = effectiveDetectorArea ./ M.^2;
%individual depth increment
dz = omega ./ specArea;

% wide angle correction
cumZ = cumsum(double(dz));
pos.z = cumZ + zP;
clear ang rad d zP




%% co-reconstruction of objects
if exist('objects','var')
    for o = 1:length(objects)
        % detector coordinates in polar form
        [ang, rad] = cart2pol(objects(o).vertices(:,1),objects(o).vertices(:,2));
        
        % launch angle relaive to specimen axis
        thetaP = atan(rad / flightLength); % mm/mm
        theta = thetaP + asin((ICF - 1) * sin(thetaP));
        
        % distance from axis and z shift of each hit
        [zP, d] = pol2cart(theta, Rspec(round(objects(o).vertices(:,3)))); % nm
        
        % x and y coordinates from the angle on the detector and the distance to
        % the specimen axis.
        [objectsR(o).vertices(:,1), objectsR(o).vertices(:,2)] = pol2cart(ang, d); % nm
        
        
        % the z shift with respect to the top of the cap is Rspec - zP
        zP = Rspec(round(objects(o).vertices(:,3))) - zP;
        
        objectsR(o).vertices(:,3) = cumZ(round(objects(o).vertices(:,3))) + zP;
        
        objectsR(o).faces = objects(o).faces;
        
    end
end
