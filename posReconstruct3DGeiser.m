function pos = posReconstruct3DGeiser(pos,flightPathLength,detectorEfficiency,effectiveDetectorArea,ionVolume,radiusEvolution,ICF,objects)
%atom probe reconstruction after: Gault et al., Ultramicroscopy 111 (2011) 448 - 457
%detx, dety are the detector hit coordinates in mm
%kf is the field factor and ICF is the image compression factor
% for reflectron data an image correction is required!

%% fundamental data reconstruction
% detector coordinates in polar form
[ang, rad] = cart2pol(detx, dety);

% launch angle relaive to specimen axis
thetaP = atan(rad / flightLength); % mm/mm

% image compression correction
theta = thetaP + asin((ICF - 1) * sin(thetaP));

% distance from axis and z shift of each hit
[zP, d] = pol2cart(theta, Rspec); % nm

% x and y coordinates from the angle on the detector and the distance to
% the specimen axis.
[x, y] = pol2cart(ang, d); % nm

%% calculate z coordinate
% the z shift with respect to the top of the cap is Rspec - zP
zP = Rspec - zP;

% accumulative part of z
omega = 1 ./ avgDens; % atmic volume in nm^3

dz = omega * flightLength^2 * kf^2 * Fevap^2 / (detEff * Adet * ICF^2) * V.^-2; % nm^3 * mm^2 * V^2/nm^2 / (mm^2 * V^2)

% wide angle correction
cumZ = cumsum(double(dz));
z = cumZ + zP;

%% co-reconstruction of objects