function [x, y, z, objectsR] = posGeiserReconstructFromDetector(detx, dety, V, kf, ICF,objects)
%atom probe reconstruction after: Gault et al., Ultramicroscopy 111 (2011) 448 - 457
%detx, dety are the detector hit coordinates in mm
%kf is the field factor and ICF is the image compression factor
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg



%% constants and variable setup

% instrument parameters
flightLength = 110; % flight path length in mm
detEff = 0.82; % detector efficiency

% specimen parameters
avgDens = 60.2; % atomic density in atoms / nm3
Fevap = 65; % evaporation field in V/nm

%V = V * 1000; %voltage provided in kV

% detector coordinates in polar form
[ang, rad] = cart2pol(detx, dety);

% calcualting effective detector area:
Adet = (max(rad))^2 * pi();

% radius evolution from voltage curve (in nm)
Rspec = V/(kf * Fevap);



%% calcualte x and y coordinates

% launch angle relaive to specimen axis
thetaP = atan(rad / flightLength); % mm/mm
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



%% object reconstruction
clear ang rad d zP

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


%% end of reconstruction


