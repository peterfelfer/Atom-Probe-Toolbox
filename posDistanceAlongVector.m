function distance = posDistanceAlongVector(pos,origin,direction)
% posDistanceAlongVector calculates the distance of entries in a pos file along a specified
% vector with an origin of 'origin' and a direction of 'direction'
%
% distance = posDistanceAlongVector(pos,origin,direction)
%
% INPUT
% pos =         pos file with the important data 
%
% origin =      origin of the Vector
%
% direction =   direction of the Vector
%
% OUTPUT
% distance =    distance between entries in a pos file along a specified
%               vector

points = [pos.x pos.y pos.z];
points = points - repmat(origin,[height(pos),1]); % re-center origin

direction = direction/norm(direction); %make axis unit length

distance = sum(points.*repmat(direction,[height(pos),1]),2);% dot product with scale axis for projected length