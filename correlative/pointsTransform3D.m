function ptsOut = pointsTransform3D(pointsIn, translation, scale, rotation, rotationOrigin)
% function to shift, rotate and scale any object based on point
% coordinates. Can be an Nx3 array, a pos table (no rotation of raw data)
% or a patch (mesh object)

if istable(pointsIn)
    coords.x = pointsIn.x;
    coords.y = pointsIn.y;
    coords.z = pointsIn.z;

end

if isstruct(pointsIn)
    coords.x = pointsIn.vertices(:,1);
    coords.y = pointsIn.vertices(:,2);
    coords.z = pointsIn.vertices(:,3);
end

end