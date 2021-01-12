function [fv energy] = patchSnapToAtomPositions(fv,pos,rad)
%   NEEDS DOCUMENTATION AND UPDATES
%performs a DCOM step with influence radius rad on the surface fv relative
%to the positions in pos for an object of dimension d. If the noromals are
%parsed, fv is just a vertex list.

%In a DCOM step,each vertex v of a surface fv is moved to the center of
%mass of all points in pos that are closer to v than rad.

%the DCOM energy is the sum of the square of all displacements - for
%convergence measurments.


addpath('utilities_patch_normals');

if istable(pos)
    pos = [pos.x, pos.y, pos.z];
end

debug = false;
noTris = false;

if debug 
    fvOld = fv;
    figure('Name','surface evolution');
    patch(fv,'FaceColor',[0 1 1],'FaceAlpha',0);
    axis equal; rotate3d on;   
    tic;
end


if ~isstruct(fv)
    normals = estimateNormals(fv,pos,rad);
    fv.vertices = fv;
    noTris = true;
else
    normals = patchnormals(fv);
end


%% discard all points that cannot contibute to solution (outside bounding box of vertices + rad)
pos = pos(:,1:3);
mx = max(fv.vertices) + rad;
mn = min(fv.vertices) - rad;

in = (pos > repmat(mn,length(pos(:,1)),1)) & (pos < repmat(mx,length(pos(:,1)),1));
in = sum(in,2) == 3;

% evaluating how many atoms we have clipped out
numAtomRaw = length(pos(:,1));
discarded = sum(~in);


pos(~in,:) = [];
clear in





%% calculate distance between verts and points and DCOM step
numAtoms = length(pos(:,1));
rad2 = rad^2;




fvOld = fv;

for vert = 1:length(fv.vertices(:,1))
    % calculate relative coordinates to vertex
    relCoord = pos - repmat(fv.vertices(vert,:),numAtoms,1);
    
    % calculate if dist < rad
    in = sum(relCoord.^2,2) <= rad2;
    
    % center of mass of these points (unweighed)
    cent(vert,:) = mean(relCoord(in,:),1);
    numAt(vert) = length(relCoord(in,:));
    
    % project onto surface normal
    shift(vert) = dot(cent(vert,:),normals(vert,:));
    
    %shift vertex coordinate
    fv.vertices(vert,:) = fv.vertices(vert,:) + normals(vert,:) * shift(vert);
    
end

if numAt == 0
    warning(['Vertices ', num2str(find(numAt == 0)), ' have no associated atoms. Position calculation fails']);
end


%% calculation of DCOM energy
energy = sum(sqrt(sum((fv.vertices - fvOld.vertices).^2,2)));

if noTris
    fv = fv.vertices;
end


if debug
    hold on;
    patch(fv,'FaceColor',[0 1 1],'FaceAlpha',0.6);
    toc;
    disp(['discarded atoms (box clipping): ' num2str(discarded) ', (' num2str(discarded/numAtomRaw) '%)']);
    
    displace = fv.vertices - fvOld.vertices;
    disp(['average displacement vector:   ' num2str(mean(displace,1))]);
    disp(['average displacement distance: ' num2str(norm(mean(displace,1)))]);
    
    figure('Name','displacement vectors');
    scatter3(displace(:,1),displace(:,2),displace(:,3),'k.'); axis equal; rotate3d on;
end




end