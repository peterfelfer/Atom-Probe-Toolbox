function distance = posDistanceToLine(pos,line)
% posDistanceToLine calculates the distance of individual atoms to a line 
% object representing e.g. a dislocation or triple line.
% line is a line object consisting of vertices and line segments between
% them 
% line is parsed as a struct with line.vertices (Nx3) and line.edges (Mx2)
%
%
%
% INPUT
% pos:          pos file 
%
% line:         is a line object consisting of
%               line.vertices (Nx3) and line.edges (Mx2)
%
% OUTPUT
% distance:     distance between entries in a pos file along a specified
%               line    

%% distances are calculated perpendicular to the line elements.
[lineVector, ~] = lineVectors(line.vertices,line.edges);


%% tessellation and distance calculation
% for overall pos file
% finding closest point for each atomic position
closest = dsearchn(line.vertices,delaunayn(line.vertices),[pos.x, pos.y, pos.z]);

% vector from atom to closest vertex
vec = [pos.x, pos.y, pos.z] - line.vertices(closest,1:3); 
vecLen = sqrt(sum(vec.^2,2));

% distance along line vector
distLV = dot(vec, lineVector(closest,:), 2);

% distance normal to line vector
distance = sqrt(vecLen.^2 - distLV.^2); 