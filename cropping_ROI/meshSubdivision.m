%Copyright (c) 2016, Ankur Pawar
%All rights reserved.
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.


function [W, G] = meshSubdivision(V, F)
%function divides a triangular mesh represented by V, F
%linearly. Output vertices W are unique
%Output faces G are unique.
% V is n by 3 array of input vertices  
% F is m by 3 array of input faces
% W is array of output vertices 
% G is array of output faces
%
%      v1                   v1
%     / \                  / \
%    /   \      ->        a---c
%   /     \              / \ / \
% v2 ----- v3          v2-- b --v3
%
% start with checking input arguments
if nargin ~= 2
    error('wrong number of arguments');
end

if (size(V, 2) ~= 3) 
    error('vertices should contain 3 columns');
elseif (size(F, 2) ~= 3)
    error('faces should contain 3 columns');
end

nV = size(V,1); %number of vertices
nF = size(F,1); %number of faces

Edges = zeros(3 * nF, 2); %assuming 3 edges on each face
edgeMid = false(3 * nF,1);

e = [0 0;0 0;0 0];
r = [0 0 0];
triMesh = triangulation(F,V);
Edges = edges(triMesh);
MidPoints = (V(Edges(:,1),:) + V(Edges(:,2),:))/2; %mid points

nMidPoints = size(MidPoints, 1);

nW = nV + nMidPoints; %number of new vertices
nG = 4 * nF; %number of new faces
G = zeros(nG, 3); %new faces array
W = zeros(nW, 3); %new vertices array

W(1:nV,:) = V;
W((nV+1):(nV+nMidPoints),:) = MidPoints;


k = [1 1 1];
ind = [1 1 1];
for n = 1 : nF
    aface = F(n,:);
    %find and sort edged in a face
    if aface(1) < aface(2)
        e(1,1) = aface(1);
        e(1,2) = aface(2);
    else
        e(1,2) = aface(1);
        e(1,1) = aface(2);
    end
    
    if aface(2) < aface(3)
        e(2,1) = aface(2);
        e(2,2) = aface(3);
    else
        e(2,2) = aface(2);
        e(2,1) = aface(3);
    end
    
    if aface(1) < aface(3)
        e(3,1) = aface(1);
        e(3,2) = aface(3);
    else
        e(3,2) = aface(1);
        e(3,1) = aface(3);
    end
    
    k(1) = findaRow(Edges, e(1,:));
    k(2) = findaRow(Edges, e(2,:));
    k(3) = findaRow(Edges, e(3,:));
    
    ind(1:3) = F(n,:);
    
    G(4*(n-1)+1,1) = ind(1);
    G(4*(n-1)+1,2) = nV + k(1);
    G(4*(n-1)+1,3) = nV + k(3);
    
    G(4*(n-1)+2,1) = nV + k(1);
    G(4*(n-1)+2,2) = ind(2);
    G(4*(n-1)+2,3) = nV + k(2);
    
    G(4*(n-1)+3,1) = nV + k(2);
    G(4*(n-1)+3,2) = ind(3);
    G(4*(n-1)+3,3) = nV + k(3);
    
    G(4*(n-1)+4,1) = nV + k(1);
    G(4*(n-1)+4,2) = nV + k(2);
    G(4*(n-1)+4,3) = nV + k(3);
end
end

function k = findaRow(arr1, arr2)
%find arr2 that is row vector in another vector arr1
rows = size(arr1, 1);
k = 1;
indexFound = false;
n = 1;
while (n <= rows) && (indexFound == false)
    if (arr1(n,1) == arr2(1)) && (arr1(n,2) == arr2(2))
        k = n;
        indexFound = true;
    end
    n = n + 1;
end

end