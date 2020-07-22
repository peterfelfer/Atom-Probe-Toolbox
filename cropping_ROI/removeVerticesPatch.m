function fvcOut = removeVerticesPatch(fvcIn,removeVerticeList)
%% Remove Vertices from FV structure containing patch data. 
% This function will take in a FV structure of patch data and remove any
% repeated or unused vertices in fvcIn.vertices. It will find vertices in
% fvcIn that match points in the removeVerticeList, and will also remove
% any faces in fvcIn that reference vertices to be removed. It then check
% to see if any further vertices became unused after faces were removed,
% and remove them as well. If fvcIn contains field 'facevertexcdata', that
% field will be also updated. This function was created to modify patches
% before using code to convert patch to STL file; ran into problems with
% STL files generated from patches with unused vertices. (Note, developed
% in MATLAB R2013a)
%
% function fvcOut = removeVerticesPatch(fvcIn,removeVerticeList)
%
% Input: fvcIn - Structure with fields 'vertices' and 'faces', 
%                 optionally 'facevertexcdata'
%        removeVerticeList - Nx3 matrix of (x,y,z) coordinates, columns 
%                 [X Y Z]
%
% Output: fvcOut - Structure same as fvcIn with modified 'vertices','faces'
%                   and optionally 'facevertexcdata' field values
%
% Example: Remove points from a patch created from a surf plot of peaks
%     [X,Y,Z] = peaks(100);
%     % R = sqrt(X.^2 + Y.^2 + Z.^2); % distance from origin 
%     R = sqrt(X.^2 + Y.^2); % distance from origin on (x,y) plane
%     cutOff = 1;
%     fvc = surf2patch(X,Y,Z,Z);
%     ptsRemove = [X(R<cutOff),Y(R<cutOff),Z(R<cutOff)];
%     fvc = surf2patch(X,Y,Z,Z)
%     fvc1 = removeVerticesPatch(fvc,ptsRemove)
% 
%     figure
%     subplot(1,2,1)
%     h1 = patch(fvc);
%     shading faceted
%     view(3)
%     axis tight
%     daspect([1,1,2])
%     grid on
% 
%     subplot(1,2,2)
%     h2 = patch(fvc1);
%     shading faceted
%     view(3)
%     axis tight
%     daspect([1,1,2])
%     grid on
%
% Lane Foulks - Nov 2015

%Copyright (c) 2015, Lane Foulks
%All rights reserved.
%
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

%% Input Validation
assert(isstruct(fvcIn),'First input must be a structure.')
assert(all(isfield(fvcIn,{'vertices','faces'})),'First input must have both ''vertices'' and ''faces'' as fields.');
assert(isempty(removeVerticeList)||size(removeVerticeList,2)==3,'Second input must be empty or Nx3')

fvcOut = fvcIn; % copy all other fields of fvcIn to the output fvcOut
v = fvcIn.vertices;
f = fvcIn.faces;

if isempty(removeVerticeList)
    removeVerticeList = [nan nan nan]; % nan == nan always returns false, will find no matching points
end

removeCdata = false; % patch created by surf2patch may have field 'facevertexcdata'
if isfield(fvcIn,'facevertexcdata')
    removeCdata = true;
    fc = fvcIn.facevertexcdata;
end


%% Remove any duplicate vertices in patch

[vUnique, ~, indexn] =  unique(v, 'rows'); % [c,ia,ic] = unique(a), c = a(ia), a = c(ic)
fnew = indexn(f); % f (faces) ARE indices of v (vertices), maps f to unique vertices
v = vUnique;
f = fnew;   


%% find the position (rows) of the vertices to be deleted

[~,IndexVerticesDelete] = ismember(removeVerticeList,v,'rows'); % row numbers in v to remove
logicalRemoveVertices = ismember([1:length(v)],IndexVerticesDelete); % logical vector: true for rows to be removed

%%  Create new vertice reference tags

vnew = v;
tagsOld = 1:length(v);

tagsNew = tagsOld; % copy original tags
newCount = cumsum(~logicalRemoveVertices); % counts the vertices that are to remain
tagsNew(~logicalRemoveVertices) = newCount(~logicalRemoveVertices); % the newCount are the new reference tags
tagsNew(logicalRemoveVertices) = nan; % place nans at loctions that shouldn't be used

% tagMap = [tagsOld1' tagsNew1']; % map face references to new tags (remaining vertice numbering)

%% Delete vertices

vnew(logicalRemoveVertices,:) = []; 
if removeCdata
    fcnew = fc;
    fcnew(logicalRemoveVertices,:) = []; 
end

%% Delete faces

fnew = f; 
[IndexFacesRowDelete,~] = find(ismember(f,IndexVerticesDelete)); % find the position (rows) of the faces to delete
fnew(IndexFacesRowDelete,:) = []; % deletes faces that reference any vertice removed

%% Renumber faces

% fnewRenumbered = tagMap(fnew(:),2); % 2nd column of tagMap is newTags
% fnewRenumbered = reshape(fnewRenumbered,size(fnew));

fnewRenumbered = tagsNew(fnew); % maintains shape of fnew, same as above 2 lines

%% New Patch

fvcOut.faces = fnewRenumbered;
fvcOut.vertices = vnew;

if removeCdata
    fvcOut.facevertexcdata = fcnew;
end

%% Recursively call to Remove additional Vertices that are NOW unreferenced by the removal of faces

tagsNew = 1:length(vnew); % same as removing the nans - much faster than: tagsNew1(isnan(tagsNew1)) = [];
unusedVerticeList = tagsNew(~ismember(tagsNew,fnewRenumbered));
if ~isempty(unusedVerticeList)
%     fprintf(1,'Unused vertices = %d\n',length(unusedVerticeList));
    fvcOut = removeVerticesPatch(fvcOut,fvcOut.vertices(unusedVerticeList,:)); % recursive function call!!
end

