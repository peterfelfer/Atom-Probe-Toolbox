function fv = patchCreateSampledAlphaHull(pos,alpha,sample)
%   MIGRATE TO MATLAB BUILTIN ALPHA HULL FUNCTION
% calculates the alpha hull of a set of coordinates. The alpha hull is a
% hull where concavities up to a radius of alpha are allowed. In order to
% have a mesh that is computed in a reasonable time and with reasonable
% speed, a sample of the points is taken. The default is 200k points and a
% radius of 20. The output is a Matlab 'patch' (a struct with fv.vertices
% and fv.faces).
%
% fv = patchCreateSampledAlphaHull(pos)
% alpha defaults to 20
% sample defaults to 2E5
%
% fv = patchCreateSampledAlphaHull(pos, alpha)
% sample dafaults to 2E5
%
% fv = patchCreateSampledAlphaHull(pos, alpha, sample)
%
% INPUTS:
%       pos: table with reconstructed atom positions (x, y, z) 
%       
%       alpha: value of the radius, up to which concavities are allowed 
%
%       sample: number of points of the sample
%
% OUTPUT:
%       fv: structure with fv.vertices and fv.faces (triangulated)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

if ~exist('alpha','var')
    alpha = 20;
end

if ~exist('sample','var')
    sample = 2E5;
end


vis = false;

if istable(pos)
    pos = [pos.x, pos.y, pos.z, pos.mc];
end

if length(pos(:,1)) < sample
    sample = 1:length(pos(:,1));
else
    sample = randsample(length(pos(:,1)),sample);
end

pos = pos(:,1:3);

samplePos = pos(sample,:);

% calculating boundary
[S V] = alphavol(double(samplePos),alpha,vis);

disp(['Volume of dataset: ' num2str(S) ' nm^3']);

% boundary triangles
bnd = V.bnd;

% recalculating indices
[C, ia, ic] = unique(bnd);

idx = bnd(ia);
% idx(:,2) = 1:length(idx);

for i = 1:length(bnd(:))
    bnd(i) = find(idx == bnd(i));
end    
    


% defining mesh object
fv.vertices = samplePos(idx(:,1),:);
fv.faces = bnd;

% f = figure('Name', 'Dataset hull');
% patch(fv,'FaceColor',[0 1 1],'FaceAlpha',.2); rotate3d on; axis equal;



function [V,S] = alphavol(X,R,fig)

if nargin < 2 || isempty(R), R = inf; end
if nargin < 3, fig = 0; end

% check coordinates
dim = size(X,2);
if dim < 2 || dim > 3
    error('alphavol:dimension','X must have 2 or 3 columns.')
end

% check probe radius
if ~isscalar(R) || ~isreal(R) || isnan(R)
    error('alphavol:radius','R must be a real number.')
end

% unique points
[X,imap] = unique(X,'rows');

% Delaunay triangulation
T = delaunay(X);

% remove zero volume tetrahedra since
% these can be of arbitrary large circumradius
if dim == 3
    n = size(T,1);
    vol = volumes(T,X);
    epsvol = 1e-12*sum(vol)/n;
    T = T(vol > epsvol,:);
    holes = size(T,1) < n;
end

% limit circumradius of simplices
[~,rcc] = circumcenters(TriRep(T,X));
T = T(rcc < R,:);
rcc = rcc(rcc < R);

% volume/area of alpha shape
vol = volumes(T,X);
V = sum(vol);

% return?
if nargout < 2 && ~fig
    return
end

% turn off TriRep warning
warning('off','MATLAB:TriRep:PtsNotInTriWarnId')

% alpha shape boundary
if ~isempty(T)
    % facets referenced by only one simplex
    B = freeBoundary(TriRep(T,X));
    if dim == 3 && holes
        % the removal of zero volume tetrahedra causes false boundary
        % faces in the interior of the volume. Take care of these.
        B = trueboundary(B,X);
    end
else
    B = zeros(0,dim);
end

% plot alpha shape
if fig
    if dim == 2
        % plot boundary edges and point set
        x = X(:,1);
        y = X(:,2);
        plot(x(B)',y(B)','r','linewidth',2), hold on
        plot(x,y,'k.'), hold off
        str = 'Area';
    elseif ~isempty(B)
        % plot boundary faces
        trisurf(TriRep(B,X),'FaceColor','red','FaceAlpha',1/3);
        str = 'Volume';
    else
        cla
        str = 'Volume';
    end
    axis equal
    str = sprintf('Radius = %g,   %s = %g',R,str,V);
    title(str,'fontsize',12)
end

% turn on TriRep warning
warning('on','MATLAB:TriRep:PtsNotInTriWarnId')

% return structure
if nargout == 2
    S = struct('tri',imap(T),'vol',vol,'rcc',rcc,'bnd',imap(B));
end



%--------------------------------------------------------------------------
function vol = volumes(T,X)
% VOLUMES volumes/areas of tetrahedra/triangles

% empty case
if isempty(T)
    vol = zeros(0,1);
    return
end

% local coordinates
A = X(T(:,1),:);
B = X(T(:,2),:) - A;
C = X(T(:,3),:) - A;
    
if size(X,2) == 3
    % 3D Volume
    D = X(T(:,4),:) - A;
    BxC = cross(B,C,2);
    vol = dot(BxC,D,2);
    vol = abs(vol)/6;
else
    % 2D Area
    vol = B(:,1).*C(:,2) - B(:,2).*C(:,1);
    vol = abs(vol)/2;
end


%--------------------------------------------------------------------------
function B = trueboundary(B,X)
% TRUEBOUNDARY True boundary faces
%   Remove false boundary caused by the removal of zero volume
%   tetrahedra. The input B is the output of TriRep/freeBoundary.

% surface triangulation
facerep = TriRep(B,X);

% find edges attached to two coplanar faces
E0 = edges(facerep);
E1 = featureEdges(facerep, 1e-6);
E2 = setdiff(E0,E1,'rows');

% nothing found
if isempty(E2)
    return
end

% get face pairs attached to these edges
% the edges connect faces into planar patches
facelist = edgeAttachments(facerep,E2);
pairs = cell2mat(facelist);

% compute planar patches (connected regions of faces)
n = size(B,1);
C = sparse(pairs(:,1),pairs(:,2),1,n,n);
C = C + C' + speye(n);
[~,p,r] = dmperm(C);

% count planar patches
iface = diff(r);
num = numel(iface);

% list faces and vertices in patches
facelist = cell(num,1);
vertlist = cell(num,1);
for k = 1:num

    % neglect singel face patches, they are true boundary
    if iface(k) > 1
        
        % list of faces in patch k
        facelist{k} = p(r(k):r(k+1)-1);
        
        % list of unique vertices in patch k
        vk = B(facelist{k},:);
        vk = sort(vk(:))';
        ik = [true,diff(vk)>0];
        vertlist{k} = vk(ik);
        
    end
end

% sort patches by number of vertices
ivert = cellfun(@numel,vertlist);
[ivert,isort] = sort(ivert);
facelist = facelist(isort);
vertlist = vertlist(isort);

% group patches by number of vertices
p = [0;find(diff(ivert));num] + 1;
ipatch = diff(p);

% initiate true boundary list
ibound = true(n,1);

% loop over groups
for k = 1:numel(ipatch)

    % treat groups with at least two patches and four vertices
    if ipatch(k) > 1 && ivert(p(k)) > 3

        % find double patches (identical vertex rows)
        V = cell2mat(vertlist(p(k):p(k+1)-1));
        [V,isort] = sortrows(V);
        id = ~any(diff(V),2);
        id = [id;0] | [0;id];
        id(isort) = id;

        % deactivate faces in boundary list
        for j = find(id')
            ibound(facelist{j-1+p(k)}) = 0;
        end
        
    end
end

% remove false faces
B = B(ibound,:);


