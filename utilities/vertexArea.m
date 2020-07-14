function A = vertexArea(fv,pos,clipD)
%calculate the area each vertex occupies. If a pos file is given, a prism
%approximation with a clipping distance clipD is used (better for atom probe data with density
%fluctuations and / or a boundary)


DEBUG = true;

addpath('PointTool','toolbox_graph','fitNormal','patch_normals');

remesh.vertices = fv.vertices;
remesh.faces = fv.faces;

numVerts = length(remesh.vertices(:,1));



if ~exist('pos','var')
    
    %% calcualtion of vertex area based on mesh
    %if DEBUG
    %    scatter3(remesh.vertices(:,1),remesh.vertices(:,2),remesh.vertices(:,3),'MarkerFaceColor','r','MarkerEdgeColor','k','SizeData',72);
    %end
    
    % recalculating triangulation
    
    %% detect boundary and resample boundary
    edgeList = tri2edgeList(remesh.faces);
    
    
    if DEBUG
        beg = 53;%randsample(60,1)
        numVerts = beg;
        
        CVTsteps = 1
    else
        beg = 1;
    end
    
    
    if DEBUG
        disp('debugging mode');
        ph = patch(fv,'CData',[0 1 0],'FaceAlpha',0.6,'FaceColor',[0 1 0],'EdgeColor',[0 0 0]);
        ah = gca;
        axis equal
        rotate3d on
        hold on
    end
    
    
    normals = patchnormals(remesh);
    
    wb = waitbar(0,['calculating vertex areas']);
    
    
    for vert = beg:numVerts
        
        % get polygon around that vertex (all edges that contain the point)
        currentVert = remesh.vertices(vert,:);
        
        neighbours = sum(edgeList == vert,2);
        neighbours = edgeList(logical(neighbours),:);
        neighbours = unique(neighbours);
        neighbours = neighbours(neighbours~=vert);
        
        neighbourVerts = remesh.vertices(neighbours,:);
        
        %poly = remesh.vertices(neighbours,:);
        
        % project polygon into 2D space
        normal = normals(vert,:);
        
        if DEBUG
            
            scatter3(ah,currentVert(:,1),currentVert(:,2),currentVert(:,3),'r','SizeData',128,'MarkerFaceColor','r');
            scatter3(ah,neighbourVerts(:,1),neighbourVerts(:,2),neighbourVerts(:,3),'b','SizeData',128,'MarkerFaceColor','b');
            
        end
        
        R = null(normal);
        
        uv = [neighbourVerts; currentVert] * R;
        
        % to be fixed (hull is not generally convex)
        
        [edg Atemp] = convhull(uv);
        
        uv = uv(edg,:);
        
        % plot all 2d properties
        if DEBUG
            %figure();
            %voronoi([uv(:,1); orig(1)],[uv(:,2); orig(2)]);
            disp(['area = ' num2str(Atemp)]);
            figure();
            plot(uv(:,1),uv(:,2),'-o');
            hold on
            axis equal
        end
        
        A(vert) = Atemp;
        
        
        waitbar(vert/numVerts,wb);
        
    end
    
    close(wb);



else
    %% calcualtion of vertex area based on pos file (prism approximation)
    LIM = 1E5;%maximum amount of atoms for the alpha hull 
    RAD = 10; %nm radius for alpha hull the volume is calcualted for

    if ~exist('clipD','var')
        clipD = 2;%nm
    end

    normals = patchnormals(remesh);
    closest = dsearchn(fv.vertices,delaunayn(fv.vertices),pos(:,1:3));
    distVec = pos(:,1:3) - fv.vertices(closest,:);
    %distance through dot product
    dist = abs(sum(normals(closest,:) .* distVec,2));

    closest = closest(dist <= clipD);




    % calculating the atomic density of the dataset
    numAtom = length(pos);
    if numAtom > LIM
        datasetVol = alphavol(pos(randsample(1:length(pos),LIM),1:3),RAD);
    else
        datasetVol = alphavol(pos(:,1:3));
    end

    density = numAtom/datasetVol;

    
    for vert = 1:numVerts
        A(vert) = sum(closest == vert) / (density * clipD);
    end

    
    
end


end





