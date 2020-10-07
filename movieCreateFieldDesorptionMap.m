function fdmStack = movieCreateFieldDesorptionMap(detx,dety,frames,sample,RES,fileName)
% movieCreateFieldDesportionMap produces a movie of the field desorption 
% map as the experiment progresses, based on detector coordinates detx and dety. 
% It will be 'frames' frames and each frame will used 'sample' number of 
% hits to calculate the FDM.
%
% fdmStack = movieCreateFieldDesorptionMap(detx,dety)
% frames defaults to floor(numAtoms/sample), therefore depends on the 
% size/length of the pos variable (--> numAtoms) and the sample size, 
% which defaults to 1M ions; RES defaults to 128; fileName defaults to FDM
%
% fdmStack = movieCreateFieldDesorptionMap(detx,dety,frames)
% sample defaults to 1M ions; RES defaults to 128; fileName defaults to FDM
%
% fdmStack = movieCreateFieldDesorptionMap(detx,dety,frames,sample)
% RES defaults to 128; fileName defaults to FDM
%
% fdmStack = movieCreateFieldDesorptionMap(detx,dety,frames,sample,RES)
% fileName defaults to FDM
%
% fdmStack = movieCreateFieldDesorptionMap(detx,dety,frames,sample,RES,fileName)
%
%
% INPUTS:
%       detx: x coordinate of each ion hit on the detector, e.g., epos.detx
%       
%       dety: y coordinate of each ion hit on the detector, e.g., epos.dety
%
%       frames: number of created field desoprtion map frames, which assemble the movie, optional 
%
%       sample: number of ions for each field desorption map frame; limited
%               by the resolution of the field desorption map, optional
%
%       RES: resolution of the FDM, optional
%
%       fileName: string, name of the movie file with an .avi suffix; saved in the
%                 current folder, optional
% 
% OUTPUTS:
%       fdmStack: matrix of size (frames)x(resolution)x(resolution)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg



numAtoms = length(detx);

if ~exist('sample','var')
    % the sample will be 1M ions
    sample = 1E6;
end

if ~exist('frames','var')
    frames = floor(numAtoms/sample);
end

if ~exist('RES','var')
    RES = 128;      
end

if ~exist('fileName', 'var')
    fileName = 'FDM';
end

%% calculating bin centers for FDM

mi = min(min(detx),min(dety));
mx = max(max(detx),max(dety));

mi = mi * 1.1;
mx = mx * 1.1;

ctr{1} = linspace(mi,mx,RES);
ctr{2} = ctr{1};




%% determining sequence ranges for FDM
step = numAtoms/frames;

for f = 1:frames
    range(f,1) = step * (f-1) +1;
    range(f,2) = range(f,1) + sample;
    if range(f,2) > numAtoms
        range(f,2) = numAtoms;
    end
end

range = round(range);


%% calculating individual frames

fig = figure;

% aviobj = avifile('FDM.avi','fps',FPS);
v = VideoWriter(fileName);
open(v);
for f = 1:frames
    FDM = hist3([detx(range(f,1):range(f,2)), dety(range(f,1):range(f,2))],ctr);
    fdmStack(f,:,:) = FDM;
    im = imagesc(FDM); axis equal;
%     mov(f) = getframe(gca);
    frame = getframe(gca);
%     aviobj = addframe(aviobj,frame);
    writeVideo(v,frame);
end
close(v);
% aviobj = close(aviobj);

% movie(mov);