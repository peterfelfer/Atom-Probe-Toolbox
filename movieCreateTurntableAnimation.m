function mov = movieCreateTurntableAnimation(deg,frameRate,fileName)
% movieCreateTurntableAnimation creates a turntable animation of the current 
% figure and returns movie variable. The movie is saved in an *.avi format.
%
% mov = movieCreateTurntableAnimation(deg,frameRate,fileName)
% mov = movieCreateTurntableAnimation(deg,frameRate)
% mov = movieCreateTurntableAnimation(deg)
%
% INPUTS
% deg:       step size in degree in which the animation rotates. determines
%            the length and speed of the film.
% 
% frameRate: frame rate of the movie in frames per second
% 
% fileName:  name under which the film is saved.
%
%
%
% OUTPUT
%           movie of the turning atom probe tip.

if ~exist('frameRate','var')
    frameRate = 30;     % default value set to 30
end

if ~exist('fileName','var')
    [file path] = uiputfile('*.avi','export movie');
    fileName = [path file];
end

set(gcf,'Color','w');
axis vis3d

nFrames = 360/deg;

for f = 1:nFrames
    frame = getframe(gcf);
    mov(f) = frame;
    camorbit(deg,0);
end

myVideo = VideoWriter(fileName);

myVideo.FrameRate = frameRate;  
myVideo.Quality = 100;          % default set to 100

open(myVideo);
writeVideo(myVideo, mov);
close(myVideo);
