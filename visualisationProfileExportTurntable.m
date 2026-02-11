function exportInfo = visualisationProfileExportTurntable(pos, colorScheme, visualisationProfile, options)
% VISUALISATIONPROFILEEXPORTTURNTABLE Scriptable turntable video export.
%
% exportInfo = visualisationProfileExportTurntable(pos, colorScheme, visualisationProfile, ...)

arguments
    pos table
    colorScheme
    visualisationProfile (1,1) struct
    options.outputFile (1,:) char = 'apt_turntable.mp4'
    options.stepDeg (1,1) double = NaN
    options.frameRate (1,1) double = NaN
    options.totalAngleDeg (1,1) double = NaN
    options.quality (1,1) double = 95
    options.axes = []
    options.speciesPolicy (1,1) string = ""
    options.showWidget (1,1) logical = false
    options.closeControlFig (1,1) logical = true
end

profile = visualisationProfileMigrate(visualisationProfile);
visualisationProfileValidate(profile);

if isnan(options.stepDeg)
    stepDeg = double(profile.export.turntableStepDeg);
else
    stepDeg = options.stepDeg;
end
if isnan(options.frameRate)
    frameRate = double(profile.export.turntableFrameRate);
else
    frameRate = options.frameRate;
end
if isnan(options.totalAngleDeg)
    totalAngleDeg = double(profile.export.turntableTotalAngleDeg);
else
    totalAngleDeg = options.totalAngleDeg;
end

stepDeg = max(0.01, stepDeg);
frameRate = max(1, frameRate);
totalAngleDeg = max(stepDeg, totalAngleDeg);

[~, ax, controlFig, infoApply] = visualisationProfileApply( ...
    pos, colorScheme, profile, ...
    'axes', options.axes, ...
    'showWidget', options.showWidget, ...
    'speciesPolicy', options.speciesPolicy, ...
    'restoreStateFromAxis', false, ...
    'persistStateToAxis', true);

fig = ancestor(ax, 'figure');
if isempty(fig) || ~isgraphics(fig)
    error('visualisationProfileExportTurntable:noFigure', ...
        'Could not resolve a figure for turntable export.');
end

if ~strcmp(ax.XLimMode, 'manual') || ~strcmp(ax.YLimMode, 'manual') || ~strcmp(ax.ZLimMode, 'manual')
    axis(ax, 'vis3d');
    ax.XLimMode = 'manual';
    ax.YLimMode = 'manual';
    ax.ZLimMode = 'manual';
end

[outDir, ~, ext] = fileparts(options.outputFile);
if ~isempty(outDir) && ~isfolder(outDir)
    mkdir(outDir);
end
if isempty(ext)
    ext = '.mp4';
    options.outputFile = [options.outputFile, ext];
else
    ext = lower(ext);
end

switch ext
    case '.avi'
        writer = VideoWriter(options.outputFile, 'Motion JPEG AVI');
    case '.mp4'
        writer = VideoWriter(options.outputFile, 'MPEG-4');
    otherwise
        error('visualisationProfileExportTurntable:unsupportedExt', ...
            'Unsupported video extension ''%s''. Use .mp4 or .avi.', ext);
end
writer.FrameRate = frameRate;
if isprop(writer, 'Quality')
    writer.Quality = max(1, min(100, options.quality));
end

nFrames = max(1, round(totalAngleDeg / stepDeg));

camState = captureCameraState(ax);
cleanupCam = onCleanup(@() restoreCameraState(ax, camState)); %#ok<NASGU>

open(writer);
cleanupWriter = onCleanup(@() closeVideoWriterSafe(writer)); %#ok<NASGU>

for k = 1:nFrames
    frame = getframe(fig);
    writeVideo(writer, frame);
    rotateAxesAroundZ(ax, stepDeg);
    drawnow;
end

exportInfo = struct();
exportInfo.outputFile = string(options.outputFile);
exportInfo.nFrames = nFrames;
exportInfo.stepDeg = stepDeg;
exportInfo.frameRate = frameRate;
exportInfo.applyInfo = infoApply;

if options.closeControlFig && isgraphics(controlFig)
    close(controlFig);
end

end

function state = captureCameraState(ax)
state = struct();
props = {'XLim','YLim','ZLim','DataAspectRatio','PlotBoxAspectRatio', ...
    'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle', ...
    'XLimMode','YLimMode','ZLimMode', ...
    'DataAspectRatioMode','PlotBoxAspectRatioMode', ...
    'CameraPositionMode','CameraTargetMode','CameraUpVectorMode','CameraViewAngleMode'};
for i = 1:numel(props)
    p = props{i};
    try
        state.(p) = ax.(p);
    catch
    end
end
end

function restoreCameraState(ax, state)
if ~isgraphics(ax)
    return;
end
props = fieldnames(state);
for i = 1:numel(props)
    p = props{i};
    if isprop(ax, p)
        try
            ax.(p) = state.(p);
        catch
        end
    end
end
end

function rotateAxesAroundZ(ax, stepDeg)
if ~isgraphics(ax)
    return;
end
try
    camorbit(ax, stepDeg, 0, 'data', [0 0 1]);
catch
    axes(ax);
    camorbit(stepDeg, 0);
end
end

function closeVideoWriterSafe(writer)
if isempty(writer)
    return;
end
try
    close(writer);
catch
end
end
