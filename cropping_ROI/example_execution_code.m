%example code for post move execution

% example ROI object: roi
% example data: sol
% example plot: roiplot

% crops data to ROI object (convex only!)
in = inhull(sol(:,1:3),roi.Vertices);

% transforms data to ROI object's coordinate system
ptsIn = coordSystemTransform(roi, sol(in,:));

% updates plot data
roiplot.XData = ptsIn(:,1); roiplot.YData = ptsIn(:,2); roiplot.ZData = ptsIn(:,3);

% updates 1D concentration profile in plot 'lp'
 [profile, binVec] = calculate1Dprofile(roi.UserData.ROIzaxis,bin,solutes,allAtoms); 
 lp.XData = binVec; lp.YData = profile;