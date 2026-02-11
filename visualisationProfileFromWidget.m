function visualisationProfile = visualisationProfileFromWidget(source)
% VISUALISATIONPROFILEFROMWIDGET Build a visualisation profile from widget/axes.
%
% visualisationProfile = visualisationProfileFromWidget()
% visualisationProfile = visualisationProfileFromWidget(controlFig)
% visualisationProfile = visualisationProfileFromWidget(ax)

if nargin < 1
    source = [];
end

state = scatterPlotPosWidgetGetState(source);
visualisationProfile = visualisationProfileMigrate(state);
visualisationProfileValidate(visualisationProfile);

end
