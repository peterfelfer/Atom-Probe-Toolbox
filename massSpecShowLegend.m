function massSpecShowLegend(spec,items)
% massSpecShowLegende shows the legend of a mass spectrum plot, 
% with selected item(s) in focus as shown below, order of items is editable
% unknown ions will not be shown in the legend
%
% massSpecShowLegend(spec)
% massSpecShowLegend(spec,"item1")
% massSpecShowLegend(spec,["item1","item2"])
%
% INPUT
% spec:     spec is the name of the areaplot that contains the
%           mass spectrum, area
%
% items:    defines which information is written down in the legend
%           "ion": list of ions by type (isotopes not listed)
%           "range": list of ions and ranged isotopes
%           "massSpectrum": only the grey mass spectrum is listed in the
%           legend
%           default item is "ion"
%
% OUTPUT
%           shows the legend of the ranged ions in the mass spectrum
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%%
if ~exist('items','var')
    items = ["ion", "unknown"]; % "ion" "text" "massSpectrum" "range"
end

% get plot type
plots = spec.Parent.Children;
for pl = 1:length(plots)
    try
        plotType(pl,:) = plots(pl).UserData.plotType;
    catch
        plotType(pl,:) = "unknown";
    end
end

% determine which one is a selected item
isLegend = any(plotType == items,2);

legend(plots(isLegend));
    