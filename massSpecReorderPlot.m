function massSpecReorderPlot(spec,order)
% massSpecReorderPlot rearranges the order of the various components of the
% mass spectrum in the selected order (e.g., text, ion, range,...)
% order is a cell array of strings as below
%
% massSpecReorderPlot(spec)
% massSpecReorderPlot(spec,order)
% massSpecReorderPlot(spec,["order1","order2",etc])
%
% INPUT
% spec:     spec is the name of the area plot that contains the
%           mass spectrum, area
%
% order:    set what kind of information stands in focus. if more
%           information is given, order can contain following types in any order:
%           "text": the text of the ion
%           "backgroundEstimate": background estimation line (e.g., from concentration calculation)
%           "ion": stem plot of the ion, which shows the exact position of the peak
%           "range": ranged area of the peaks
%           "background": background ranges for the mass spectrum
%           "massSpectrum": the neutral mass spectrum
%           "unknown": everything else
%           default order is (first = foreground, last = background):
%           ["text","backgroundEstimate","ion","range","background","massSpectrum","unknown"]
%
% OUTPUT
%          mass spectrum with rearranged order
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg


%% default if no input order is given
if ~exist('order','var')
    order = ["text","backgroundEstimate","ion","range","background","massSpectrum","unknown"];
end


plots = spec.Parent.Children;

for pl = 1:length(plots)
    mcbegin(pl,:) = 0;
    try
        plotType(pl,:) = plots(pl).UserData.plotType;
        
        if any(plotType(pl) == ["range", "ion"])
            mcbegin(pl,:) = plots(pl).XData(1);
        end
        
    catch
        plotType(pl,:) = "unknown";
    end
end

% ordinal categorical array is used to sort
plotType = categorical(plotType,order,'Ordinal',true);
sortTable = table(plotType,mcbegin);

[~, idx] = sortrows(sortTable);

spec.Parent.Children = spec.Parent.Children(idx);