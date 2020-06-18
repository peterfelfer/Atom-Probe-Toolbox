function massSpecReorderPlot(spec,order)
% rearranges the mass spectrum visiblity in the selected order
% order is a cell array of strings as below. 
%
% massSpecReorderPlot(spec)
% massSpecReorderPlot(spec,order)
% massSpecReorderPlot(spec,["order1","order2",etc])
%
% INPUTS
% spec:     spec is the name of the areaplot that contains the
%           massspectrum, area.
% order:    set what kind of information stands in focus. if more
%           informations are given order can contain follwoing types in any order:
%           "text": the text in which the art of the ion is written.
%           "ion": deshed line which shows exact position of element peak.
%           "range": ranged area of the peaks.
%           "background": background ranges for the mass spectrum
%           "massSpectrum": the neutral massspectrum.
%           "unknown": everything else.
%           default order is
%           ["text","ion","range","massSpectrum","unknown"]
%           
% OUTPUTS
%          Massspectrum with focus on wanted informations 



if ~exist('order','var')
    order = ["text","ion","range","background","massSpectrum","unknown"];
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