function rangeTable = rangesExtractFromMassSpec(spec)
% pulls all ranges and additional information from a mass spectrum plot
% gets all plots connected to the mass spectrum
% 
% INPUT: spec, area plot that displays the mass spectrum (histogram of m/c frequencies) 
%        either in raw counts or normalised to bin width and total ion count
%
% OUTPUT: rangeTable, table with allocated ranges of the ions and additional information
%         (charge state, corresponding color code)
%
plots = spec.Parent.Children;



idx = 1;
for pl = 1:length(plots)
    
    % find all the ones that are ranges
    try
        type = plots(pl).UserData.plotType;
    catch
        type = "unknown";
    end
    
    if type == "range"
        mcbegin(idx,:) = plots(pl).XData(1);
        mcend(idx,:) = plots(pl).XData(end);
        if istable(plots(pl).UserData.ion)
            rangeName{idx,:} = ionConvertName(plots(pl).UserData.ion.element);
            ion{idx,:} = plots(pl).UserData.ion;
            chargeState(idx,:) = plots(pl).UserData.chargeState;
        else
            rangeName{idx,:} = plots(pl).UserData.ion;
            element = categorical(string(plots(pl).UserData.ion));
            isotope = NaN;
            ion{idx,:} = table(element,isotope);
            chargeState(idx,:) = NaN;
        end
        volume(idx,:) = 0;
        color(idx,:) = plots(pl).FaceColor;

        idx = idx +1;
    end
    
end

rangeName = categorical(rangeName);
rangeTable = table(rangeName,chargeState,mcbegin,mcend,volume,ion,color);
rangeTable = sortrows(rangeTable,'mcbegin','ascend');