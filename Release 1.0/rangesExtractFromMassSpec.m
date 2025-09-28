function rangeTable = rangesExtractFromMassSpec(spec)
% rangesExtractFromMassSpec pulls all ranges and additional information 
% from a mass spectrum plot and gets all plots connected to the mass 
% spectrum
%
% rangeTable = rangesExtractFromMassSpec(spec)
%
% INPUT
% spec:         area plot that displays the mass spectrum (histogram of m/c
%               frequencies)either in raw counts or normalised to bin width
%               and total ion count
%
% OUTPUT 
% rangeTable:   table with allocated ranges of the ions and additional 
%               information(charge state, corresponding color code)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

plots = spec.Parent.Children;



rngIdx = 0;
for pl = 1:length(plots)
    
    % find all the ones that are ranges
    try
        type = plots(pl).UserData.plotType;
    catch
        type = "unknown";
    end
    
    if type == "range"
        
        rngIdx = rngIdx +1;
        
        mcbegin(rngIdx,:) = plots(pl).XData(1);
        mcend(rngIdx,:) = plots(pl).XData(end);
        % errors in rangeAdd lead to area plots without allocated name
        if sum(strcmp(fieldnames(plots(pl).UserData),'ion')) == 0
            rangeName{rngIdx,:} = 'not assigned';
            chargeState(rngIdx,:) = NaN;
            ion{rngIdx,:} = table({'not assigned'});
            
            display(['Range between ' num2str(round(mcbegin./0.01).*0.01) ' and ' num2str(round(mcend./0.01).*0.01) ...
                    ' cannot be allocated to ion/atom!']);
        elseif istable(plots(pl).UserData.ion) 
            rangeName{rngIdx,:} = ionConvertName(plots(pl).UserData.ion.element);
            ion{rngIdx,:} = plots(pl).UserData.ion;
            chargeState(rngIdx,:) = plots(pl).UserData.chargeState;
        else
            rangeName{rngIdx,:} = plots(pl).UserData.ion;
            element = categorical(string(plots(pl).UserData.ion));
            isotope = NaN;
            ion{rngIdx,:} = table(element,isotope);
            chargeState(rngIdx,:) = NaN;
        end
        volume(rngIdx,:) = 0;
        color(rngIdx,:) = plots(pl).FaceColor;
        
        
    end
    
end

if rngIdx == 0
    rangeTable = [];
else
    rangeName = categorical(rangeName);
    rangeTable = table(rangeName,chargeState,mcbegin,mcend,volume,ion,color);
    rangeTable = sortrows(rangeTable,'mcbegin','ascend');
end