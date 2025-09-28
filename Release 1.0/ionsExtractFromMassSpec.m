function ionTable = ionsExtractFromMassSpec(spec)
% ionsExtractFromMassSpec pulls all ions and corresponding information 
% from a mass spectrum plot and gets all plots connected to the mass 
% spectrum
%
% ionTable = ionsExtractFromMassSpec(spec)
%
% INPUT
% spec: area plot that displays the mass spectrum (histogram of 
%       m/c frequencies)either in raw counts or normalised to bin width 
%       and total ion count
%
% OUTPUT 
% ionTable: table with allocated ions and additional information
%           (charge state, corresponding color code)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%%
plots = spec.Parent.Children;

idx = 0;


for pl = 1:length(plots)
    
    % find all the ones that are ranges
    try
        type = plots(pl).UserData.plotType;
    catch
        type = "unknown";
    end
    
    if type == "ion"
        idx = idx + 1;
        ion{idx,:} = plots(pl).UserData.ion{1}.element; % do not extract isotope information
        chargeState(idx,:) = plots(pl).UserData.chargeState(1);
        ionName{idx,:} = ionConvertName(plots(pl).UserData.ion{1}.element);
        color(idx,:) = plots(pl).Color;
        
    end
    
end

if idx > 0
    ionName = categorical(ionName);
    ionTable = table(ionName,chargeState,ion,color);
else
    ionName = [];
    ionTable = [];
end
