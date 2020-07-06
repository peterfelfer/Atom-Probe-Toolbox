function [h, txt] = rangeAdd(spec,colorScheme,manualName)
% adds a range to a mass spectrum using graphical input
% output is the handle to the area plot and the corresponding text
%
% if multiple isotopic combinations of the same element are within the range,
% automatically the one with the higher abundance (peak height) will be taken
%
% if various ions of different element species are within one range,
% the user can choose the desired ion out of a list in a pop-up window;
% default selection is set to the ion with higher natural abundance
%
% if the new range covers an existing range, the command gets abborted
%
% if the new range is entirely covered by an existing range, the command
% gets abborted
%
% if the new range partially overlaps with an existing range, the new range
% gets clipped to the adjacent range
%
% if no inserted ion is located within the range, the user can manually enter a
% range name (must be an ion, which is included in colorScheme) in a pop-up window
%
% possible use without specific outputs: rangeAdd(spec,colorScheme)
% generates only area (ans) but no handle or text in the workspace
%
% [h, txt] = rangeAdd(spec,colorScheme)
%
%
% INPUT
% spec:         area plot that displays the mass spectrum (histogram of
%               m/c frequencies)either in raw counts or normalised to bin
%               width and total ion count
%
% colorScheme:  table with elements assigned to color code
%
% manualName:   name of range if no ion is defined, following keywords
%               are possible: 'background'. This defines a background range
%
% OUTPUT
% h:            handle to the area plot of the range
% txt:          corresponding text

%% set current axes
ax = spec.Parent;
axes(ax);

%% user input
lim = ginput(2);

% check how many points were marked
isValid = false;
if ~isempty(lim)
    lim = lim(:,1);
    lim = sort(lim);
    isValid = length(lim) == 2;
end

if isValid
    %% check for manual range name input
    isManual = exist('manualName','var'); % is a manual name present?
    isBackground = false;
    if isManual
        isBackground = strcmp(manualName,'background'); % is it a background range?
    end
    
    %% check for overlap with already existing peak range
    plots = spec.Parent.Children;
    
    idx = 1;
    for pl = 1:length(plots)
        
        % find all the ones that are ranges
        try
            type = plots(pl).UserData.plotType;
        catch
            type = "unknown";
        end
        
        % case #1: new range is entirely part of existing range
        if type == "range"
            if lim(1) > plots(pl).XData(1) & lim(2) < plots(pl).XData(end)
                error 'Total overlap with existing range.'
                return
            end
        end
        
        % case #2: existing range is entirely covered by new range
        if type == "range"
            if lim(1) < plots(pl).XData(1) & lim(2) > plots(pl).XData(end)
                error 'Input range completely covers already existing range.'
            end
        end
        
        % case #3: new range partially overlaps with existing range
        %          new range gets clipped to the adjacent range
        
        % case #3a: partial overlap on right side of existing range
        if type == "range"
            if lim(1) < plots(pl).XData(end) & lim(2) > plots(pl).XData(end)
                lim(1) = plots(pl).XData(end) + (spec.XData(end) - spec.XData(end-1));
                f = msgbox({'Due to a partial overlap,'; 'the new range was clipped to the adjacent range.'},'Notification');
            end
        end
        
        
        % case #3b: partial overlap on left side
        if type == "range"
            if lim(1) < plots(pl).XData(1) & lim(2) > plots(pl).XData(1)
                lim(2) = plots(pl).XData(1) - (spec.XData(end) - spec.XData(end-1));
                f = msgbox({'Due to a partial overlap,'; 'the new range was clipped to the adjacent range.'},'Notification');
            end
        end
    end
    %%
    
    isIn = (spec.XData > lim(1)) & (spec.XData < lim(2));
    h = area(spec.XData(isIn),spec.YData(isIn));
    h.FaceColor = [1 1 1];
    h.UserData.plotType = "range"; %is overwritten for background ranges
    
    % take care of manual ranges
    if isManual & ~isBackground
        [ion, chargeState] = ionConvertName(manualName);
        if ~isnan(chargeState)
            h.UserData.ion = ion;
            h.UserData.chargeState = chargeState;
            h.DisplayName = ionConvertName(h.UserData.ion,h.UserData.chargeState);
            color = colorScheme.color(colorScheme.ion == ionConvertName(h.UserData.ion.element),:);
        else
            h.UserData.ion = manualName;
            h.UserData.chargeState = NaN;
            h.DisplayName = manualName;
            color = colorScheme.color(colorScheme.ion == manualName,:);
        end
        
        if isempty(color)
            name = h.DisplayName;
            delete(h);
            error(['color for atom/ion ' name ' undefined']);
        end
        
        h.FaceColor = color;
        
        
    elseif isManual & isBackground
        h.UserData.plotType = "background";
        h.DisplayName = "background";
        h.FaceColor = colorScheme.color(colorScheme.ion == 'background',:);
    end
    
    
    
    
    %% search for ions in mass spectrum plot
    plots = ax.Children;
    isIon = false(length(plots),1);
    for pl = 1:length(plots)
        try
            isIon(pl) = strcmp(plots(pl).UserData.plotType,"ion");
        end
    end
    ionPlots = plots(isIon);
    
    % find ions in range (if there are any)
    potentialIon = {};
    potentialIonChargeState = [];
    potentialIonPeakHeight = [];
    if ~isempty(ionPlots)
        for pl = 1:length(ionPlots)
            isIn = (ionPlots(pl).XData > lim(1)) & (ionPlots(pl).XData < lim(2));
            if any(isIn)
                % if multiple isotopic combinations of the ion are within the range,
                % the most abundant one is automatically chosen
                isIn = (ionPlots(pl).YData == max(ionPlots(pl).YData(isIn))) & isIn;
                
                potentialIon{end+1} = ionPlots(pl).UserData.ion{isIn};
                if isscalar(ionPlots(pl).UserData.chargeState)
                    potentialIonChargeState(end+1) = ionPlots(pl).UserData.chargeState;
                else
                    potentialIonChargeState(end+1) = ionPlots(pl).UserData.chargeState(isIn);
                end
                potentialIonPeakHeight(end+1) = ionPlots(pl).YData(isIn);
            end
        end
        
    end
    
    %% select which ion it is if necessary
    
    % manual input
    if isempty(potentialIon) & ~isManual
        error('no ion defined in this range: provide name manually');
        
        
        % clear choice
    elseif length(potentialIon) == 1 & ~isManual
        h.UserData.ion = potentialIon{1};
        h.UserData.chargeState = potentialIonChargeState(1);
        h.DisplayName = ionConvertName(h.UserData.ion,h.UserData.chargeState);
        h.FaceColor = colorScheme.color(colorScheme.ion == ionConvertName(h.UserData.ion.element),:);
        
        
    elseif ~isManual % selection
        numPotIon = length(potentialIon);
        for i = 1:numPotIon
            names{i} = [ionConvertName(potentialIon{i}, potentialIonChargeState(i)) '   ' num2str(potentialIonPeakHeight(i))];
        end
        
        % select the ion, defaulting to most abundant
        [~, maxIdx] = max(potentialIonPeakHeight);
        [idx, isSelection] = listdlg('ListString',names,'PromptString','Select ion species','SelectionMode','single',...
            'InitialValue',maxIdx);
        
        if ~isSelection
            delete(h);
            return
        end
        
        h.UserData.ion = potentialIon{idx};
        h.UserData.chargeState = potentialIonChargeState(idx);
        h.DisplayName = ionConvertName(h.UserData.ion,h.UserData.chargeState);
        h.FaceColor = colorScheme.color(colorScheme.ion == h.UserData.ion.element(1),:);
    end
    
    % define for all hit multiplicities
    h.UserData.hitMultiplicities = [0 Inf];
    
    % add text to denote range if it's not a background range
    if ~isBackground
        if isManual
            txt = text(h.XData(1),max(h.YData)*1.4,manualName,'clipping','on');
            txt.DisplayName = manualName;
        else
            txt = text(h.XData(1),max(h.YData)*1.4,ionConvertName(h.UserData.ion,h.UserData.chargeState,'LaTeX'),'clipping','on');
            txt.DisplayName = ionConvertName(h.UserData.ion,h.UserData.chargeState,'plain');
        end
        txt.UserData.plotType = "text";
        
        
        % delete function for ion text and corresponding range
        h.DeleteFcn = @(~,~) delete(txt);
    end
end