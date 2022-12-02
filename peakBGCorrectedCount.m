function [peakData] = peakBGCorrectedCount(massSpec, pos, options)
% peakBGCorrectedCount fits a linear background (least squares) to a peak in a mass spectrum
% based on 'brushed' background bins. For this, the mass spectrum has to be
% in the active axes, with the background bins marked. Also plots the
% results.
% I
% INPUT
% massSpec   a massSpec plot - created with massSpecPlot.m
% pos        pos file - in raw format
% options    'r': add a fictional range and calculate the signal to background for that
%            range
%            'a': same as 'r', but automatically created ranges by optimising the
%            ratio of missed atoms / background atoms = 1
%
% OUTPUT
% peakData   struct array with the counts of the ions in the peak, the
%            percentage and the location of the peak
% figure     figure with the peak and the correction


rngLabelHeight = 0.65; % height of the stem plot delineating the range



%% get the X axis of the massSpec
if exist('massSpec','var')
    handles = get(gca,'Children');
    mcScale = get(handles(end),'Xdata'); % spec X achse - als vektor
else
    error("please create a massSpectra with the function massSpecPlot")
end
%% check for RAW pos file
isTableCol = @(pos, atom) ismember(atom, pos.Properties.VariableNames);
rawCheck = isTableCol(pos, 'atom');
if rawCheck == 1     
     error("please use a raw pos variable")
end

%% get baseline input

% get the region before and after the peak
    [xLim, ~, ~] = ginput(4);
    xLim = sort(xLim); % sort it
% Line before peak
    inBefore = mcScale>xLim(1) & mcScale<xLim(2);
% Line after peak
    inAfter = mcScale>xLim(3) & mcScale<xLim(4);

% Creates a logical - 1 - in range for background correction - 0 not in
% range for the entire dataset
in = or(inBefore,inAfter);

% just the range for the background correction
    mcData = mcScale(in);
%% get the counts per bin    
% neues Count - wie viele ionen hab ich pro bin?
counts = hist(pos.mc,mcScale); % for the entire dataset
%brauch ich vllt nicht wegen entire range 
cntData = counts(in); % just select the bins that are in the range 
% calculation Range - range based for the calculation with the adjacent
% counts of each bin 
calcRange = table(mcData', cntData', 'VariableNames', {'mcRange' 'counts'});

% find the start and end point of the range for background correction of the entire data set 
idxBeg = find(mcScale == min(calcRange.mcRange));
idxEnd = find(mcScale == max(calcRange.mcRange));

% get the entire range for background correction + peak in between both
entireRange = table ((mcScale(idxBeg:idxEnd))', (counts(idxBeg:idxEnd))', in(idxBeg:idxEnd)', 'VariableNames', {'mcRange', 'counts', 'corrRange'});


% Peak location - find Peak position
pkloc = entireRange.mcRange(entireRange.counts == max(entireRange.counts));


%% calculate linear fit of baseline
lin_reg = polyfit(entireRange.mcRange(entireRange.corrRange == 1),entireRange.counts(entireRange.corrRange == 1),1);

a1 = lin_reg(1);
b1 = lin_reg(2);

fitCts = a1*entireRange.mcRange+b1;

sumFitCnt = sum(fitCts);
sumCnt = sum(entireRange.counts);

pkcnt = sumCnt - sumFitCnt;

%% plotting of results
numAtoms = sum(counts);
pct = pkcnt/numAtoms * 100;
sym = ' %';


if pct < 0.1
    pct = pct/100 * 1E6;
    sym = ' ppm';
end


f = figure();
plot(entireRange.mcRange,entireRange.counts);
hold on
plot(entireRange.mcRange,fitCts,'r','LineWidth',2);


set(gcf, 'Color', [1 1 1]);
set(gcf, 'Name', 'Background fit for peak');
set(get(gca,'XLabel'),'String','mass-to-chargestate [Da]');
set(get(gca,'YLabel'),'String','frequency');

legend('mc spectrum','bg fit');
xLim = get(gca,'XLim');
yLim = get(gca,'YLim');
txtPos = [xLim(1) + 0.02 * (xLim(2) - xLim(1)); ...
    yLim(1) + 0.8 * (yLim(2) - yLim(1))];
txt = {['ions in peak: ' num2str(round(pkcnt))],...
    ['pct of all ions: ' num2str(pct,3) sym],...
    ['peak location: ' num2str(pkloc) ' Da']};


%% executing optional commands
if exist('options','var')
    
    switch options
        case 'r'
            
            %signal/background for range
            rec = getrect(gca);
            
            mcmin = rec(1);
            mcmax = mcmin + rec(3);
            
            rngCtsPk = sum(cnt((mc>mcmin) & (mc<mcmax)));
            rngCtsBg = sum(fitCts((mc>mcmin) & (mc<mcmax)));
            
            BGfraction = rngCtsBg/(rngCtsPk+rngCtsBg)*100;
            
            rngBgGlobal = rngCtsBg/(numAtoms * (mcmax - mcmin));
            
            sym = ' %';
            if BGfraction < 0.1
                BGfraction = pct/100 * 1E6;
                sym = ' ppm';
            end
            
            txt{end+1} = ' ';
            txt{end+1} = ['range: ' num2str(mcmin) ' - ' num2str(mcmax) ' Da'];
            txt{end+1} = ['range background: ' num2str(BGfraction,3) sym];
            txt{end+1} = ['range background: ' num2str(rngBgGlobal*1E6,3) ' ppm/Da'];
            
            
            missedAt = 100 - (rngCtsPk - rngCtsBg)/pkcnt *100;
            txt{end+1} = ['missed atoms: ' num2str(missedAt) ' %'];
            
            
            stem([mcmin mcmax],[yLim(2)*rngLabelHeight yLim(2)*rngLabelHeight],'k','Marker','none','LineWidth',2);
            
            
        case 'a'
            
            % auto range optimisation. moves borders of range until the
            % peak background matches the amount of ions missed.
            
            % pkloc = peak location
            % cnt = counts per bin
            % mc = mass to charge center of bin
            
            % divide into before and after peak
            
            % before
            ctsBeforePk = cnt(mc<pkloc);
            
            mcBeforePk = mc(mc<pkloc);
            fitCtsBeforePk = fitCts(mc<pkloc);
            
            rngCtsPk = cumsum(ctsBeforePk);
            rangeCtsBg = cumsum(fitCtsBeforePk);
            
            BGincluded = rangeCtsBg(end) - rangeCtsBg;
            missedAtoms = rngCtsPk - rangeCtsBg;
            
            bal = BGincluded - missedAtoms;
            
            mcBegin = min(mcBeforePk(bal<0));
            
            
            
            
            %after
            ctsAfterPk = cnt(mc>=pkloc);
            ctsAfterPk = fliplr(ctsAfterPk);
            
            mcAfterPk = mc(mc>=pkloc);
            mcAfterPk = fliplr(mcAfterPk);
            
            fitCtsAfterPk = fitCts(mc>=pkloc);
            fitCtsAfterPk = fliplr(fitCtsAfterPk);
            
            rngCtsPk = cumsum(ctsAfterPk);
            rangeCtsBg = cumsum(fitCtsAfterPk);
            
            BGincluded = rangeCtsBg(end) - rangeCtsBg;
            missedAtoms = rngCtsPk - rangeCtsBg;
            
            bal = BGincluded - missedAtoms;
            
            mcEnd = min(mcAfterPk(bal>0));
            
            
            
            % determination of mislabeled atoms
            mcmin = mcBegin;
            mcmax = mcEnd;
            
            rngCtsPk = sum(cnt((mc>mcmin) & (mc<mcmax)));
            rngCtsBg = sum(fitCts((mc>mcmin) & (mc<mcmax)));
            
            BGfraction = rngCtsBg/(rngCtsPk+rngCtsBg)*100;
            
            rngBgGlobal = rngCtsBg/(numAtoms * (mcmax - mcmin));
            
            peakData.mcbegin = mcBegin;
            peakData.mcend = mcEnd;
            
            sym = ' %';
            if BGfraction < 0.1
                BGfraction = pct/100 * 1E6;
                sym = ' ppm';
            end
            
            txt{end+1} = ' ';
            txt{end+1} = ['range: ' num2str(mcmin) ' - ' num2str(mcmax) ' Da'];
            txt{end+1} = ['range background: ' num2str(BGfraction,3) sym];
            txt{end+1} = ['range background: ' num2str(rngBgGlobal*1E6,3) ' ppm/Da'];
            
            
            
            missedAt = 100 - (rngCtsPk - rngCtsBg)/pkcnt *100;
            txt{end+1} = ['missed atoms: ' num2str(missedAt) ' %'];
            
            
            stem([mcmin mcmax],[yLim(2)*rngLabelHeight yLim(2)*rngLabelHeight],'k','Marker','none','LineWidth',2);
    end
end



%% plotting numerical results

text(txtPos(1),txtPos(2),txt);


peakData.counts = round(pkcnt);
peakData.pct = pct;
peakData.loc = pkloc;


% exports the figure to the clipboard on Windows
if ispc
    hgexport(f,'-clipboard');
end

end