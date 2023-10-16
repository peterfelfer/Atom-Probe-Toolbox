function [peakData] = peakBGCorrectedCount(pos, varargin)
% peakBGCorrectedCount fits a linear background (least squares) to a peak in a mass spectrum
% based on 'brushed' background bins. For this, the mass spectrum has to be
% in the active axes, with the background bins marked. Also plots the
% results.
% 
% [peakData] = peakBGCorrectedCount(pos, rangeTable, ionName, boundary)
% 
% [peakData] = peakBGCorrectedCount( pos, rangeTable, ionName)
% [peakData] = peakBGCorrectedCount( pos, rangeTable, ionName, boundary)
% [peakData] = peakBGCorrectedCount( pos, rangeTable, ionName, boundary, options)
%
% INPUT
% massSpec   a massSpec plot - created with massSpecPlot.m
% pos        pos file - in raw format
% ionName    ion Name in this form '14N 12C2 +' it is important to write the
%            specific isotopes
% rangeTable if no range Table is parsed - the gcf has to be the massSpec
%               of the pos variable
% boundary   range in Da before and after the peak that should be used for
%            background correction either one value or two [3 3] are
%            allowed. If no value is parsed, the default value of 0.5 Da is used 
% figure     true - figure window is created (default)
%            false - figure window is suppressed
% offset     if you want to have an offset between the range from the
%            rangeTable and the background correction range just possible
%            with boundaries
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
%% generelle Änderungen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Ich hab das Massenspektrum entfernt. im einfachsten Fall parst
%%%%%%%%%% man nur die pos variable - die funktion springt auf die gcf und
%%%%%%%%%% man kann seine ranges anlegen. Wenn man einen bestimmtne peak
%%%%%%%%%% mit range Table korrigieren mag - dann ist der default wert für
%%%%%%%%%% die bin width vom massenspektrum 0.01 - soll ich den variabel
%%%%%%%%%% machen ? - wenn ja müssen für boundary 2 werte eingegeben werden
%%%%%%%%%% ansonsten kann man es ja nicht unterscheiden was für was ist. 



%% check for RAW pos file
if ismember('atom',pos.Properties.VariableNames)
    pos = posUnDecompose(pos);
end


%% Check for variable input and divide it to ionNameTable, boundary and options

% preset
figOut = 1;
t = 0;
offset = 0;


for k = 1:length(varargin)
    if istable(varargin{k})
        rangeTable = varargin{k};
    elseif ischar(varargin{k}) && length(varargin{k})>1
        [ionNameTable, chargeState] = ionConvertName(varargin{k});
        if isnan(chargeState)
            error('Please add the charge state');
        end
    elseif ~ischar(varargin{k})&& length(varargin{k})>1
        boundary = varargin{k};
        t = 1; % count index - first numbers are for boundary second for offset
    elseif ~ischar(varargin{k})&& length(varargin{k})>1 && t == 1;
        offset = varargin{k};
    elseif ischar(varargin{k}) && length(varargin{k})==1
        options = varargin{k};
    elseif islogical(varargin{k})
        if varargin{k} == 1
            figOut = 1;
        else
            figOut =0;
        end
    end
        
end

%% Check if manual input needed or if variable input with ranges is there

if exist('boundary', 'var') && exist('ionNameTable', 'var')
    % find ranges From range Table
    if istable(ionNameTable)
        for i = 1:height(rangeTable)
         ionRow = ismember(rangeTable.ion{i}, ionNameTable);
         if sum(ionRow) == height(ionNameTable) && chargeState == rangeTable.chargeState(i)
             break
         elseif i == height(rangeTable)
             error('Check if the isotope of the desired ion is correct')
         end
        end
        mcBegin = rangeTable.mcbegin(i);
        mcEnd = rangeTable.mcend(i);      
    end

    if exist('boundary', 'var') && length(boundary) == 1
        boundary = [boundary boundary]; 
    end
    if exist('offset', 'var') && length(offset) == 1
        offset = [offset offset];
    end
    
    % use boundary to get xLim
    if ~exist('boundary', 'var') % no boundary - use default value 1Da
        xLim = [mcBegin-0.5, mcBegin, mcEnd, mcEnd+0.5];
    elseif length(boundary)== 2 % specific background calculation width before and after peak range
       xLim = [mcBegin-offset(1)-boundary(1), mcBegin-offset(1), mcEnd+offset(2), mcEnd+offset(2)+ boundary(2)];
    end
else
    % get baseline input
    % get the region before and after the peak
    [xLim, ~, ~] = ginput(4);
    xLim = sort(xLim); % sort it
end

%% get the X axis of the massSpec
%%%%%%%%%%%%% ich mach das jetzt komplett anders - vorher hat man das
%%%%%%%%%%%%% Massenspektrum gebraucht um den x vektor zu bekommen - der
%%%%%%%%%%%%% war aber allerdings immer in der bin width in der das
%%%%%%%%%%%%% massenspektrum erstellt wurde - ich würde einen default value
%%%%%%%%%%%%% von 0.01 nehmen für alle massenspektren
%%%%%%%%%%%%%         Ich würde auch das massenspektrum komplett rausnehmen
%%%%%%%%%%%%%         und einfach wenn dann nur den rangeTable für die
%%%%%%%%%%%%%         Funktion nehmen



% calculates mcScale x vector for interesting range
mcScale = linspace(xLim(1), xLim(4), round((xLim(4)-xLim(1))/0.01)+1);

%% Calculate lines for background

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
% counts = hist(pos.mc,mcScale); 
counts = histcounts(pos.mc,mcScale);
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

%% calculating atom percentage
% peakcount is diveded by all atoms of whole dataset
numAtoms = max(pos.atomNum);
pctAtoms = pkcnt/numAtoms * 100;
symAtoms = ' %';

if pctAtoms < 0.1
    pctAtoms = pctAtoms/100 * 1E6;
    symAtoms = ' ppm';
end
%% calculating ion percentage
ions = pos.ion;
ions(isundefined(ions))='0';
numIons = sum(ions~='0');
pctIons = pkcnt/numIons * 100;
symIons = ' %';

if pctIons < 0.1
    pctIons = pctIons/100 * 1E6;
    symIons = ' ppm';
end
%% plotting of results

if figOut == 1
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
        ['pct of all atoms: ' num2str(pctAtoms,3) symAtoms],...
        ['pct of all ions: ' num2str(pctIons,3) symIons],...
        ['peak location: ' num2str(pkloc) ' Da']};
    %%%%%%%%%%%%%%%%% ab hier eigentlich nach 'r' und 'a' und vor OUtput
    %%%%%%%%%%%%%%%%% table
    % plotting numerical results
    text(txtPos(1),txtPos(2),txt);
    % exports the figure to the clipboard on Windows
    if ispc
        hgexport(f,'-clipboard');
    end
end

%% executing optional commands
if exist('options','var')

    switch options           

        case 'a'

            % auto range optimisation. moves borders of range until the
            % peak background matches the amount of ions missed.

            % pkloc = peak location
            % cnt = counts per bin
            % mc = mass to charge center of bin

            % divide into before and after peak

            % Counts before Peak 
            % mc enspricht mcScale pkloc - peak location
            
            %% check if fitCts and mcScale have the same size
            if length(mcScale) > length(fitCts)
                mcScale(1)=[];
                mcScale(end) = [];
            end

            
            % Before peak
            ctsBeforePk = counts(mcScale<pkloc);
            mcBeforePk = mcScale(mcScale<pkloc);
            fitCtsBeforePk = fitCts(mcScale<pkloc);
            % Sum before the peak of the normal counts and the fitted
            rngCtsPk = cumsum(ctsBeforePk);
            rngCtsBg = cumsum(fitCtsBeforePk);

           
            BGincluded = rngCtsBg(end) - rngCtsBg;
            missedAtoms = rngCtsPk' - rngCtsBg;


            bal = BGincluded - missedAtoms;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mit Peter
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% absprechen!!!!!!!!!!!!!!!!!!!!!!!!!
            if sum(bal<0) == 0
                bal = BGincluded + missedAtoms;
            end

            mcBegin = min(mcBeforePk(bal<0));



            % mcBefore Pk - sind alle werte bevor der Peak losgeht - also
            
            % die mc die x achse
            % bal = 


            %after
            ctsAfterPk = counts(mcScale>=pkloc);
            ctsAfterPk = fliplr(ctsAfterPk);

            mcAfterPk = mcScale(mcScale>=pkloc);
            mcAfterPk = fliplr(mcAfterPk);

            fitCtsAfterPk = fitCts(mcScale>=pkloc);
            fitCtsAfterPk = fliplr(fitCtsAfterPk);

            rngCtsPk = cumsum(ctsAfterPk);
            rngCtsBg = cumsum(fitCtsAfterPk);

            BGincluded = rngCtsBg(end) - rngCtsBg;
            missedAtoms = rngCtsPk' - rngCtsBg;

            bal = BGincluded - missedAtoms;

            mcEnd = min(mcAfterPk(bal>0));



            % determination of mislabeled atoms
            mcmin = mcBegin;
            mcmax = mcEnd;
            mc = mcScale;
            rngCtsPk = sum(counts((mc>mcmin) & (mc<mcmax)));
            rngCtsBg = sum(fitCts((mc>mcmin) & (mc<mcmax)));

            BGfraction = rngCtsBg/(rngCtsPk+rngCtsBg)*100;

            rngBgGlobal = rngCtsBg/(numAtoms * (mcmax - mcmin));

            peakData.mcbegin = mcBegin;
            peakData.mcend = mcEnd;

            % calculation of the amount of hydrogen of all ions
            ctsPk = rngCtsPk - rngCtsBg;
            corrPct = ctsPk/height(pos) * 100;
            ppm = corrPct/100 * 1E6;

            peakData.countsNew = ctsPk;
            peakData.ppmNew = ppm;
            peakData.pctNew = corrPct;

            sym = ' %';
            if BGfraction < 0.1
                BGfraction = pct/100 * 1E6;
                sym = ' ppm';
            end
           
            if figOut ==1
                txt{end+1} = ' ';
                txt{end+1} = ['range: ' num2str(mcmin) ' - ' num2str(mcmax) ' Da'];
                txt{end+1} = ['range background: ' num2str(BGfraction,3) sym];
                txt{end+1} = ['range background: ' num2str(rngBgGlobal*1E6,3) ' ppm/Da'];
                missedAt = 100 - (rngCtsPk - rngCtsBg)/pkcnt *100;
                txt{end+1} = ['missed atoms: ' num2str(missedAt) ' %'];
                stem([mcmin mcmax],[yLim(2)*rngLabelHeight yLim(2)*rngLabelHeight],'k','Marker','none','LineWidth',1, 'DisplayName', 'optimised range' );
            end
    end
end
%% Output Table
peakData.counts = round(pkcnt);
peakData.pctAtoms = pctAtoms;
peakData.pctIons = pctIons;
peakData.loc = pkloc;

end








