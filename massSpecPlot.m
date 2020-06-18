function handle = massSpecPlot(mc, bin, mode)
% plotMassSpec plots the data from pos to get a Massspectrum
%
% handle = plotMassSpec(mc, bin, mode)
% handle = plotMassSpec(mc, bin)
%
% INPUTS
% mc:       is the mass-to-charge(mc)-ratio [Da] of the events in the
%           APT-Measurement stored in pos, table
% bin:      is the width of the steps in which the plot is performed,
% mode:     Specifies the way the counts are applied
%           'count' records the number of counts
%           'normalised' records the number of counts if the bin was one Da
%           wide over the overall number of counts
%           default mode is 'count'
%
% OUTPUTS
%           handle to the plot that contains counts or
%           (counts/Dalton)/totalCounts over Dalton. Used in further
%           analysis to find new ions


% count is default mode
if ~exist('mode','var')
    mode = 'count';
end

if istable(mc)
    mc = mc.mc;
end

if length(mc(1,:)) > 1
    mc = mc(:,4);
end


f = figure('Name','mass spectrum');
ax = axes(f);

mcmax = max(mc);

x = linspace(0,mcmax,round(mcmax/bin));

if  strcmp(mode,'count')
    y = hist(mc,x);
elseif strcmp(mode,'normalised')
    % calculate as counts/(Da * totalCts) so that mass spectra with different
    % count numbers are comparable
    y = hist(mc,x) / bin / length(mc);
    %med = median(y);
else y = hist(mc,x);
end

% plot all mass spectrum
handle = area(x,y,'FaceColor',[.9 .9 .9]);
handle.UserData.plotType = "massSpectrum";
hold on;
ax = get(handle,'Parent');

set(gca,'YScale','Log');
set(gcf, 'Name', 'Mass spectrum');
set(gcf, 'Color', [1 1 1]);
set(get(gca,'XLabel'),'String','mass-to-chargestate [Da]');

if strcmp(mode,'count')
    ylabel('frequency [counts]');
elseif strcmp(mode,'normalised')
    ylabel('frequency [cts / Da / totCts]');
end


%% annotation with range stats
t = annotation('textbox');
% determining the background at 4Da
upperLim = 4.5; %Da
lowerLim = 3.5; %Da
BG4 = sum(y((x >= lowerLim) & (x <= upperLim)))/(upperLim-lowerLim);
BG4 = BG4/length(mc) * 1E6;
t.String = {['bin width: ' num2str(bin) ' Da'], ['num atoms: ' num2str(length(mc)) ], ['backG @ 4Da: ' num2str(BG4,3) ' ppm/Da']};
t.BackgroundColor = 'w';
t.FaceAlpha = 0.8;
t.Position = [.15 .8 .27 .1];
pan xon
zoom xon

handle.DisplayName = 'mass spectrum';
