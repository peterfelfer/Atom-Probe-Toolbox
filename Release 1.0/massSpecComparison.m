function specComp = massSpecComparison(varargin)
% massSpecComparison plots multiple massSpec in one figure
% 
% specComp = massSpecComparison(pos1, pos2);
% specComp = massSpecComparison(pos1, pos2, mode);
% specComp = massSpecComparison(pos1, pos2, mode, posNames);
% specComp = massSpecComparison(pos1, pos2, mode, posNames, binWidth, lineSize, colors);
%
%
% INPUT
% mode = 'normalised' or 'count'
% pos = different pos variables up to a number of 7 (because the preset
% colour vector has just 7 colors - if more needed please use the color
% Vector)
% posNames = Struct array with the name of the pos variables {'data1',
% 'data2'}
% bindWidth = number between 0 and 1
% lineSize = number starting from 1
% colors = if specific colors are needed, please parse colors as vectors [0 0.4470 0.7410; 0.8500 0.3250
% 0.0980]
%
% OUTPUT
% specComp = with all the pos variables plotted
% 
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg

%% default

mode = 'normalised';
binWidth = 0.01;
lineSize = 1;
t = 1;

% Standard Matlab Colors
edgeColorVec = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];

%% check input variables 

for k = 1:length(varargin)
    if ischar(varargin{k})
        mode = varargin{k};
    elseif iscell(varargin{k})
        posNames = varargin{k};
    elseif istable(varargin{k})
        posVec{t,1} = varargin{k};
        t = t+1;
    elseif isnumeric(varargin{k}) && length(varargin{k})==1
        if  varargin{k}>=1
            lineSize = varargin{k};
        elseif varargin{k}<1
            binWidth = varargin{k};
        end      
    elseif isnumeric(varargin{k}) && length(varargin{k})==3
        edgeColorVec = varargin{k};
    end
end

% create pos names if variable doesn't exist
if ~exist('posNames', 'var')
    for p = 1:length(posVec)
        posNames(1,p) = {"pos" + p};
    end
end


%% check for Raw pos file input
if length(posVec)<2
    error('There is just one pos Variable parsed')
end

for i = 1:length(posVec)
    pos = posVec{i,1};
    isTableCol = @(pos, atom) ismember(atom, pos.Properties.VariableNames);
    rawCheck = isTableCol(pos, 'atom');
    if rawCheck == 1     
         pos = posUnDecompose(pos);
    end
    posVec{i,1} = pos;
end

%% calculate x and y for pos Variable

for j = 1:length(posVec)
    pos = posVec{j,1};
    mc = pos.mc;
    mcmax = max(mc);

    x = linspace(0,mcmax,round(mcmax/binWidth));

    if  strcmp(mode,'count')
        y = hist(mc,x);
    elseif strcmp(mode,'normalised')
        % calculate as counts/(Da * totalCts) so that mass spectra with different
        % count numbers are comparable
        y = hist(mc,x) / binWidth / length(mc);
        %med = median(y);
    else y = hist(mc,x);
    end
    
    xyVec{j,1} = x;
    xyVec{j,2} = y;
    
end

%% plot figure
specComp = figure('Name','mass spectrum comparison');
ax = axes(specComp);




for h = 1:length(posVec)
    % plot massSpec
    handleIn = area(xyVec{h,1},xyVec{h,2},'FaceAlpha',0, 'EdgeColor', edgeColorVec(h,:), 'EdgeAlpha', 0.5, 'DisplayName', posNames{h}, 'LineWidth', lineSize);
    hold on
end
set(ax,'YScale','Log');

%% legend for figure
legend

ax = get(handleIn,'Parent');
% make axis

set(gcf, 'Color', [1 1 1]);
set(get(gca,'XLabel'),'String','mass-to-charge-state ratio [Da]');

if strcmp(mode,'count')
    ylabel('frequency [counts]');
elseif strcmp(mode,'normalised')
    ylabel('frequency [cts / Da / totCts]');
end


end

