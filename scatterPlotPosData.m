function [p, ax] = scatterPlotPosData(pos,species,sample,colorScheme,size,plotAxis)
% scatterPlotPosData plots APT data in the typical APT style.
% 
% [p, ax] = scatterPlotPosData(pos,species,sample,colorScheme,size,plotAxis)
% [p, ax] = scatterPlotPosData(pos,species,sample,colorScheme,size)
% [p, ax] = scatterPlotPosData(pos,species,sample,colorScheme)
% 
% INPUT
% pos:            decomposed or allocated pos file with pos.ion column      
% 
% species:        Which ion/atom species should be plotted
%                 one single ion {'C'}
%                 more ions {'H','C','N'}
%                 if multiple species are given and no plotAxis, array of 
%                 axes is plotted with synched axis movement. 
%                 Otherwise they all go into the same axis.
% 
% sample:         allows the specification of a subset to plot. 
%                 If sample <1, it is a fraction of the overall number of
%                 atoms, if >1, it is a fixed number. Can be a scalar or a 
%                 vector with same length as species.
% 
% colorScheme:    a colorScheme can be parsed or a RGB color vector
%
% size:           size of the single atoms in the plot
%                 can be a scalar or a vector with same length as species.
%
% plotAxis:       an axis for the plot can be parsed. If no plotAxis is 
%                 parsed, a new axis is created. If you want to plot
%                 mulitple Ions in one plot write axes(figure()) for
%                 plotAxis.
% 
% OUTPUT
% p:              plot handle
%
% ax:             axis from the plot, important for tiled plot
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


%%
% find how many plots need to be done
if iscell(species)
    numPlots = length(species);
else
    numPlots = 1;
    species = {species}; % convert to cell array to be consistant with input of multiple species
end

% default size is 15, if size is not given
if ~exist('size','var')
    size = 15;
end


% check if the pos file is decomposed
isDecomposed = any(string(pos.Properties.VariableNames) == "atom");


%% set up the sample and size variables for multiple plots, if needed also color
if numPlots > 1
    if isscalar(sample); sample = repmat(sample,numPlots,1); end
    if isscalar(size); size = repmat(size,numPlots,1); end
    if ~istable(colorScheme)
        if size(colorScheme,1) == 1
            color = repmat(colorScheme,numPlots,1);
        end
    end
    
end

%% set up plot figure and axis
isTilePlot = false;

if exist('plotAxis','var'); isExistingAxis = true; % so that all plots go in existing axis 
else % plots in new axis/figure
    fig = figure();
    fig.Color = 'w';
    isExistingAxis = false;
    if numPlots > 1
        isTilePlot = true; % if multiple plots are made in a new figure, we use tiles
        tileLayoutHandle = tiledlayout(fig,'flow','TileSpacing','normal');
        fig.Name = 'atom maps';
        
    else
        fig.Name = [species{1} ' atom map'];
    end
end



%% loop through plots
for pl = 1:numPlots % number of species
    
    
    
    % check if color already exist in colorScheme, if not, an error will
    % appear
    iT = ionConvertName(species{1,pl});
    iN = ionConvertName(iT);
    iC = cellstr(iN);
    iCat = categorical(iC);
    colorExist = (colorScheme.ion==iCat);
    if sum(colorExist)==0 
       error('color of the molecule is not in the colorScheme. Please add color to your colorScheme with colorSchemeIonAdd.m')
    end
    
    
    % finding the atoms in the pos file
    % element: 'Fe'
    % elemental isotope: '56Fe'
    % ion: 'Fe+' or 'Fe' if pos variable is not decomposed
    [ionTable, ionChargeState] = ionConvertName(species{pl});
    displayName = ionConvertName(ionTable, ionChargeState, 'LaTeX');
    
    % finds chemical elements
    if any(isnan(ionTable.isotope)) && isDecomposed && isnan(ionChargeState) && height(ionTable) == 1
        % use atoms
        speciesIdx = find(pos.atom == species{pl});
        color = colorScheme.color(colorScheme.ion == species{pl},:);
        
    % finds isotopes (only single atom basis)    
    elseif ~any(isnan(ionTable.isotope)) && isDecomposed && isnan(ionChargeState) 
        if height(ionTable) > 1
            error('plotting of ions based on isotopic information not currently supported');
        end
        speciesIdx = find(pos.atom == ionConvertName(ionTable.element) & pos.isotope == ionTable.isotope);
        color = colorScheme.color(colorScheme.ion == ionConvertName(ionTable.element),:);
        
    elseif ~any(isnan(ionTable.isotope)) && ~isnan(ionChargeState)
        error('plotting of ions based on isotopic information not currently supported');
    
    % finds specific ions with charge state
    elseif ~isnan(ionChargeState)
        speciesIdx = find(pos.ion == ionConvertName(ionTable.element) & pos.chargeState == ionChargeState);
        color = colorScheme.color(colorScheme.ion == ionConvertName(ionTable.element),:);
    % finds specific ions without chargestate
    else
        speciesIdx = find(pos.ion == ionConvertName(ionTable.element));
        if height(ionTable) == 1
            displayName = [displayName ' (ion)'];
        end
        color = colorScheme.color(colorScheme.ion == ionConvertName(ionTable.element),:);
        
    end
    
    
    
    
    
    numAtoms = length(speciesIdx);
    
    % checking if number of atoms/ions requested is greater than number
    % available
    if sample(pl) > numAtoms
        sample(pl) = numAtoms;
        warning(['number of requested ' species{pl}...
            ' atoms/ions greater than number of atoms in dataset. All atoms plotted']);
    elseif sample(pl) <= 1
        sample(pl) = round(numAtoms * sample(pl));
    end
    % calculating sample indices
    plotIndices = speciesIdx(randsample(numAtoms,sample(pl)));
    
    
    % set up axis for plotting if required
    if isExistingAxis
        ax(pl) = plotAxis;
        hold on
    else
        if isTilePlot
            ax(pl) = nexttile(tileLayoutHandle);
        else % new plot in new axis
            ax(pl) = axes(fig);
        end
    end
    axisSpatialAptify;
    
    % displayName = species{pl};
    

    
    % marker setup depending on marker size
    if size(pl) > 35
        edgeColor = [0 0 0];
        fillColor = color;
        markerStyle = 'o';
    else
        edgeColor = color;
        fillColor = color;
        markerStyle = '.';
    end
    
    % actual scatter plot
    p(pl) = scatter3(pos.x(plotIndices),pos.y(plotIndices),pos.z(plotIndices),...
        markerStyle,...
        'MarkerEdgeColor',edgeColor,...
        'MarkerFaceColor',fillColor,...
        'SizeData',size(pl),...
        'DisplayName',displayName,...
        'Parent',ax(pl));
    
    
    if ~isExistingAxis % condition new axes
        ax(pl).Box = 'on';
        ax(pl).BoxStyle = 'full';
        ax(pl).XColor = [1 0 0];
        ax(pl).YColor = [0 1 0];
        ax(pl).ZColor = [0 0 1];
        ax(pl).ZDir = 'reverse';
        if isTilePlot; title(ax(pl),displayName); end
        axis equal;
    end
    
end


%% link rotation between axes
if isTilePlot
    rotLink = linkprop(ax,{'CameraUpVector', 'CameraPosition',...
        'CameraTarget','XLim','YLim','ZLim'}); % links rotation
    setappdata(fig, 'StoreTheLink', rotLink);
end

rotate3d on;
ax = ax';
pl = pl';

