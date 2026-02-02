function h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme, isTracer)
% ionStemPlot plots the relative abundances of an ion in a stem plot
% the information about the ion is given in h.UserData
%
% h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme)
% h = ionStemPlot(ax, weight, abundance, ionList, chargeStates, colorScheme, isTracer)
%
% INPUT
% ax:           axis of the current mass spectrum
%
% weight:       weight of the ion in amu
%
% adundance:    adundance for the chosen ion
%
% ionList:      list of all possible ions
%
% chargeStates: charge states for the given ions
%
% colorScheme:  color scheme as provided or self made
%
% isTracer:     (optional) logical flag indicating if this is a tracer ion
%               Default: false. If true, uses diamond markers instead of
%               circles and appends " (tracer)" to display name.
%
% OUTPUT
% h:            handle to the stem plot
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

% Handle optional isTracer argument
if ~exist('isTracer', 'var') || isempty(isTracer)
    isTracer = false;
end

% generate ion name
ion = ionConvertName(ionList{1}.element);

h = stem(ax,weight, abundance);
try
    h.Color = colorScheme.color(colorScheme.ion == ion,:);
catch
    warning('ion color undefined');
end

if length(chargeStates) == 1
    displayName = [ion repmat('+',1,chargeStates)];
else
    displayName = ion;
end

% Add tracer suffix to display name if applicable
if isTracer
    displayName = [displayName ' (tracer)'];
end

h.DisplayName = displayName;
h.LineWidth = 2;
h.UserData.plotType = "ion";
h.UserData.ion = ionList;
h.UserData.chargeState = chargeStates;
h.UserData.isTracer = isTracer;  % Store tracer flag for downstream use
h.ButtonDownFcn = @(~,~) disp([h.DisplayName ' RGB color: ' num2str(h.Color)]);

% Set marker style: diamond for tracers, circle for regular ions
if isTracer
    h.Marker = 'd';  % Diamond marker for tracers
    h.MarkerSize = 8;
    h.MarkerFaceColor = h.Color;  % Filled marker for visibility
else
    h.Marker = 'o';  % Circle marker for regular ions
end

% change stem line depending on charge state if only one charge state is given
if length(chargeStates) == 1
    switch chargeStates
        case 1
            h.LineStyle = '--';
        case 2
            h.LineStyle = ':';
        case 3
            h.LineStyle = '-.';
        case 4
            h.LineStyle = '-';

    end
end
