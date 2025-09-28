function graphVMC = voltageCurveMC(epos)
% voltageCurveMC produces a diagram with both the voltage and the
% mass-to-charge values on the left and right y-axis, respectively. On the
% x-axis the ion sequence number is plotted.
% The mass-to-charge values are clipped to a reasonable range starting from
% 0 Da.
%
% INPUT
% epos:         epos file, table with experiment data
%
% OUTPUT
% graphVMC:     figure, graph with voltage (V) and mass-to-charge ratio
%               (MC) plotted over the ion sequence number
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg


%% create diagram with two graphs (base voltage curve and mc ratio)
% create new figure
graphVMC = figure;
xlabel('ion index [-]');

% define x and both y values (y1 and y2)
x = epos.ionIdx;
y1 = epos.VDC;
y2 = epos.mc;

% right y-axis: plot mc ratio values over ion index number
yyaxis right;

    % plot mc histogram so that the user can select the desired maximum mc value
    hist(epos.mc,height(epos)/4000);
    set(gca, 'YScale', 'log');
    limits = ginput(1);
    ylim([0,limits(1)]);

scatter(x(1:500:end),y2(1:500:end),1,'.');
ylim([0,limits(1)]); % clip right y axis to previously determined maximum mc value
ylabel('mass-to-charge-ratio [Da]');

% left y-axis: plot base voltage values over ion index number
yyaxis left;
plot(x(1:5000:end),y1(1:5000:end)); % line plot, only every 5000th point is considered 
ylabel('base voltage [V]');
xlim([0,max(x)]);

end

