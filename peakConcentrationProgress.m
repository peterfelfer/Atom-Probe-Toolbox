function p = peakConcentrationProgress(epos,bin,startPeak,endPeak)
% peakConcentrationProgress plots the evolution of the concentration of
% the assigned mc range over the progressing measurement.

% The bin specifies the number of hits taken into account for each
% concentration calculation.
%
%

%% asdg
eposMax = height(epos);
x = linspace(0,eposMax,round(eposMax/bin)-2);

for i = 1:((round(height(epos)/bin))-2);
    start = bin*i;
    ending = bin*i+bin;
    eposPart = epos(start:ending,:);
    contentPeak = sum(eposPart.mc >= startPeak & eposPart.mc <= endPeak);
    contentRest = sum(eposPart.mc >= 0 & eposPart.mc < startPeak) + sum(eposPart.mc > endPeak);
    conc(i) = (contentPeak / contentRest)*100;
end

figure;
p = area(x,conc,'FaceColor',[.9 .9 .9]);
xlabel(['ion hit index / -']);
ylabel(['concentration in the range of ',num2str(startPeak) ' and ' num2str(endPeak),' Da / at.%']);
hold on;
yyaxis right
plot(epos.ionIdx,epos.VDC);
ylabel(['base voltage / V'])