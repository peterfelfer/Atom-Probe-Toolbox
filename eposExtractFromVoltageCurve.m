function eposExtract = eposExtractFromVoltageCurve(epos)
% dataExtractFromVoltageCurve extracts data on the basis of a voltage curve
% 
% INPUT
% epos:        epos file, table with experiment data 
%
% OUTPUT
% eposExtract: table with experiment data, subset of initial epos file,
%              clipped on the basis of the user's selection
%

%%   plot voltage curve
figure;
plot(epos.ionIdx, epos.VDC);
title('voltage curve');
xlabel('ion index');
ylabel('base voltage [V]');

%%   select start and end value of desired region of voltage curve
voltageValues = ginput(2);
voltageValueStart = voltageValues(1);               % for extraction not the voltage values (y coordinates of the points)
voltageValueEnd = voltageValues(2);                 % but the x coordinates of the points are used

%%   generate clipped epos table
eposExtract = epos((epos.ionIdx >= voltageValueStart & epos.ionIdx <= voltageValueEnd),:);

end

