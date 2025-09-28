% paraprobe-transcoder
% script to read Cameca/AMETEK *.APT file from APSuite6/IVAS4
% Markus K\"uhbach, m.kuehbach@mpie.de, 2020/03/19

%clearvars;

%specify the location of the *.apt input file
%%fn = '\\maaptsrv.mpie.de\Atomprobe\LEAP5076XS\data_analyse\Felipe Oliveira\Cameca\None\Reconstructions\b4980dbf-5574-4bfb-8842-04ce801dbd59.apt';

%create a class object containing all variables of the *.apt file
apt = PARAPROBE_Transcoder2(fn); 
dx = apt.XDet_mm;
dy = apt.YDet_mm;
x = apt.Position(1,:);
y = apt.Position(2,:);
z = apt.Position(3,:);
m = apt.Mass;
t = apt.tof;
vdc = apt.Vref;
vp = apt.Vap;
nat_pulse = apt.Multiplicity;

% exemplarily call a method of the internal classes inside the file object
% here to show what the individual section of the file contain
apt.header.print();

% Example

dat_ions=[dx; dy]';

%Apply routine, need to download kde2d -> https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
[bandwidth,density,X,Y] = kde2d(dat_ions);
%Plot the data and the density estimate
figure()
surf(X,Y,density,'LineStyle','none'), view([0,90])
colormap hot, hold on, alpha(.8)
set(gca, 'color', 'blue');
axis equal
axis([min(apt.XDet_mm) max(apt.XDet_mm) min(apt.YDet_mm) max(apt.YDet_mm)])

