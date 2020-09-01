function eposExtract = eposFromFDM(epos)
% eposFromFDM lets the user define an area on the field desorption map
% (FDM) and extracts the data and generates a new epos (eposExtract)
%
% INPUT
% epos:             epos file, table with experiment data
%
% OUTPUT
% eposExtract:      table with experiment data, subset of initial epos file,
%                   clipped on the basis of the user's selection

%   create rudimentary field desorption map (FDM) (with IMELLIPSE)
FDM = histogram2(epos.detx,epos.dety,100,'FaceColor','flat');
axis equal;
view(90,90);
camroll(90);
xlim([-20,20])
ylim([-20,25])
set(gca,'ZDir','reverse');
xlabel('x-coordinate of detector');
ylabel('y-coordinate of detector');

%   user input of desired area on FDM (ellipse)
annotation('textbox',[0.29, 0.9, 0, 0],'String','Press any key if selection is completed','FitBoxToText','on');
h = imellipse;
pause;

% convert ellipse vertices into table with columns x and y
vert = h.getVertices;
vert = array2table(vert);
vert.Properties.VariableNames{1} = 'x';
vert.Properties.VariableNames{2} = 'y';

%   create new epos with data of selected area on FDM
eposExtract = epos(inpolygon(epos.detx, epos.dety, vert.x, vert.y),:);

% % with DRAWELLIPSE
% FDM = histogram2(epos.detx,epos.dety,100,'FaceColor','flat');
% axis equal;
% view(90,90);
% camroll(90);
% xlim([-20,20])
% ylim([-20,25])
% set(gca,'ZDir','reverse');
% xlabel('x-coordinate of detector');
% ylabel('y-coordinate of detector');
% 
% %   user input of desired area on FDM (ellipse)
% annotation('textbox',[0.29, 0.9, 0, 0],'String','Press any key if selection is completed','FitBoxToText','on');% 
% h = drawellipse;
% pause;
%
% % convert ellipse vertices into table with columns x and y
% vert = h.Vertices;
% vert = array2table(vert);
% vert.Properties.VariableNames{1} = 'x';
% vert.Properties.VariableNames{2} = 'y';
% 
% %   create new epos with data of selected area on FDM
% eposExtract = epos(inpolygon(epos.detx, epos.dety, vert.x, vert.y),:);

close(gcf);
end

