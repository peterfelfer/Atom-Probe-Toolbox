function colorScheme = colorSchemeCreate(ionTable)
% colorSchemeCreate creates a table with RGB color codes for each ion; 
% the colors are generated in a way, that no two colors are perceived as similar
%
% INPUT: 
% ionTable:     table with variable ion (ionTable.ionName)
%               example: ionTable = table(categorical({'Fe'; 'C'; 'H'; 'N'}));
%                        ionTable.Properties.VariableNames{1} = 'ionName';
%
% OUTPUT:
% colorScheme:  table with ion field and corresponding color field
 
%% Check for ions that are twice in the ionTable because of their charge state

Tnew = unique(ionTable.ionName);
ionTable = table(Tnew);
ionTable.Properties.VariableNames{1} = 'ionName';

%% determine value for Value (V) of the HSV color code
V = 1.0;                % Value, value between 0 and 1 [1 means 100 %]
 
%% generate evenly spaced values for H
lsH = table(linspace(0, 0.8, height(ionTable)));     % value of 0 is equivalent to value of 1, therefore upper limit is set to 0.8
lsH.Properties.VariableNames{1} = 'Hue';
 
%% generate evenly spaced values for S
lsS = table(linspace(0.2, 1, height(ionTable)));     % start value set to 0.2 so the color 'white' cannot be generated
lsS.Properties.VariableNames{1} = 'Saturation';
 
%% generate colorScheme table
colorScheme = table('Size',[height(ionTable) 2],'VariableTypes',{'categorical','double'});
colorScheme.Properties.VariableNames{1} = 'ion';
colorScheme.Properties.VariableNames{2} = 'color';
 
%% fill colorScheme table with HSV values
for i = 1:height(ionTable);
    colorScheme.ion(i) = ionTable.ionName(i);
    colorScheme.color(i,1) = lsH.Hue(1,i);
    colorScheme.color(i,2) = lsS.Saturation(1,i);
    colorScheme.color(i,3) = V;
end
 
%% create colorScheme table with RGB values
 cst = table(colorScheme.color);
 HSV = table2array(cst);           % insert HSV values in array
 
 % convert HSV color code to RGB color code
 RGB = hsv2rgb(HSV);               
 
 colorScheme.color = RGB;
end
