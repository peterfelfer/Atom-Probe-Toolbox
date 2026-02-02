%[text] # **colorScheme**
%[text] **creation of ion table and colorScheme, addition of ion and corresponding color to existing colorScheme** 
%[text] The *colorScheme* is a table with a three element color vector for each ion or molecule.This list holds all of the elements of the periodic table and a few molecules. However, the list makes no claim to be exhaustive. If the user detects new molecules, the user can add these molecules to the list and send a note to us. We will expand the list with more ions in the future.
%[text] In addition, the user can create and modify the *colorScheme*. It is possible to create a complete *colorScheme* for specific data sets or material classes.
%%
%[text] ## Creation of a new colorScheme 
%[text] With the function *colorSchemeCreate* it is possible to create a new color scheme with colors adjusted to the user's data set. The colors will be created in a first step as HSV codes (with fixed values of S = 0.8 and V = 1.0) with equidistant hue values (between 0 and 0.8, corresponding to degree values of 0 to about 290°). The created colors will then be converted into RGB codes.
%[text] The created color scheme, however, may not be ideal since the human color perception is not linear. Therefore, subsequent changes within the *colorScheme* may be necessary.
%[text] To create a new color scheme, the user needs an *ionTable*, which can be created in two ways. 
%[text] One way is to extract the ions from the mass spectrum that contains all of the possible ions.
% ionTable = ionsExtractFromMassSpec(spec);
%[text] The second way is to create a new *ionTable* based on a list of elements. Here, the user hast to put in all wanted ion symbols in the way {'Fe'; 'C'; 'Mn'}.
ions = {'Fe'; 'C'; 'H'; 'N';'O'}; % cell array of all the possible ions of the data set %[control:editfield:3e24]{"position":[8,33]}
ionTable = table(categorical(ions));
ionTable.Properties.VariableNames{1} = 'ion';
  %[control:button:70c9]{"position":[1,2]}
%%
%[text] With this generated *ionTable* it is possible to create a custom *colorScheme*.
colorScheme_new = colorSchemeCreate(ionTable)
  %[control:button:7b00]{"position":[1,2]}
%%
%[text] ## Add a new ion to the colorScheme
%[text] It is possible to add a new (molecular) ion with a new color to the *colorScheme*. By choosing *selection = 'select'* the user can choose a specific color for the ion or molecule. With *selection = 'create'* the function creates a new color by picking the one that is the farthest from all already existing colors.
% load colorScheme.mat;
newIon = '1865Da'; % name of the new molecule %[control:editfield:2bc7]{"position":[10,18]}
selection = 'create'; % selection mode %[control:dropdown:6328]{"position":[13,21]}
colorScheme = colorSchemeIonAdd(colorScheme, newIon, selection) %[output:28a9f542]
  %[control:button:148f]{"position":[1,2]}
%%
%[text] ## Adding new ions to the colorScheme with the ionTable
%[text] When a lot of new ions are added to the mass spectra that are not in the colorScheme, the following code can be executed to complete the colorScheme. The output "ion already exists in colorScheme" can be ignored
ionTable = ionsExtractFromMassSpec(spec); % create the ionTable with all the detected ions %[output:4622aff2]
newName = char(ionTable.ionName); % get the names and convert them into characters
for i = 1:size(newName)
    colorScheme = colorSchemeIonAdd(colorScheme, newName(i,:)); % add each of them to the colorScheme
end 
      %[control:button:9352]{"position":[5,6]}

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":8.7}
%---
%[control:editfield:3e24]
%   data: {"defaultValue":"{'Fe'; 'C'; 'H'; 'N';'O'}","label":"ions","run":"Nothing","valueType":"MATLAB code"}
%---
%[control:button:70c9]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:7b00]
%   data: {"label":"Run","run":"Section"}
%---
%[control:editfield:2bc7]
%   data: {"defaultValue":"'1865Da'","label":"Edit field","run":"Nothing","valueType":"Char"}
%---
%[control:dropdown:6328]
%   data: {"defaultValue":"'create'","itemLabels":["'select'","'create'"],"items":["'select'","'create'"],"label":"selection","run":"Nothing"}
%---
%[control:button:148f]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:9352]
%   data: {"label":"Run","run":"Section"}
%---
%[output:28a9f542]
%   data: {"dataType":"tabular","outputData":{"columnNames":["ion","color"],"columns":2,"dataTypes":["categorical","double"],"groupedColumnIndices":[null,null],"header":"146×2 table","name":"colorScheme","rows":146,"type":"table","value":[[["H"],["0.6981","0.6665","0.1781"]],[["He"],["0.1280","0.9991","0.1711"]],[["Li"],["0.0326","0.5612","0.8819"]],[["Be"],["0.6692","0.1904","0.3689"]],[["B"],["0.4607","0.9816","0.1564"]],[["C"],["0.8555","0.6448","0.3763"]],[["N"],["0.1909","0.4283","0.4820"]],[["O"],["0.1206","0.5895","0.2262"]],[["F"],["0.3846","0.5830","0.2518"]],[["Ne"],["0.2904","0.6171","0.2653"]],[["Na"],["0.8244","0.9827","0.7302"]],[["Mg"],["0.3439","0.5841","0.1078"]],[["Al"],["0.9063","0.8797","0.8178"]],[["Si"],["0.2607","0.5944","0.0225"]]]}}
%---
%[output:4622aff2]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Array indices must be positive integers or logical values.\n\nError in <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('ionConvertName>sortedFindgroups', 'Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionConvertName.m', 294)\" style=\"font-weight:bold\">ionConvertName>sortedFindgroups<\/a> (<a href=\"matlab: opentoline('Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionConvertName.m',294,0)\">line 294<\/a>)\n    isotopeGroup(i) = idx(isotopeGroup(i));\n\nError in <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('ionConvertName', 'Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionConvertName.m', 181)\" style=\"font-weight:bold\">ionConvertName<\/a> (<a href=\"matlab: opentoline('Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionConvertName.m',181,0)\">line 181<\/a>)\n    isotopeGroup = sortedFindgroups(ionTable);\n\nError in <a href=\"matlab:matlab.internal.language.introspective.errorDocCallback('ionsExtractFromMassSpec', 'Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionsExtractFromMassSpec.m', 38)\" style=\"font-weight:bold\">ionsExtractFromMassSpec<\/a> (<a href=\"matlab: opentoline('Z:\\GitHub-Matlab-Atomsonde\\Atom-Probe-Toolbox\\ionsExtractFromMassSpec.m',38,0)\">line 38<\/a>)\n        ionName{idx,:} = ionConvertName(plots(pl).UserData.ion{1}.element);"}}
%---
