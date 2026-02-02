%[text] # Analysis of metadata in APT datasets
%[text] If the user has stored all of his APT measurments in seperate folders containing the metadata list, the .epos or .pos file and/or a mass spectra with a ranged figure, Data analysis over all measurements is possible. 
%%
%[text] ## Find all folders
%[text] First of all, the user has to find all of the APT data folders in the folder structure. If there are also folders that do not contain data or that does not interest you, you may have to remove these folders. Just have a look at the variable 'folders' after running the section. The variable folder contain the name, the stored folder, the date and bytes ... ??
%path = '/Volumes/apt-data';
path = 'L:\';
% extract all subfolders of the path
contents = dir(path);
folders = contents([contents(:).isdir]);
folders = folders(~ismember({folders(:).name},{'.','..'})); % remove hidden folders
  %[control:button:8406]{"position":[1,2]}
%%
%[text] ## Load metadata for all folders
%[text] All of your metadata file information of all of your datasets is stored in one variable called metaDataList
numDatasets = length(folders);

metaDataList = {};
for ds = 1:numDatasets %[output:group:8f271d87]
    try
        disp(folders(ds).name); %[output:68d230e2] %[output:5306653f] %[output:9481c78c] %[output:55af6516] %[output:48c7ba05] %[output:28f0fe4f] %[output:4efdb7e9] %[output:940b65f0] %[output:0caa6970] %[output:76ea124a] %[output:3d0e6b0a] %[output:13680c1a] %[output:954ccde1] %[output:1ffde8bd] %[output:59b956c1] %[output:7692b052] %[output:5599cfe7] %[output:8dfcac42] %[output:88fb106a] %[output:6990bac3] %[output:321fd6c4] %[output:8f63d0ac] %[output:00ac8e4f] %[output:92245aa0] %[output:569e681a] %[output:0653b743] %[output:350ebc2e] %[output:1b8f4ec8] %[output:06002130] %[output:8da7af50] %[output:1be49a22] %[output:75d5136e] %[output:646094a1] %[output:2564afd7] %[output:91d0c8ee] %[output:89dfb9c8] %[output:3ce91317] %[output:7a713631] %[output:9d4f8109] %[output:9b4080db] %[output:16ce2af9] %[output:62552850] %[output:2b60dc7b] %[output:63e5a11b] %[output:96aaa77a] %[output:26d45876] %[output:3a142d02] %[output:654d5778] %[output:249f44b0] %[output:5b945e1b] %[output:5f8f1b03] %[output:24d20803] %[output:370a26c2] %[output:9b7267bb] %[output:7f9c8c28] %[output:78520967] %[output:313ea599] %[output:8ae247e2] %[output:4d2abc5c] %[output:5bf4dd0c] %[output:038e74cb] %[output:9b515301] %[output:5fdb23d5] %[output:3ef9abad] %[output:5638629e] %[output:8bed0e7f] %[output:8e596f47] %[output:5e6b1012] %[output:8ac17693] %[output:779c9df5] %[output:1de604e1] %[output:27f719ba] %[output:3a560479] %[output:5255a2ef] %[output:4fda71c6] %[output:4e834c5b] %[output:6d5af7ae] %[output:32b33c27] %[output:7b6252c0] %[output:0bc197ae] %[output:3eae81dc] %[output:7b6a6c82] %[output:1c89fb78] %[output:03ecb636] %[output:1623a0d0] %[output:1f49a2e0] %[output:29213785] %[output:49e5665e] %[output:161f109d] %[output:63837842] %[output:82e0445a] %[output:9d171763] %[output:551918c7] %[output:0ca43578] %[output:579fd785] %[output:4c5bcb03] %[output:4886bb3f] %[output:1356f23f] %[output:01ff7d9c] %[output:74847e87] %[output:7ce5af0d] %[output:42933222] %[output:29f0aa33] %[output:103e1fc9] %[output:279abdfc] %[output:9f9d5e2d] %[output:416594ed] %[output:43d464e5] %[output:06c22b32] %[output:35283f5e] %[output:03d9da65] %[output:7fd99ec1] %[output:2308259a] %[output:0df3c767] %[output:6fe13e16] %[output:7c91d561] %[output:265bbd92] %[output:72bc1e81] %[output:7639c7a8] %[output:2132a507] %[output:807e3719] %[output:7d983c8f] %[output:03a32c0c] %[output:02151087] %[output:5230dd18] %[output:9fdc90e2] %[output:74d6da98] %[output:47cb410f] %[output:11db43f0] %[output:7104c10f] %[output:07301523] %[output:565158d4] %[output:8127d8c4] %[output:2bfd026b] %[output:3e1ffa50] %[output:39ab5672] %[output:27592aeb] %[output:58b8a750] %[output:36053a74] %[output:29e06349] %[output:2cab0ae8] %[output:26ca7546] %[output:6eccdf08] %[output:500a5c1e] %[output:98df6923] %[output:1810de4b] %[output:4229fb09] %[output:3928dc91] %[output:3dbdf65f] %[output:1624d433] %[output:674ccd0f] %[output:98724927] %[output:226cc14b] %[output:9506fee2] %[output:26d53b12] %[output:7eb109e3] %[output:28546498] %[output:17230cab] %[output:69845175] %[output:9692d9d4] %[output:7d7fbc6e] %[output:60f38922] %[output:964706a5] %[output:7ae2cb99] %[output:10f0ef8e] %[output:95d71a40] %[output:5d054139] %[output:53783bb9] %[output:7ae874f8] %[output:1c8db6f5] %[output:67f2a07e] %[output:849897ba] %[output:760fa6ec] %[output:735fda90] %[output:1533646d] %[output:4ff11220] %[output:4e3ca744] %[output:45ea8f76] %[output:215d95c7] %[output:6e55758d] %[output:7c2feee2] %[output:282f93ed] %[output:05176a03] %[output:98fbcd28] %[output:85398972] %[output:9f7a7ea9] %[output:4bab330b] %[output:34ec4be8] %[output:86ee2c1b] %[output:525a0095] %[output:78215ae1] %[output:636e7ca4] %[output:585844ca] %[output:3a9417ac] %[output:22e43bf8] %[output:9ccc26e5] %[output:39cbbcfd] %[output:1426e71a] %[output:25d9bab5] %[output:2e4cd0fe] %[output:7437fc6c] %[output:909cbf79]
        dirName = [folders(ds).folder '/' folders(ds).name];
        fileList = dir(dirName);
        fileList = struct2cell(fileList);
        fileList = string(fileList(1,:));
        isFileName = contains(fileList,"txt");
        
        % find metadata file in directory
        metaDataFileName = string(dirName) + "/" + fileList(isFileName) %[output:941da4ba] %[output:4c1cfafc] %[output:5175b38c] %[output:316b3252] %[output:9cfc5872] %[output:07da783d] %[output:6b086aa6] %[output:10ad5301] %[output:91d942dc] %[output:0fccdc50] %[output:7d56cc60] %[output:0f259f13] %[output:500f8464] %[output:0d63183c] %[output:7a033b9f] %[output:03bda3e3] %[output:5285c257] %[output:9d3c5398] %[output:210ce62a] %[output:629f15f8] %[output:2b66e96c] %[output:601acca7] %[output:4cffc2be] %[output:87973830] %[output:420f047e] %[output:391a901e] %[output:6b1880c5] %[output:72240231] %[output:6c379866] %[output:3f6e6dca] %[output:49fe02d4] %[output:89a3b41e] %[output:1aaffd15] %[output:6451e192] %[output:613830e9] %[output:7b9cd4a4] %[output:91fd0c15] %[output:32170351] %[output:2d242547] %[output:49ff0b7c] %[output:7812a37e] %[output:23b0c380] %[output:7df0ff6e] %[output:1cf188a2] %[output:0a8f5764] %[output:77f8086a] %[output:1609f1ea] %[output:7c789464] %[output:1c409805] %[output:0aaa1fdc] %[output:2ce55eb8] %[output:70cd1102] %[output:542010fe] %[output:994f3653] %[output:4bad8d2b] %[output:1025e476] %[output:31154ac8] %[output:52833ed5] %[output:5e96c5f8] %[output:0b6cc115] %[output:66699fd3] %[output:74a44538] %[output:63709a6d] %[output:77e7ecf4] %[output:1ef60dec] %[output:7b756827] %[output:6840eb79] %[output:6cfc9ec6] %[output:439ce37e] %[output:3c2e5962] %[output:629279f2] %[output:2dd6dafd] %[output:5016a5c4] %[output:9e6e0ba8] %[output:1bcff863] %[output:91150468] %[output:4b55018c] %[output:48992f9e] %[output:17c6860b] %[output:9f2e4da4] %[output:1fd880d7] %[output:05b3be01] %[output:4b943ff8] %[output:1788f1ee] %[output:9a9a79a6] %[output:3ab68b1f] %[output:098a4576] %[output:15cfd33e] %[output:9a55bcea] %[output:94ad4389] %[output:2a63886f] %[output:16fbc7fd] %[output:4b7434ee] %[output:77fa5c9b] %[output:43ff6615] %[output:5836a739] %[output:0e0f5d16] %[output:525ddcaa] %[output:2e455ded] %[output:6035cdff] %[output:0dce7c04] %[output:7fbedd2f] %[output:32ce9706] %[output:9f3df855] %[output:2c365825] %[output:46d3cec8] %[output:905d6a8e] %[output:2dd90390] %[output:1963b6a5] %[output:606003fd] %[output:32e557fc] %[output:016f9ca7] %[output:5ff82d1c] %[output:3da7b489] %[output:4247064e] %[output:1a391164] %[output:43f38553] %[output:7f6eb220] %[output:7ec8fbb0] %[output:64ca3f73] %[output:194db0df] %[output:77a046bc] %[output:21dcd96a] %[output:24e2a784] %[output:3c47db5e] %[output:7ee9a9cf] %[output:223b53c9] %[output:2ecf0bc7] %[output:8cc1ec6b] %[output:41503874] %[output:703a60dc] %[output:342e5579] %[output:87116081] %[output:1e635fcf] %[output:3a9b7694] %[output:2979a345] %[output:23cdf59d] %[output:37109786] %[output:9bdc3d12] %[output:77aada5f] %[output:683d5eba] %[output:2d54ee24] %[output:0c58963f] %[output:62dfa76d] %[output:8bf79fd0] %[output:33876e28] %[output:87a882f0] %[output:386610e5] %[output:8a60523f] %[output:22676016] %[output:19c74582] %[output:2bc1bba4] %[output:17904a45] %[output:709db2ec] %[output:360eba51] %[output:6b1e556d] %[output:9827d00d] %[output:4a308f8d] %[output:929478a1] %[output:8f3a0541] %[output:3b04b451] %[output:8f6c6345] %[output:7f3f078d] %[output:2faf6203] %[output:6cfdc143] %[output:91b4bb51] %[output:6aebaf1c] %[output:4669d46e] %[output:022c8df0] %[output:35cf3fe9] %[output:21d4b765] %[output:1863283e] %[output:771495c3] %[output:6be9981c] %[output:0332419f] %[output:3e0df357] %[output:11463770] %[output:74e0d079] %[output:78c920ab] %[output:42f736be] %[output:598a6cac] %[output:2f56854d] %[output:1ff45d55] %[output:6c9aea91] %[output:45ac1728] %[output:6169765d] %[output:3f8876b2] %[output:9f079394] %[output:45a82f46] %[output:936e584a] %[output:691b70d5] %[output:0a762cf6] %[output:5594db7b] %[output:139729d5] %[output:334426c1] %[output:7392833f] %[output:82074216] %[output:4c2d4728] %[output:5445ac26] %[output:7e32feb7] %[output:965403e0] %[output:6334f0af]
        
        % loading of file
        metaDataList{ds,1} = metaDataReadTextFile(metaDataFileName);
    catch
        
    end
end %[output:group:8f271d87]

% vector that tells us if metadata is present
hasMetaData = not(cellfun(@(x) isempty(x), metaDataList));
  %[control:button:6c81]{"position":[1,2]}

%%
%[text] ## finding a certain property in a cell array of metadata cell arrays
%[text] the variable propertyValue contains all of the certain property parameters for each measurement.
%property = '/atomProbeTomography/experiment/measurementChamberPressure';
%property = '/atomProbeTomography/experiment/specimenTemperatureActual';
%property = '/atomProbeTomography/experiment/pulseType';
%property = '/atomProbeTomography/specimen/preparation/method';
property = '/atomProbeTomography/experiment/laserPulseEnergy';
%property = '/project/contactPerson';

for m=1:length(metaDataList)
    if hasMetaData(m)
        propertyRowNumber = cellfun(@(x) strcmp(x,property), metaDataList{m}(:,1));
        
        if not(isempty(metaDataList{m}{propertyRowNumber,2}))
            
            propertyValue(m,1) = metaDataList{m}{propertyRowNumber,2};
            %propertyUnit(m,1) = metaDataList{m}{propertyRowNumber,3};
            
        end
        
    end
end
  %[control:button:758b]{"position":[1,2]}
%%
%[text] ## Visualize the property cell array
%[text] After the user has searched for a certain property, the variable propertyValue can be displayed as a histogram and the user can choose 
histogram(propertyValue)
  %[control:button:5c69]{"position":[1,2]}
%%
%[text] ## Finding the (e)pos file path for each dataset folder
for ds = 1:numDatasets
    if hasMetaData(ds)
        dirName = [folders(ds).folder '/' folders(ds).name];
        fileList = dir(dirName);
        fileList = struct2cell(fileList);
        fileList = string(fileList(1,:));
        isFileName = contains(fileList,[".pos",".POS",".epos",".EPOS"]);
        
        if not(isFileName)
            disp([dirName ' does not contain a dataset'] );
        else
            posFileNames(ds,1) = string(dirName) + "/" + fileList(isFileName);
        end
        
    end
end
  %[control:button:003f]{"position":[1,2]}
%%
%[text] ## Loading ranges and ions from a figure file
[file, path] = uigetfile('*.fig','select figure file with ranges');
% loading of ranges from figure
rangeFigureVisibility = 'invisible'; % or 'visible'
rngFig = openfig([path '/' file],rangeFigureVisibility);
spec = findobj( get(rngFig,'Children'), '-depth', 2, 'DisplayName', 'mass spectrum');
rangeTable = rangesExtractFromMassSpec(spec);
ionTable = ionsExtractFromMassSpec(spec);
close(rngFig);
clear file path rangeFigureVisibility;
  %[control:button:12dd]{"position":[1,2]}

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":27.3}
%---
%[control:button:8406]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:6c81]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:758b]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:5c69]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:003f]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:12dd]
%   data: {"label":"Run","run":"Section"}
%---
%[output:68d230e2]
%   data: {"dataType":"text","outputData":{"text":"Analysis results\n","truncated":false}}
%---
%[output:941da4ba]
%   data: {"dataType":"text","outputData":{"text":"\nmetaDataFileName = \n\n  1×0 empty <a href=\"matlab:helpPopup string\" style=\"font-weight:bold\">string<\/a> array\n\n","truncated":false}}
%---
%[output:5306653f]
%   data: {"dataType":"text","outputData":{"text":"R56_00912\n","truncated":false}}
%---
%[output:4c1cfafc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_00912\/R56_00912 experiment metadata list.txt\""}}
%---
%[output:9481c78c]
%   data: {"dataType":"text","outputData":{"text":"R56_00997\n","truncated":false}}
%---
%[output:5175b38c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_00997\/R56_00997 experiment metadata list.txt\""}}
%---
%[output:55af6516]
%   data: {"dataType":"text","outputData":{"text":"R56_01173\n","truncated":false}}
%---
%[output:316b3252]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01173\/R56_01173 experiment metadata list.txt\""}}
%---
%[output:48c7ba05]
%   data: {"dataType":"text","outputData":{"text":"R56_01174\n","truncated":false}}
%---
%[output:9cfc5872]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01174\/R56_01174 experiment metadata list.txt\""}}
%---
%[output:28f0fe4f]
%   data: {"dataType":"text","outputData":{"text":"R56_01274\n","truncated":false}}
%---
%[output:07da783d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01274\/R56_01274 metadata list.txt\""}}
%---
%[output:4efdb7e9]
%   data: {"dataType":"text","outputData":{"text":"R56_01296\n","truncated":false}}
%---
%[output:6b086aa6]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01296\/R56_01296 metadata list.txt\""}}
%---
%[output:940b65f0]
%   data: {"dataType":"text","outputData":{"text":"R56_01329\n","truncated":false}}
%---
%[output:10ad5301]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01329\/R56_01329 metadata list.txt\""}}
%---
%[output:0caa6970]
%   data: {"dataType":"text","outputData":{"text":"R56_01346\n","truncated":false}}
%---
%[output:91d942dc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01346\/R56_01346 metadata list.txt\""}}
%---
%[output:76ea124a]
%   data: {"dataType":"text","outputData":{"text":"R56_01347\n","truncated":false}}
%---
%[output:0fccdc50]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01347\/R56_01347 metadata list.txt\""}}
%---
%[output:3d0e6b0a]
%   data: {"dataType":"text","outputData":{"text":"R56_01357\n","truncated":false}}
%---
%[output:7d56cc60]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01357\/R56_01357 experiment metadata list.txt\""}}
%---
%[output:13680c1a]
%   data: {"dataType":"text","outputData":{"text":"R56_01360\n","truncated":false}}
%---
%[output:0f259f13]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01360\/R56_01360 metadata list.txt\""}}
%---
%[output:954ccde1]
%   data: {"dataType":"text","outputData":{"text":"R56_01371\n","truncated":false}}
%---
%[output:500f8464]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01371\/R56_01371 metadata list.txt\""}}
%---
%[output:1ffde8bd]
%   data: {"dataType":"text","outputData":{"text":"R56_01378\n","truncated":false}}
%---
%[output:0d63183c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01378\/R56_01378 metadata list.txt\""}}
%---
%[output:59b956c1]
%   data: {"dataType":"text","outputData":{"text":"R56_01384\n","truncated":false}}
%---
%[output:7a033b9f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01384\/R56_01384 metadata list.txt\""}}
%---
%[output:7692b052]
%   data: {"dataType":"text","outputData":{"text":"R56_01395\n","truncated":false}}
%---
%[output:03bda3e3]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01395\/R56_01395 metadata list.txt\""}}
%---
%[output:5599cfe7]
%   data: {"dataType":"text","outputData":{"text":"R56_01399\n","truncated":false}}
%---
%[output:5285c257]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01399\/R56_01399 metadata list.txt\""}}
%---
%[output:8dfcac42]
%   data: {"dataType":"text","outputData":{"text":"R56_01410\n","truncated":false}}
%---
%[output:9d3c5398]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01410\/R56_01410 metadata list.txt\""}}
%---
%[output:88fb106a]
%   data: {"dataType":"text","outputData":{"text":"R56_01417\n","truncated":false}}
%---
%[output:210ce62a]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01417\/R56_01417 metadata list.txt\""}}
%---
%[output:6990bac3]
%   data: {"dataType":"text","outputData":{"text":"R56_01420\n","truncated":false}}
%---
%[output:629f15f8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01420\/R56_01420 metadata list.txt\""}}
%---
%[output:321fd6c4]
%   data: {"dataType":"text","outputData":{"text":"R56_01429\n","truncated":false}}
%---
%[output:2b66e96c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01429\/R56_01429 metadata list.txt\""}}
%---
%[output:8f63d0ac]
%   data: {"dataType":"text","outputData":{"text":"R56_01430\n","truncated":false}}
%---
%[output:601acca7]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01430\/R56_01430 metadata list.txt\""}}
%---
%[output:00ac8e4f]
%   data: {"dataType":"text","outputData":{"text":"R56_01431\n","truncated":false}}
%---
%[output:4cffc2be]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01431\/R56_01431 metadata list.txt\""}}
%---
%[output:92245aa0]
%   data: {"dataType":"text","outputData":{"text":"R56_01446\n","truncated":false}}
%---
%[output:87973830]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01446\/R56_01446 experiment metadata list.txt\""}}
%---
%[output:569e681a]
%   data: {"dataType":"text","outputData":{"text":"R56_01448\n","truncated":false}}
%---
%[output:420f047e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01448\/R56_01448 metadata list.txt\""}}
%---
%[output:0653b743]
%   data: {"dataType":"text","outputData":{"text":"R56_01452\n","truncated":false}}
%---
%[output:391a901e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01452\/R56_01452 metadata list.txt\""}}
%---
%[output:350ebc2e]
%   data: {"dataType":"text","outputData":{"text":"R56_01453\n","truncated":false}}
%---
%[output:6b1880c5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01453\/R56_01453 experiment metadata list.txt\""}}
%---
%[output:1b8f4ec8]
%   data: {"dataType":"text","outputData":{"text":"R56_01456\n","truncated":false}}
%---
%[output:72240231]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01456\/R56_01456 experiment metadata list.txt\""}}
%---
%[output:06002130]
%   data: {"dataType":"text","outputData":{"text":"R56_01515\n","truncated":false}}
%---
%[output:6c379866]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01515\/R56_01515 metadata list.txt\""}}
%---
%[output:8da7af50]
%   data: {"dataType":"text","outputData":{"text":"R56_01518\n","truncated":false}}
%---
%[output:3f6e6dca]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01518\/R56_01518 metadata list.txt\""}}
%---
%[output:1be49a22]
%   data: {"dataType":"text","outputData":{"text":"R56_01543\n","truncated":false}}
%---
%[output:49fe02d4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01543\/R56_01543 metadata list.txt\""}}
%---
%[output:75d5136e]
%   data: {"dataType":"text","outputData":{"text":"R56_01544\n","truncated":false}}
%---
%[output:89a3b41e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01544\/R56_01544 metadata list.txt\""}}
%---
%[output:646094a1]
%   data: {"dataType":"text","outputData":{"text":"R56_01554\n","truncated":false}}
%---
%[output:1aaffd15]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01554\/R56_01554 metadata list.txt\""}}
%---
%[output:2564afd7]
%   data: {"dataType":"text","outputData":{"text":"R56_01557\n","truncated":false}}
%---
%[output:6451e192]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01557\/R56_01557 experiment metadata list.txt\""}}
%---
%[output:91d0c8ee]
%   data: {"dataType":"text","outputData":{"text":"R56_01577\n","truncated":false}}
%---
%[output:613830e9]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01577\/R56_01577 metadata list.txt\""}}
%---
%[output:89dfb9c8]
%   data: {"dataType":"text","outputData":{"text":"R56_01578\n","truncated":false}}
%---
%[output:7b9cd4a4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01578\/R56_01578 metadata list.txt\""}}
%---
%[output:3ce91317]
%   data: {"dataType":"text","outputData":{"text":"R56_01597\n","truncated":false}}
%---
%[output:91fd0c15]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01597\/R56_01597 metadata list.txt\""}}
%---
%[output:7a713631]
%   data: {"dataType":"text","outputData":{"text":"R56_01599\n","truncated":false}}
%---
%[output:32170351]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01599\/R56_01599 experiment metadata list.txt\""}}
%---
%[output:9d4f8109]
%   data: {"dataType":"text","outputData":{"text":"R56_01600\n","truncated":false}}
%---
%[output:2d242547]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01600\/R56_01600 experiment metadata list.txt\""}}
%---
%[output:9b4080db]
%   data: {"dataType":"text","outputData":{"text":"R56_01604\n","truncated":false}}
%---
%[output:49ff0b7c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01604\/R56_01604 experiment metadata list.txt\""}}
%---
%[output:16ce2af9]
%   data: {"dataType":"text","outputData":{"text":"R56_01616\n","truncated":false}}
%---
%[output:7812a37e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01616\/R56_01616 metadata list.txt\""}}
%---
%[output:62552850]
%   data: {"dataType":"text","outputData":{"text":"R56_01618\n","truncated":false}}
%---
%[output:23b0c380]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01618\/R56_01618 metadata list.txt\""}}
%---
%[output:2b60dc7b]
%   data: {"dataType":"text","outputData":{"text":"R56_01660\n","truncated":false}}
%---
%[output:7df0ff6e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01660\/R56_01660 experiment metadata list.txt\""}}
%---
%[output:63e5a11b]
%   data: {"dataType":"text","outputData":{"text":"R56_01665\n","truncated":false}}
%---
%[output:1cf188a2]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01665\/R56_01665 experiment metadata list.txt\""}}
%---
%[output:96aaa77a]
%   data: {"dataType":"text","outputData":{"text":"R56_01685\n","truncated":false}}
%---
%[output:0a8f5764]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01685\/R56_01685 metadata list.txt\""}}
%---
%[output:26d45876]
%   data: {"dataType":"text","outputData":{"text":"R56_01688\n","truncated":false}}
%---
%[output:77f8086a]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01688\/R56_01688 metadata list.txt\""}}
%---
%[output:3a142d02]
%   data: {"dataType":"text","outputData":{"text":"R56_01700\n","truncated":false}}
%---
%[output:1609f1ea]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01700\/R56_01700 experiment metadata list.txt\""}}
%---
%[output:654d5778]
%   data: {"dataType":"text","outputData":{"text":"R56_01715\n","truncated":false}}
%---
%[output:7c789464]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01715\/R56_01715 experiment metadata list.txt\""}}
%---
%[output:249f44b0]
%   data: {"dataType":"text","outputData":{"text":"R56_01785\n","truncated":false}}
%---
%[output:1c409805]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01785\/R56_01785 experiment metadata list.txt\""}}
%---
%[output:5b945e1b]
%   data: {"dataType":"text","outputData":{"text":"R56_01787\n","truncated":false}}
%---
%[output:0aaa1fdc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01787\/R56_01787 experiment metadata list.txt\""}}
%---
%[output:5f8f1b03]
%   data: {"dataType":"text","outputData":{"text":"R56_01820\n","truncated":false}}
%---
%[output:2ce55eb8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01820\/R56_01820 experiment metadata list.txt\""}}
%---
%[output:24d20803]
%   data: {"dataType":"text","outputData":{"text":"R56_01827\n","truncated":false}}
%---
%[output:70cd1102]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01827\/R56_01827 experiment metadata list.txt\""}}
%---
%[output:370a26c2]
%   data: {"dataType":"text","outputData":{"text":"R56_01863\n","truncated":false}}
%---
%[output:542010fe]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01863\/R56_01863 experiment metadata list.txt\""}}
%---
%[output:9b7267bb]
%   data: {"dataType":"text","outputData":{"text":"R56_01865\n","truncated":false}}
%---
%[output:994f3653]
%   data: {"dataType":"matrix","outputData":{"columns":2,"header":"1×2 string array","name":"metaDataFileName","rows":1,"type":"string","value":[["L:\\\/R56_01865\/R56_01908 experiment metadata list.txt","L:\\\/R56_01865\/R56_01965 experiment metadata list.txt"]]}}
%---
%[output:7f9c8c28]
%   data: {"dataType":"text","outputData":{"text":"R56_01908\n","truncated":false}}
%---
%[output:4bad8d2b]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01908\/R56_01908 experiment metadata list.txt\""}}
%---
%[output:78520967]
%   data: {"dataType":"text","outputData":{"text":"R56_01943\n","truncated":false}}
%---
%[output:1025e476]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01943\/R56_01943 experiment metadata list.txt\""}}
%---
%[output:313ea599]
%   data: {"dataType":"text","outputData":{"text":"R56_01964\n","truncated":false}}
%---
%[output:31154ac8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01964\/R56_01964 metadata list.txt\""}}
%---
%[output:8ae247e2]
%   data: {"dataType":"text","outputData":{"text":"R56_01974\n","truncated":false}}
%---
%[output:52833ed5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_01974\/R56_01974 metadata list.txt\""}}
%---
%[output:4d2abc5c]
%   data: {"dataType":"text","outputData":{"text":"R56_02020\n","truncated":false}}
%---
%[output:5e96c5f8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02020\/R56_02020 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:5bf4dd0c]
%   data: {"dataType":"text","outputData":{"text":"R56_02039\n","truncated":false}}
%---
%[output:0b6cc115]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02039\/R56_02039 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:038e74cb]
%   data: {"dataType":"text","outputData":{"text":"R56_02044\n","truncated":false}}
%---
%[output:66699fd3]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02044\/R56_02044 metadata list.txt\""}}
%---
%[output:9b515301]
%   data: {"dataType":"text","outputData":{"text":"R56_02078\n","truncated":false}}
%---
%[output:74a44538]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02078\/R56_02078 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:5fdb23d5]
%   data: {"dataType":"text","outputData":{"text":"R56_02081\n","truncated":false}}
%---
%[output:63709a6d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02081\/R56_02081 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:3ef9abad]
%   data: {"dataType":"text","outputData":{"text":"R56_02095\n","truncated":false}}
%---
%[output:77e7ecf4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02095\/R56_02095 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:5638629e]
%   data: {"dataType":"text","outputData":{"text":"R56_02099\n","truncated":false}}
%---
%[output:1ef60dec]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02099\/R56_02099 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:8bed0e7f]
%   data: {"dataType":"text","outputData":{"text":"R56_02203\n","truncated":false}}
%---
%[output:7b756827]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02203\/R56_02203 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:8e596f47]
%   data: {"dataType":"text","outputData":{"text":"R56_02210\n","truncated":false}}
%---
%[output:6840eb79]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02210\/R56_02210 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:5e6b1012]
%   data: {"dataType":"text","outputData":{"text":"R56_02213\n","truncated":false}}
%---
%[output:6cfc9ec6]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02213\/R56_02213 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:8ac17693]
%   data: {"dataType":"text","outputData":{"text":"R56_02217\n","truncated":false}}
%---
%[output:439ce37e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02217\/R56_02217 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:779c9df5]
%   data: {"dataType":"text","outputData":{"text":"R56_02220\n","truncated":false}}
%---
%[output:3c2e5962]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02220\/R56_02220 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1de604e1]
%   data: {"dataType":"text","outputData":{"text":"R56_02289\n","truncated":false}}
%---
%[output:629279f2]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02289\/R56_02289 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:27f719ba]
%   data: {"dataType":"text","outputData":{"text":"R56_02297\n","truncated":false}}
%---
%[output:2dd6dafd]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02297\/R56_02297 metadata list.txt\""}}
%---
%[output:3a560479]
%   data: {"dataType":"text","outputData":{"text":"R56_02303\n","truncated":false}}
%---
%[output:5016a5c4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02303\/R56_02303 metadata list.txt\""}}
%---
%[output:5255a2ef]
%   data: {"dataType":"text","outputData":{"text":"R56_02304\n","truncated":false}}
%---
%[output:9e6e0ba8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02304\/R56_02304 metadata list.txt\""}}
%---
%[output:4fda71c6]
%   data: {"dataType":"text","outputData":{"text":"R56_02306\n","truncated":false}}
%---
%[output:1bcff863]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02306\/R56_02306 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4e834c5b]
%   data: {"dataType":"text","outputData":{"text":"R56_02317\n","truncated":false}}
%---
%[output:91150468]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02317\/R56_02317 metadata list.txt\""}}
%---
%[output:6d5af7ae]
%   data: {"dataType":"text","outputData":{"text":"R56_02319\n","truncated":false}}
%---
%[output:4b55018c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02319\/R56_02319 metadata list.txt\""}}
%---
%[output:32b33c27]
%   data: {"dataType":"text","outputData":{"text":"R56_02322\n","truncated":false}}
%---
%[output:48992f9e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02322\/R56_02322 metadata list.txt\""}}
%---
%[output:7b6252c0]
%   data: {"dataType":"text","outputData":{"text":"R56_02325\n","truncated":false}}
%---
%[output:17c6860b]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02325\/R56_02325 metadata list.txt\""}}
%---
%[output:0bc197ae]
%   data: {"dataType":"text","outputData":{"text":"R56_02329\n","truncated":false}}
%---
%[output:9f2e4da4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02329\/R56_02329 metadata list.txt\""}}
%---
%[output:3eae81dc]
%   data: {"dataType":"text","outputData":{"text":"R56_02331\n","truncated":false}}
%---
%[output:1fd880d7]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02331\/R56_02331 metadata list.txt\""}}
%---
%[output:7b6a6c82]
%   data: {"dataType":"text","outputData":{"text":"R56_02346\n","truncated":false}}
%---
%[output:05b3be01]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02346\/R56_02346 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1c89fb78]
%   data: {"dataType":"text","outputData":{"text":"R56_02366\n","truncated":false}}
%---
%[output:4b943ff8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02366\/R56_02366 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:03ecb636]
%   data: {"dataType":"text","outputData":{"text":"R56_02374\n","truncated":false}}
%---
%[output:1788f1ee]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02374\/R56_02374 metadata list.txt\""}}
%---
%[output:1623a0d0]
%   data: {"dataType":"text","outputData":{"text":"R56_02379\n","truncated":false}}
%---
%[output:9a9a79a6]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02379\/R56_02436 metadata list.txt\""}}
%---
%[output:1f49a2e0]
%   data: {"dataType":"text","outputData":{"text":"R56_02436\n","truncated":false}}
%---
%[output:3ab68b1f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02436\/R56_02436 metadata list.txt\""}}
%---
%[output:29213785]
%   data: {"dataType":"text","outputData":{"text":"R56_02461\n","truncated":false}}
%---
%[output:098a4576]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02461\/R56_02461 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:49e5665e]
%   data: {"dataType":"text","outputData":{"text":"R56_02464\n","truncated":false}}
%---
%[output:15cfd33e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02464\/R56_02464 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:161f109d]
%   data: {"dataType":"text","outputData":{"text":"R56_02472\n","truncated":false}}
%---
%[output:9a55bcea]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02472\/R56_02472 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:63837842]
%   data: {"dataType":"text","outputData":{"text":"R56_02476\n","truncated":false}}
%---
%[output:94ad4389]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02476\/R56_02476 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:82e0445a]
%   data: {"dataType":"text","outputData":{"text":"R56_02489\n","truncated":false}}
%---
%[output:2a63886f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02489\/R56_02489 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:9d171763]
%   data: {"dataType":"text","outputData":{"text":"R56_02492\n","truncated":false}}
%---
%[output:16fbc7fd]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02492\/R56_02492 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:551918c7]
%   data: {"dataType":"text","outputData":{"text":"R56_02507\n","truncated":false}}
%---
%[output:4b7434ee]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02507\/R56_02507 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:0ca43578]
%   data: {"dataType":"text","outputData":{"text":"R56_02508\n","truncated":false}}
%---
%[output:77fa5c9b]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02508\/R56_02508 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:579fd785]
%   data: {"dataType":"text","outputData":{"text":"R56_02509\n","truncated":false}}
%---
%[output:43ff6615]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02509\/R56_02509 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4c5bcb03]
%   data: {"dataType":"text","outputData":{"text":"R56_02533\n","truncated":false}}
%---
%[output:5836a739]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02533\/R56_02533 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4886bb3f]
%   data: {"dataType":"text","outputData":{"text":"R56_02536\n","truncated":false}}
%---
%[output:0e0f5d16]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02536\/R56_02536 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1356f23f]
%   data: {"dataType":"text","outputData":{"text":"R56_02544\n","truncated":false}}
%---
%[output:525ddcaa]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02544\/R56_02544 metadata list.txt\""}}
%---
%[output:01ff7d9c]
%   data: {"dataType":"text","outputData":{"text":"R56_02551\n","truncated":false}}
%---
%[output:2e455ded]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02551\/R56_02551 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:74847e87]
%   data: {"dataType":"text","outputData":{"text":"R56_02580\n","truncated":false}}
%---
%[output:6035cdff]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02580\/R56_02580 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7ce5af0d]
%   data: {"dataType":"text","outputData":{"text":"R56_02584\n","truncated":false}}
%---
%[output:0dce7c04]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02584\/R56_02584 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:42933222]
%   data: {"dataType":"text","outputData":{"text":"R56_02591\n","truncated":false}}
%---
%[output:7fbedd2f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02591\/R56_02591 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:29f0aa33]
%   data: {"dataType":"text","outputData":{"text":"R56_02592\n","truncated":false}}
%---
%[output:32ce9706]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02592\/R56_02592 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:103e1fc9]
%   data: {"dataType":"text","outputData":{"text":"R56_02596\n","truncated":false}}
%---
%[output:9f3df855]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02596\/R56_02596 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:279abdfc]
%   data: {"dataType":"text","outputData":{"text":"R56_02597\n","truncated":false}}
%---
%[output:2c365825]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02597\/R56_02597 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:9f9d5e2d]
%   data: {"dataType":"text","outputData":{"text":"R56_02614\n","truncated":false}}
%---
%[output:46d3cec8]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02614\/R56_02614 metadata list.txt\""}}
%---
%[output:416594ed]
%   data: {"dataType":"text","outputData":{"text":"R56_02615\n","truncated":false}}
%---
%[output:905d6a8e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02615\/R56_02615 metadata list.txt\""}}
%---
%[output:43d464e5]
%   data: {"dataType":"text","outputData":{"text":"R56_02641\n","truncated":false}}
%---
%[output:2dd90390]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02641\/R56_02641 metadata list.txt\""}}
%---
%[output:06c22b32]
%   data: {"dataType":"text","outputData":{"text":"R56_02647\n","truncated":false}}
%---
%[output:1963b6a5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02647\/R56_02647 metadata list.txt\""}}
%---
%[output:35283f5e]
%   data: {"dataType":"text","outputData":{"text":"R56_02676\n","truncated":false}}
%---
%[output:606003fd]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02676\/R56_02676 metadata list.txt\""}}
%---
%[output:03d9da65]
%   data: {"dataType":"text","outputData":{"text":"R56_02713\n","truncated":false}}
%---
%[output:32e557fc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02713\/R56_02713 metadata list.txt\""}}
%---
%[output:7fd99ec1]
%   data: {"dataType":"text","outputData":{"text":"R56_02714\n","truncated":false}}
%---
%[output:016f9ca7]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02714\/R56_02714 metadata list.txt\""}}
%---
%[output:2308259a]
%   data: {"dataType":"text","outputData":{"text":"R56_02715\n","truncated":false}}
%---
%[output:5ff82d1c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02715\/R56_02715 metadata list.txt\""}}
%---
%[output:0df3c767]
%   data: {"dataType":"text","outputData":{"text":"R56_02716\n","truncated":false}}
%---
%[output:3da7b489]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02716\/R56_02716 metadata list.txt\""}}
%---
%[output:6fe13e16]
%   data: {"dataType":"text","outputData":{"text":"R56_02720\n","truncated":false}}
%---
%[output:4247064e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02720\/R56_02720 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7c91d561]
%   data: {"dataType":"text","outputData":{"text":"R56_02727\n","truncated":false}}
%---
%[output:1a391164]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02727\/R56_02727 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:265bbd92]
%   data: {"dataType":"text","outputData":{"text":"R56_02730\n","truncated":false}}
%---
%[output:43f38553]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02730\/R56_02730 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:72bc1e81]
%   data: {"dataType":"text","outputData":{"text":"R56_02734\n","truncated":false}}
%---
%[output:7f6eb220]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02734\/R56_02734 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7639c7a8]
%   data: {"dataType":"text","outputData":{"text":"R56_02735\n","truncated":false}}
%---
%[output:7ec8fbb0]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02735\/R56_02735 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:2132a507]
%   data: {"dataType":"text","outputData":{"text":"R56_02738\n","truncated":false}}
%---
%[output:64ca3f73]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02738\/R56_02738 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:807e3719]
%   data: {"dataType":"text","outputData":{"text":"R56_02750\n","truncated":false}}
%---
%[output:194db0df]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02750\/R56_02750 metadata list.txt\""}}
%---
%[output:7d983c8f]
%   data: {"dataType":"text","outputData":{"text":"R56_02758\n","truncated":false}}
%---
%[output:77a046bc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02758\/R56_02758 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:03a32c0c]
%   data: {"dataType":"text","outputData":{"text":"R56_02763\n","truncated":false}}
%---
%[output:21dcd96a]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02763\/R56_02763 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:02151087]
%   data: {"dataType":"text","outputData":{"text":"R56_02786\n","truncated":false}}
%---
%[output:24e2a784]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02786\/R56_02786 metadata list.txt\""}}
%---
%[output:5230dd18]
%   data: {"dataType":"text","outputData":{"text":"R56_02803\n","truncated":false}}
%---
%[output:3c47db5e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02803\/R56_02803 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:9fdc90e2]
%   data: {"dataType":"text","outputData":{"text":"R56_02841\n","truncated":false}}
%---
%[output:7ee9a9cf]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02841\/R56_02841 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:74d6da98]
%   data: {"dataType":"text","outputData":{"text":"R56_02843\n","truncated":false}}
%---
%[output:223b53c9]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02843\/R56_02843 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:47cb410f]
%   data: {"dataType":"text","outputData":{"text":"R56_02845\n","truncated":false}}
%---
%[output:2ecf0bc7]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02845\/R56_02845 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:11db43f0]
%   data: {"dataType":"text","outputData":{"text":"R56_02875\n","truncated":false}}
%---
%[output:8cc1ec6b]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02875\/R56_02875 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7104c10f]
%   data: {"dataType":"text","outputData":{"text":"R56_02876\n","truncated":false}}
%---
%[output:41503874]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02876\/R56_02876 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:07301523]
%   data: {"dataType":"text","outputData":{"text":"R56_02877\n","truncated":false}}
%---
%[output:703a60dc]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02877\/R56_02877 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:565158d4]
%   data: {"dataType":"text","outputData":{"text":"R56_02898\n","truncated":false}}
%---
%[output:342e5579]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_02898\/R56_02898 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:8127d8c4]
%   data: {"dataType":"text","outputData":{"text":"R56_03001\n","truncated":false}}
%---
%[output:87116081]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03001\/R56_03001 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:2bfd026b]
%   data: {"dataType":"text","outputData":{"text":"R56_03019\n","truncated":false}}
%---
%[output:1e635fcf]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03019\/R56_03019 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:3e1ffa50]
%   data: {"dataType":"text","outputData":{"text":"R56_03026\n","truncated":false}}
%---
%[output:3a9b7694]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03026\/R56_03026 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:39ab5672]
%   data: {"dataType":"text","outputData":{"text":"R56_03037\n","truncated":false}}
%---
%[output:2979a345]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03037\/R56_03037 metadata list.txt\""}}
%---
%[output:27592aeb]
%   data: {"dataType":"text","outputData":{"text":"R56_03039\n","truncated":false}}
%---
%[output:23cdf59d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03039\/R56_03039 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:58b8a750]
%   data: {"dataType":"text","outputData":{"text":"R56_03044\n","truncated":false}}
%---
%[output:37109786]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03044\/R56_03044 metadata list.txt\""}}
%---
%[output:36053a74]
%   data: {"dataType":"text","outputData":{"text":"R56_03049\n","truncated":false}}
%---
%[output:9bdc3d12]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03049\/R56_03049 metadata list.txt\""}}
%---
%[output:29e06349]
%   data: {"dataType":"text","outputData":{"text":"R56_03057\n","truncated":false}}
%---
%[output:77aada5f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03057\/R56_03057 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:2cab0ae8]
%   data: {"dataType":"text","outputData":{"text":"R56_03062\n","truncated":false}}
%---
%[output:683d5eba]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03062\/R56_03062 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:26ca7546]
%   data: {"dataType":"text","outputData":{"text":"R56_03063\n","truncated":false}}
%---
%[output:2d54ee24]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03063\/R56_03063 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:6eccdf08]
%   data: {"dataType":"text","outputData":{"text":"R56_03064\n","truncated":false}}
%---
%[output:0c58963f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03064\/R56_03064 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:500a5c1e]
%   data: {"dataType":"text","outputData":{"text":"R56_03071\n","truncated":false}}
%---
%[output:62dfa76d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03071\/R56_03071 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:98df6923]
%   data: {"dataType":"text","outputData":{"text":"R56_03073\n","truncated":false}}
%---
%[output:8bf79fd0]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03073\/R56_03073 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1810de4b]
%   data: {"dataType":"text","outputData":{"text":"R56_03074\n","truncated":false}}
%---
%[output:33876e28]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03074\/R56_03074 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4229fb09]
%   data: {"dataType":"text","outputData":{"text":"R56_03078\n","truncated":false}}
%---
%[output:87a882f0]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03078\/R56_03078 metadata list.txt\""}}
%---
%[output:3928dc91]
%   data: {"dataType":"text","outputData":{"text":"R56_03103\n","truncated":false}}
%---
%[output:386610e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03103\/R56_03103 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:3dbdf65f]
%   data: {"dataType":"text","outputData":{"text":"R56_03139\n","truncated":false}}
%---
%[output:8a60523f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03139\/R56_03139 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1624d433]
%   data: {"dataType":"text","outputData":{"text":"R56_03153\n","truncated":false}}
%---
%[output:22676016]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03153\/R56_03153 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:674ccd0f]
%   data: {"dataType":"text","outputData":{"text":"R56_03154\n","truncated":false}}
%---
%[output:19c74582]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03154\/R56_03154 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:98724927]
%   data: {"dataType":"text","outputData":{"text":"R56_03156\n","truncated":false}}
%---
%[output:2bc1bba4]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03156\/R56_03156 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:226cc14b]
%   data: {"dataType":"text","outputData":{"text":"R56_03158\n","truncated":false}}
%---
%[output:17904a45]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03158\/R56_03158 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:9506fee2]
%   data: {"dataType":"text","outputData":{"text":"R56_03160\n","truncated":false}}
%---
%[output:709db2ec]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03160\/R56_03160 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:26d53b12]
%   data: {"dataType":"text","outputData":{"text":"R56_03186\n","truncated":false}}
%---
%[output:360eba51]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03186\/R56_03186 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7eb109e3]
%   data: {"dataType":"text","outputData":{"text":"R56_03188\n","truncated":false}}
%---
%[output:6b1e556d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03188\/R56_03188 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:28546498]
%   data: {"dataType":"text","outputData":{"text":"R56_03191\n","truncated":false}}
%---
%[output:9827d00d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03191\/R56_03191 metadata list.txt\""}}
%---
%[output:17230cab]
%   data: {"dataType":"text","outputData":{"text":"R56_03216\n","truncated":false}}
%---
%[output:4a308f8d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03216\/R56_03216 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:69845175]
%   data: {"dataType":"text","outputData":{"text":"R56_03221\n","truncated":false}}
%---
%[output:929478a1]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03221\/R56_03221 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:9692d9d4]
%   data: {"dataType":"text","outputData":{"text":"R56_03228\n","truncated":false}}
%---
%[output:8f3a0541]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03228\/R56_03228 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7d7fbc6e]
%   data: {"dataType":"text","outputData":{"text":"R56_03258\n","truncated":false}}
%---
%[output:3b04b451]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03258\/R56_03258 metadata list.txt\""}}
%---
%[output:60f38922]
%   data: {"dataType":"text","outputData":{"text":"R56_03260\n","truncated":false}}
%---
%[output:8f6c6345]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03260\/R56_03260 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:964706a5]
%   data: {"dataType":"text","outputData":{"text":"R56_03269\n","truncated":false}}
%---
%[output:7f3f078d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03269\/R56_03269 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7ae2cb99]
%   data: {"dataType":"text","outputData":{"text":"R56_03275\n","truncated":false}}
%---
%[output:2faf6203]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03275\/R56_03275 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:10f0ef8e]
%   data: {"dataType":"text","outputData":{"text":"R56_03276\n","truncated":false}}
%---
%[output:6cfdc143]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03276\/R56_03276 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:95d71a40]
%   data: {"dataType":"text","outputData":{"text":"R56_03279\n","truncated":false}}
%---
%[output:91b4bb51]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03279\/R56_03279 metadata list.txt\""}}
%---
%[output:5d054139]
%   data: {"dataType":"text","outputData":{"text":"R56_03290\n","truncated":false}}
%---
%[output:6aebaf1c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03290\/R56_03290 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:53783bb9]
%   data: {"dataType":"text","outputData":{"text":"R56_03315\n","truncated":false}}
%---
%[output:4669d46e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03315\/R56_03315 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:7ae874f8]
%   data: {"dataType":"text","outputData":{"text":"R56_03337\n","truncated":false}}
%---
%[output:022c8df0]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03337\/R56_03337 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1c8db6f5]
%   data: {"dataType":"text","outputData":{"text":"R56_03339\n","truncated":false}}
%---
%[output:35cf3fe9]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03339\/R56_03339 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:67f2a07e]
%   data: {"dataType":"text","outputData":{"text":"R56_03376\n","truncated":false}}
%---
%[output:21d4b765]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03376\/R56_03376 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:849897ba]
%   data: {"dataType":"text","outputData":{"text":"R56_03382\n","truncated":false}}
%---
%[output:1863283e]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03382\/R56_03382 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:760fa6ec]
%   data: {"dataType":"text","outputData":{"text":"R56_03383\n","truncated":false}}
%---
%[output:771495c3]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03383\/R56_03383 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:735fda90]
%   data: {"dataType":"text","outputData":{"text":"R56_03384\n","truncated":false}}
%---
%[output:6be9981c]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03384\/R56_03384 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:1533646d]
%   data: {"dataType":"text","outputData":{"text":"R56_03399\n","truncated":false}}
%---
%[output:0332419f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03399\/R56_03399 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4ff11220]
%   data: {"dataType":"text","outputData":{"text":"R56_03403\n","truncated":false}}
%---
%[output:3e0df357]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03403\/R56_03403 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:4e3ca744]
%   data: {"dataType":"text","outputData":{"text":"R56_03404\n","truncated":false}}
%---
%[output:11463770]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03404\/R56_03404 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:45ea8f76]
%   data: {"dataType":"text","outputData":{"text":"R56_03405\n","truncated":false}}
%---
%[output:74e0d079]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03405\/R56_03405 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:215d95c7]
%   data: {"dataType":"text","outputData":{"text":"R56_03416\n","truncated":false}}
%---
%[output:78c920ab]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03416\/R56_03416 metadata list.txt\""}}
%---
%[output:6e55758d]
%   data: {"dataType":"text","outputData":{"text":"R56_03418\n","truncated":false}}
%---
%[output:42f736be]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03418\/R56_03418 metadata list.txt\""}}
%---
%[output:7c2feee2]
%   data: {"dataType":"text","outputData":{"text":"R56_03423\n","truncated":false}}
%---
%[output:598a6cac]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03423\/R56_03423 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:282f93ed]
%   data: {"dataType":"text","outputData":{"text":"R56_03426\n","truncated":false}}
%---
%[output:2f56854d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03426\/R56_03426 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:05176a03]
%   data: {"dataType":"text","outputData":{"text":"R56_03436\n","truncated":false}}
%---
%[output:1ff45d55]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03436\/R56_03436 experiment metadata list.txt\""}}
%---
%[output:98fbcd28]
%   data: {"dataType":"text","outputData":{"text":"R56_03462\n","truncated":false}}
%---
%[output:6c9aea91]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03462\/R56_03462 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:85398972]
%   data: {"dataType":"text","outputData":{"text":"R56_03470\n","truncated":false}}
%---
%[output:45ac1728]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03470\/R56_03470 experiment metadata list.txt\""}}
%---
%[output:9f7a7ea9]
%   data: {"dataType":"text","outputData":{"text":"R56_03474\n","truncated":false}}
%---
%[output:6169765d]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03474\/R56_03474 experiment metadata list.txt\""}}
%---
%[output:4bab330b]
%   data: {"dataType":"text","outputData":{"text":"R56_03475\n","truncated":false}}
%---
%[output:3f8876b2]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03475\/R56_03475 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:34ec4be8]
%   data: {"dataType":"text","outputData":{"text":"R56_03478\n","truncated":false}}
%---
%[output:9f079394]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03478\/R56_03478 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:86ee2c1b]
%   data: {"dataType":"text","outputData":{"text":"R56_03481\n","truncated":false}}
%---
%[output:45a82f46]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03481\/R56_03481 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:525a0095]
%   data: {"dataType":"text","outputData":{"text":"R56_03496\n","truncated":false}}
%---
%[output:936e584a]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03496\/R56_03496 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:78215ae1]
%   data: {"dataType":"text","outputData":{"text":"R56_03504\n","truncated":false}}
%---
%[output:691b70d5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03504\/R56_03504 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:636e7ca4]
%   data: {"dataType":"text","outputData":{"text":"R56_03505\n","truncated":false}}
%---
%[output:0a762cf6]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03505\/R56_03505 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:585844ca]
%   data: {"dataType":"text","outputData":{"text":"R56_03529\n","truncated":false}}
%---
%[output:5594db7b]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03529\/R56_03529 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:3a9417ac]
%   data: {"dataType":"text","outputData":{"text":"R56_03530\n","truncated":false}}
%---
%[output:139729d5]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03530\/R56_03530 experiment metadata list_2020-08-24.txt\""}}
%---
%[output:22e43bf8]
%   data: {"dataType":"text","outputData":{"text":"R56_03533\n","truncated":false}}
%---
%[output:334426c1]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03533\/R56_03533 experiment metadata list.txt\""}}
%---
%[output:9ccc26e5]
%   data: {"dataType":"text","outputData":{"text":"R56_03535\n","truncated":false}}
%---
%[output:7392833f]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03535\/R56_03535 experiment metadata list.txt\""}}
%---
%[output:39cbbcfd]
%   data: {"dataType":"text","outputData":{"text":"R56_03550\n","truncated":false}}
%---
%[output:82074216]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03550\/R56_03550 metadata list.txt\""}}
%---
%[output:1426e71a]
%   data: {"dataType":"text","outputData":{"text":"R56_03568\n","truncated":false}}
%---
%[output:4c2d4728]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03568\/R56_03568 experiment metadata list.txt\""}}
%---
%[output:25d9bab5]
%   data: {"dataType":"text","outputData":{"text":"R56_03569\n","truncated":false}}
%---
%[output:5445ac26]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03569\/R56_03569 experiment metadata list.txt\""}}
%---
%[output:2e4cd0fe]
%   data: {"dataType":"text","outputData":{"text":"R56_03570\n","truncated":false}}
%---
%[output:7e32feb7]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03570\/R56_03570 experiment metadata list.txt\""}}
%---
%[output:7437fc6c]
%   data: {"dataType":"text","outputData":{"text":"R56_03571\n","truncated":false}}
%---
%[output:965403e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03571\/R56_03571 metadata list.txt\""}}
%---
%[output:909cbf79]
%   data: {"dataType":"text","outputData":{"text":"R56_03622\n","truncated":false}}
%---
%[output:6334f0af]
%   data: {"dataType":"textualVariable","outputData":{"name":"metaDataFileName","value":"\"L:\\\/R56_03622\/R56_03622 metadata list.txt\""}}
%---
