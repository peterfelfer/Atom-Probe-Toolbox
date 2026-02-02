%[text] ## Workflow Calculate the H peak concentration with Background correction and automatic range adjustment
%[text] 
%[text] This workflow takes each HDF5 file and calculates with the function peakBGCorrectedCount the background corrected concentration of the Hydrogen peak of each dataset @1 Da. If no range is given in the range tabel for the dataset - a predefined range is used from 0.95 to 1.13 Da
%[text] 1. Find List with HDF5 files
%[text] 2. Load file 
%[text] 3. check if rangeTable is there 
%[text] 4. Perform calculation 
%[text] 5. Find List with HDF5 files \
%%
%[text] 1\. Find List with HDF5 files
fileH5 = struct2table(dir('**/*.h5'));

% make rangeTable - in Case the h5 has no range Table stored 
    rangeName = categorical({'H'});
    chargeState = 1 ;
    mcbegin = 0.95 ;
    mcend = 1.13;
    volume = 0;
    ion{1,:} = ionConvertName('1H +');
    color = [0.6981 0.6665 0.1781];
    rangeTableH = table(rangeName, chargeState, mcbegin, mcend, volume, ion, color);
    clearvars rangeName chargeState mcbegin mcend volume ion color

% make dummy table to store the calculated peak data
    dum = zeros(height(fileH5),1);
    peakDataAll = table(fileH5.name, dum, dum, dum, dum, dum, dum, dum, dum, 'VariableNames', {'name', 'mcbegin', 'mcend', 'counts','countsNew', 'pct', 'pctNew', 'ppmNew' 'loc'});

%%
%[text] 2\. Load file
for i = 1:height(fileH5)
    
    % load pos
    pos = posTableFromHDF5(fileH5.name{i});

    % check if rangeTable is there
    rangesExist = h5info(fileH5.name{i}, '/atomProbeTomography/massToChargeRange');
    if isempty(rangesExist.Groups)
        rangeTable = rangeTableH;
    else
        rangeTable = rangeTableFromHDF5(fileH5.name{i});
        % check if H is ranged
        if ~sum((ismember(rangeTable.rangeName, 'H') + any(rangeTable.chargeState(rangeTable.rangeName == 'H') == 1)) == 2) ==1
            rangeTable = rangeTableH;         
        end
    end

    % 3. Claculate the H concentration 
    [peakData] = peakBGCorrectedCount(pos, rangeTable, [0.5 0.8], '1H +', false, 'a');

    % 4. Include data for each dataset

    if peakDataAll.name{i} == fileH5.name{i}
       peakDataAll.mcbegin(i) = peakData.mcbegin;
       peakDataAll.mcend(i) = peakData.mcend;
       peakDataAll.counts(i) = peakData.counts;
       peakDataAll.pct(i) = peakData.pct;
       peakDataAll.loc(i) = peakData.loc;
       peakDataAll.countsNew(i) = peakData.countsNew;
       peakDataAll.pctNew(i) = peakData.pctNew;
       peakDataAll.ppmNew(i) = peakData.ppmNew;      
       
    else
        error('file ordering is a problem')
    end

end

% save peakDataAll

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":38}
%---
