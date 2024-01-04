function datasetSliceMassSpec(posIn,varargin)
% datasetSliceMassSpec lets you slice your dataset in a defined intervall
% size and creates for each slice a figure with the respective mass Spectra
% A specific region of interestet of the mass spectra can also be choosen
% to highlight certain areas of interest
%
% datasetSliceMassSpec(posIn,ionTable, colorScheme, isotopeTable)
% datasetSliceMassSpec(posIn,ionTable, colorScheme, isotopeTable, interval)
% datasetSliceMassSpec(posIn,ionTable, colorScheme, isotopeTable, interval, axisROI, dataName, binWidth, mode)
% 
% INPUT
% posIn = pos file
% ionTable = ion Table of the pos file if you want the ions to be displayed
% colorScheme = colorScheme
% isotopeTable = isotopic table with the natural abundances
% interval = how many ions for one massSpec (min 100 000)
% axisROI = axis of the desired region of interest of the massSpec
% dataName = Name of the dataset
% binWidth = binWidth of the massSpec
% mode = 'normalised' or 'counts'
%
% OUTPUT
% figures with the massSpec of different dataset slices

%% Basics - Define Input Variables
s = 1;
interval = 2000000;
binWidth = 0.1; 
sampleName = 'sample';
mode = 'normalised';


for k = 1:length(varargin)
     % check for the tables
    if istable(varargin{k}) == 1
        tableSort = varargin{k};
        if ismember('ionName',tableSort.Properties.VariableNames) && ismember('chargeState',tableSort.Properties.VariableNames)
            ionTable = varargin{k};
        elseif ismember('ion',tableSort.Properties.VariableNames)&& ismember('color',tableSort.Properties.VariableNames)
            colorScheme = varargin{k};
        elseif ismember('abundance',tableSort.Properties.VariableNames)
            isotopeTable = varargin{k};
        else 
            error('please add colorScheme, isotopeTable and ionTable to your Input Variables');
        end
        continue
        % find interval
    elseif ~ischar(varargin{k})&& length(varargin{k}) == 1 && varargin{k} > 100000
        interval = varargin{k};
    % find binWidth
    elseif ~ischar(varargin{k})&& length(varargin{k}) == 1 && varargin{k}<1 
        binWidth = varargin{k};
    % check for axisROI
    elseif ~ischar(varargin{k})&& length(varargin{k})==4
        axisROI = varargin{k};

    elseif ischar(varargin{1,k})&& s==1
        sampleName = varargin{k};
                s = s+1;
    elseif ischar(varargin{1,k})&& s==2
        storageFolder = varargin{k};
    else
        error('Please check your input variables again');
    end
end
% check for raw pos file
if ismember('atom',posIn.Properties.VariableNames)
    posIn = posUnDecompose(posIn);
end

%% Basics for the calculation
maxNumSlice = round(height(posIn)/(interval/2));
startIon = 1;
endIon = startIon + (interval-1);

%% create a figure with the mass Spec

for i = 1:(maxNumSlice)
    
    % get the pos for the Spec 
    posSpec = posIn(startIon:endIon, :);
    
    % plot the pos
    specSlice = massSpecPlot(posSpec,binWidth, mode);
    title((sampleName + " massSpec slice " + i));
    
    % Add ions
    if exist('ionTable', 'var')
        for j=1:height(ionTable) 
            ionAdd(specSlice,string(ionTable.ionName(j)),ionTable.chargeState(j),isotopeTable,colorScheme,0,0.01,'most abundant',0.1);
        end
    end

    
    % change the size 
    if exist('axisROI', 'var')
        axis (axisROI)
    else
        axis ([0 100 10^-4 10^2])
    end
    f = gcf;
    f.Units = 'centimeters';
    f.Position = [30 20 15 12];
    
    

    % save the spec
    if exist('storageFolder', 'var')
        nameForSave = sampleName + "_massSpec_" + startIon + "_" + endIon;
        savefig(f, nameForSave)
        saveas(f,nameForSave,'tif');
    end
    % set the next start point
    startIon = endIon - round(interval/2);
%     startIon = endIon;
    endIon = startIon + (interval-1);
    % set the next start point
    if endIon >= height(posIn)
        endIon = height(posIn);
    end
    
    
end
end

