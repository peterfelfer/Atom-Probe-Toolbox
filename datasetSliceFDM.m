function totalFDM = datasetSliceFDM(posIn, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% INPUT
% posIn = posfile with detx and dety coordinates
% sampleName = Name of the Sample
% interval = number of ions for one interval
% range = start and end ion sequence number
% storageFolder = Folder where to store the FDM images as .fig and .tiff data
% the first character arry is defined as sampleName the secon as storage folder
%
% OUTPUT
% totalFDM = array with all the FDM
%% Basics - Define variables
s = 1;
interval = 2000000;
range = [1 height(posIn)];
sampleName = 'sample';


for k = 1:length(varargin)
    % find interval
    if ~ischar(varargin{k})&& length(varargin{k})==1
        interval = varargin{k};
    % find range
    elseif ~ischar(varargin{k})&& length(varargin{k})>1 
        range = varargin{k};
    % find sampleName
    elseif ischar(varargin{1,k})&& s==1
        sampleName = varargin{k};
        s = s+1;
    elseif ischar(varargin{1,k})&& s==2
        storageFolder = varargin{k};
    else
        error('Please check your input variables again');
        
    end
end

%% check for RAW pos file
if ismember('atom',posIn.Properties.VariableNames)
    posIn = posUnDecompose(posIn);
end

%% calculate for each slice the FDM image
maxNumSlice = round(height(posIn)/interval);
startIon = range(1);
endIon = startIon + interval;

for i = 1:maxNumSlice
    
    % get the pos for the FDM
    posFDM = posIn(startIon:endIon, :);
     
    f = figure;
    ax = axes;
    FDMresolution = 512; % number of pixels for the field desorption map -> 512
    [FDM, FDMcenters] = hist3(double([posFDM.detx posFDM.dety]),'Nbins',[FDMresolution FDMresolution]);
    FDMimage = imagesc(FDM,'XData',FDMcenters{1},'YData',FDMcenters{2});
    axis equal;
    axis tight;
    xlabel 'detector x [mm]'
    ylabel 'detector y [mm]'
    ax.YDir = 'normal';
    ax.XDir = 'normal';
    title(sampleName + " slice " + i);
    
    
    % save FDM for slider
    totalFDM{i} = FDM;
    
    if exist('storageFolder', 'var')
        % Save FDM in the Folder
        nameForSave = sampleName + "FDMslice" + startIon + "_" + endIon;
        savefig(gcf, nameForSave)
        saveas(gcf,nameForSave,'tif');
    end
    
    % set the next start point
    startIon = endIon - interval/2;
    endIon = startIon + (interval -1);
    % set the next start point
    if endIon >= height(posIn)
        endIon = height(posIn);
    end
    
end

end

