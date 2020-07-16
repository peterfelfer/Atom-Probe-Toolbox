function files = hdf5findFileByAttribute(topFolder,varargin)
% HDF5deepDrill gets a list of all hdf5 files the folder topFolder and all
% its subfolders and filters them by attributes. The attributes are parsed
% by combinations of keyword and value. You can check your hdf5 files for
% those combinations by using info = h5info('test.h5');

% test if number of input argument is viable
if mod((nargin-1),2) == 1
    error('use only attribute - value combinations');
end

numAtt = (nargin-1)/2;

% get attribute value combinations from input
for att = 1:numAtt
    attribute{att,1} = varargin{att*2-1};
    attribute{att,2} = varargin{att*2};
end

% get all hdf5 files in folder structure
files = dir([topFolder '/**/*.h5']);

for fi = 1:length(files)
    % get hdf5 file info
    info = h5info([files(fi).folder files(fi).name]);
    
    % check if attributes are present in file and compare to value
    for att = 1:numAtt
        XXXXXXXXXXXX
    end
    
end