function files = hdf5FileFindByAttribute(topFolder,varargin)
% hdf5FileFindByAttribute gets a list of all hdf5 files the folder topFolder and all
% its subfolders and filters them by attributes. The attributes are parsed
% by combinations of keyword and value. You can check your hdf5 files for
% those combinations by using info = h5info('test.h5');
%
%   files = hdf5FileFindByAttribute(topFolder,varargin)
%
%   INPUTS:
%   topFolder       folder in which the search will be performed. All
%                   subfolders will be included
%
%   varargin        pairs of attributes and their values. For numeric
%                   values, comparions can be used. If comparisons are
%                   used, the value input needs to be char type
%                   examples:
%
%                   String type
%                   '/sample/materialType'  'aluminium'
%
%                   Numeric with comparison
%                   '/atomProbeTomography/experiment/specimenTemperature'
%                   '<= 40'
%
%                   Numeric with exact value
%                   '/atomProbeTomography/experiment/specimenTemperature'
%                   40
%
%   OUTPUTS:
%   files           files that fulfill the attributes as array of structs
%                   with the fields: name, folder, bytes, date, isdir,
%                   datenum
%


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
    % check if attributes are present in file and compare to value
    for att = 1:numAtt
        delimiter = find(attribute{att,1} == '/', 1, 'last' );
        try
            attributeVal = h5readatt([files(fi).folder '/' files(fi).name],attribute{att,1}(1:delimiter), attribute{att,1}(delimiter+1:end));
            
            if ischar(attributeVal) % for string type arguments
                isAttribute(fi,att) = strcmp(attributeVal,attribute{att,2});
                
            elseif isnumeric(attributeVal) % for numeric arguments
                if isnumeric(attribute{att,2})
                    isAttribute(fi,att) = attributeVal == attributeVal;
                elseif ischar(attribute{att,2})
                    eval(['isAttribute(fi,att) = ' num2str(attributeVal) attribute{att,2} ';']);
                end
            end
            
        catch
            isAttribute(fi,att) = false;
        end
    end
    
end

isAllAttributes = all(isAttribute,2);

files = files(isAllAttributes);