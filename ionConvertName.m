function varargout = ionConvertName(varargin)
% ionConvertNameTable converts an ion name to a table variable and reverse.
% As input parameters, also categorical or an array of elements is
% possible.
%
% Create ionTable
% ionTable = ionConvertName(ionName);
% [ionTable chargeState] = ionConvertName(ionName);
%
% Create ionName
% table:
%   ionName = ionConvertName(ionTable,chargeState,format);
%   ionName = ionConvertName(ionTable,chargeState);
%   ionName = ionConvertName(ionTable,NaN,format);
%   ionName = ionConvertName(ionTable);
% categorical:
%   ionName = ionConvertName(ionCategorical,chargeState,format);
%   ionName = ionConvertName(ionCategorical,chargeState);
%   ionName = ionConvertName(ionCategorical)
% array:
%   ionName = ionConvertName(ionArray,chargeState,format,isotopeTable);
%
% INPUT/OUTPUT:
%
% ionName: Name of the ion
%   (isotope element count) x N chargestate        '56Fe2 16O3 ++' 
%   (element count) x N chargestate                'Fe2 O3 ++'
%   (element count) x N                            'Fe2 O3'
%   individual nucleides will be sorted by atomic number descending 
%   e.g. 'O H2'
% ionTable:         Table that contains the element and the isotope
% ionCategorical:   ion table that is stored as an categorical array
% ionArray:         Matrix of element and isotope number 
%                   e.g. Fe2O3 = [26,56;26,56;8,16;8,16;8,16;];
% chargeState:      is the charge State of the ion 
% NaN:              if no chargeState is parsed
% format:           can be 'plain' or 'LaTeX'
% isotopeTable:     Table with all isotopes

%% conversion from array to table. The new table is saved as an input 
%% argument and is converted to ionName as a normal table

% check for input as an array 
if ismatrix(varargin{1}) & ~istable(varargin{1}) & ~iscategorical(varargin{1}) & ~ischar(varargin{1})
   %check for isotopeTable as an input variable
   if nargin == 4 & istable(varargin{4})
    ionArray = varargin{1};
    isotopeTable = varargin{4};
    % create element List
    for j = 1:length(ionArray)
        isotopeName = isotopeTable.element(isotopeTable.atomicNumber==ionArray(j,1));
        element(j,1) = isotopeName(1,1);
    end
    % create isotope List
    isotope = ionArray(:,2);
    % create ionTable with element and isotope
    ionTable = table(element, isotope);
    varargin{1} = ionTable;
   else
       error('no isotopeTable is parsed');
   end
end

%% conversion from table to name
if istable(varargin{1})
    ionTable = varargin{1};
    if nargin == 3
        format = varargin{3};
    else
        format = 'plain';
    end
    
    % add atomic number to table
    numAtom = height(ionTable);
    for at = 1:numAtom
        atomicNumber(at,:) = symbolConvertAtomicNumber(char(ionTable.element(at)));
    end
    % for unknown ions, NaNs are swapped for 0, to enable sorting
    atomicNumber(isnan(atomicNumber)) = 0;
    
    ionTable = addvars(ionTable, atomicNumber,'NewVariableNames','atomicNumber');
    
    % sort by atomic number descending
    ionTable = sortrows(ionTable,{'atomicNumber','isotope'},{'descend','descend'},'MissingPlacement','first');
    
    % get multiplicity of atom occurances in ion and write out names
    ionName = [];
    
    
    isotopeGroup = sortedFindgroups(ionTable);
    for i = 1:max(isotopeGroup)
        idx = find(isotopeGroup == i,1);
        % add isotope to the name
        if ~isnan(ionTable.isotope(idx))
            if strcmp(format,'LaTeX')
                ionName = [ionName '^{' num2str(ionTable.isotope(idx)) '}'];
            else
                ionName = [ionName num2str(ionTable.isotope(idx))];
            end
        end
        
        % add chemical element to the name
        ionName = [ionName char(ionTable.element(idx))];
        
        % add chemical element to the name
        if sum(isotopeGroup == i) > 1
            if strcmp(format,'LaTeX')
                ionName = [ionName '_{' num2str(sum(isotopeGroup == i)) '}'];
            else
                ionName = [ionName num2str(sum(isotopeGroup == i))];
            end
        end
        ionName = [ionName ' '];
    end
    
    % add + (or -) for chargestates to the name
    if nargin > 1
        chargeState = varargin{2};
        if ~isnan(chargeState) %NaN for undefined chargestate, e.g. in noise
            
            if chargeState < 0
                sym = '-';
            else
                sym = '+';
            end
            if strcmp(format,'LaTeX')
                ionName = [ionName ' ^{' repmat(sym,1,abs(chargeState)) '}'];
            else
                ionName = [ionName repmat(sym,1,abs(chargeState))];
            end
        end
    end
    
    
    varargout{1} = strtrim(ionName);
end



%% conversion from categorical to name
if iscategorical(varargin{1})
    element = varargin{1};
    
    if nargin == 3
        format = varargin{3};
    else
        format = 'plain';
    end
    
    numAtom = length(element);
    for at = 1:numAtom
        atomicNumber(at,:) = symbolConvertAtomicNumber(char(element(at)));
    end
    ionTable = table(element,atomicNumber);
    
    % sort by atomic number descending
    ionTable = sortrows(ionTable,{'atomicNumber'},{'descend'});
    
    % get multiplicity of atom occurances in ion and write out names
    ionName = [];
    isotopeGroup = sortedFindgroups(ionTable);
    for i = 1:max(isotopeGroup)
        idx = find(isotopeGroup == i,1);
        ionName = [ionName ' ' char(ionTable.element(idx))];
        
        if sum(isotopeGroup == i) > 1
            if strcmp(format,'LaTeX')
                ionName = [ionName '_{' num2str(sum(isotopeGroup == i)) '}'];
            else
                ionName = [ionName num2str(sum(isotopeGroup == i))];
            end
        end
    end
    
    
    
    % add + (or -) for chargestates
    if nargin > 1
        chargeState = varargin{2};
        if ~isnan(chargeState) %NaN for undefined chargestate, e.g. in noise
            if chargeState < 0
                sym = '-';
            else
                sym = '+';
            end
            if strcmp(format,'LaTeX')
                ionName = [ionName ' ^{' repmat(sym,1,abs(chargeState)) '}'];
            else
                ionName = [ionName repmat(sym,1,abs(chargeState))];
            end
        end
    end
    
    
    varargout{1} = strtrim(ionName);
    
end




%% conversion from name to table
if ischar(varargin{1})
    ionName = varargin{1};
    
    %check for chargestate
    if any(ionName == '+')
        chargeState = sum(ionName == '+');
        ionName(ionName == '+') = [];
    elseif any(ionName == '-')
        chargeState = uminus(sum(ionName == '-'));
        ionName(ionName == '-') = [];
    else
        chargeState = NaN;
    end
    
    ionName = strtrim(ionName); %remove any whitespace
    
    % split individual elemental / isotopic parts
    parts = strsplit(ionName);
    numElements = length(parts);
    isotope = [];
    elementOut = categorical();
    isotopeOut = [];
    
    for el = 1:numElements
        % find the chemical element
        element = parts{el}(isstrprop(parts{el},'alpha'));
        
        isDigit = isstrprop(parts{el},'digit');
        % find the isotope
        if isDigit(1) == true
            isotope = str2num(parts{el}(1:find(~isDigit,1)-1));
        else
            isotope = NaN;
        end
        
        % find the count
        if isDigit(end) == true
            count = str2num(parts{el}(find(~isDigit,1,'last')+1:end));
        else
            count = 1;
        end
        
        for cnt = 1:count
            elementOut(end+1) = element;
            isotopeOut(end+1) = isotope;
        end
        
    end
    element = elementOut';
    isotope = isotopeOut';
    
    varargout{1} = table(element,isotope);
    varargout{2} = chargeState;
    
end

function group = sortedFindgroups(ionTable)
% replace all NaN isotopes with 0, for findgroups command
% this is hacky, but we get NaN outputs if any entry in the table is NaN
if any(string(ionTable.Properties.VariableNames) == "isotope")
    ionTable.isotope(isnan(ionTable.isotope)) = 0;
end

isotopeGroup = findgroups(ionTable); % need to re-sort REALLY HOPE THAT WONT MAKE PROBLEMS
groupIdx = unique(isotopeGroup,'stable');
[~, idx] = sort(groupIdx);
for i = 1:length(isotopeGroup)
    isotopeGroup(i) = idx(isotopeGroup(i));
end

%replace 0s with NaNs again
if any(string(ionTable.Properties.VariableNames) == "isotope")
    ionTable.isotope(ionTable.isotope == 0) = NaN;
end
group = isotopeGroup;







