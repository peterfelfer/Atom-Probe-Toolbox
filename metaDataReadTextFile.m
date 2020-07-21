function meta = metaDataReadTextFile(fileStr)
% metaDataReadTextFile reads a text file with standard UTF-8 encoding into
% a matlab struct. Information in the metadata file is the general form of:
%   variableName [format] = value [unit]
% comments are denoted by a % at the beginning of the line.
% the string expression 'NULL' denotes an undefined value, resulting in an
% empty cell in the output cell arrray. Case insensitive!
%
% meta = metaDataReadTextFile(fileStr)
%
% INPUT 
% fileStr:      string content of the *.metadata file containg the 
%               information
%
% OUTPUT
% meta:         cell array with name - value combinations in the form of:
%               {name, value, unit}


% pre-formatting of the text file
%break it into lines
fileByLine = regexp(fileStr, '\n', 'split');
fileByLine = fileByLine';
%remove empty lines
fileByLine( cellfun(@isempty,fileByLine) ) = [];
% remove comments
fileByLine( cellfun(@(x) x(1) == '%',fileByLine )) = [];

%% going through each line of metadata
meta = {};
for li = 1:length(fileByLine)
    isBrackets = find(fileByLine{li} == '[' | fileByLine{li} == ']');
    isEqualSign = find(fileByLine{li} == '=');
    
    %variable format: left of =, in brackets
    varFormat = fileByLine{li}(isBrackets(1)+1:isBrackets(2)-1);
    
    %variable name
    meta{li,1} = strtrim(fileByLine{li}(1:isBrackets(1)-1));
    
    %variable value
     varStr = strtrim(fileByLine{li}(isEqualSign+1:isBrackets(3)-1));
     
     
     %units
     meta{li,3} = strtrim(fileByLine{li}(isBrackets(3)+1:isBrackets(4)-1));

     
     if strcmp(upper(varStr),'NULL')
        meta{li,2} = {};
        
     elseif strcmp(varFormat,'string')
         meta{li,2} = varStr;
         
     elseif strcmp(varFormat,'dateTime')
         meta{li,2} = datetime(varStr,'InputFormat',strtrim(fileByLine{li}(isBrackets(3)+1:isBrackets(4)-1)));
         
     elseif strcmp(varFormat,'enum')
         cat = categorical(strtrim(split(meta{li,3},',')));
         meta{li,2} = cat(cat == varStr);
         
     elseif strcmp(varFormat,'bool')
         meta{li,2} = str2num(varStr);
         
     else % numeric formats
         meta{li,2} = str2num(varStr);
     end
    
    %variable unit
    meta{li,3} = strtrim(fileByLine{li}(isBrackets(3)+1:isBrackets(4)-1));
    
end
