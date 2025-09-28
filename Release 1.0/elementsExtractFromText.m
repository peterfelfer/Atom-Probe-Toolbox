function elements = elementsExtractFromText(text)
% elementsExtractFromText searches for elemental symbols in a text variable
% and outputs them as a string array. This can be used to create ion lists.
% For this the Element needs to be present as a one or two element
% character only with whitespace or special characters around it (colons etc..).
%
% elements = elementsExtractFromText(text)
%
% INPUT
% text =  text document with the element as a one or two element character 
%         only with whitespace or special characters around it
%
% OUTPUT
% elements = string array with the elements
%%
text = string(text);

allElements = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",...
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",...
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",...
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",...
    "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",...
    "Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",...
    "Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"]';


isContained = false(size(allElements));


for e = 1:length(allElements)
    occurances = strfind(text, allElements(e));
    elementStringLength = length(char(allElements(e)));
    
    if occurances
        for o = 1:length(occurances)
            
            %check if letter combo ist jsut part of a word
            if occurances(o) == 1
                boundaryCharacters = text{1}(2);
            elseif occurances(o) == length(text{1}) - elementStringLength + 1
                boundaryCharacters = text{1}(occurances(o)-1);
                
            else
                boundaryCharacters = [text{1}(occurances(o)-1) text{1}(occurances(o)+elementStringLength)];
            end
            
            if any(isletter(boundaryCharacters)) & not(isContained(e,1))
                isContained(e,1) = false;
            else
                isContained(e,1) = true;
            end
        end
    end
end

elements = allElements(isContained);