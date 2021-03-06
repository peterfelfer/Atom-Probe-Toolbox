function colorScheme = colorSchemeIonAdd(colorScheme, newIon, selection)
% colorSchemeIonAdd adds new ions to the colorScheme table and adds a new
% color which was not used before.
% 
% colorScheme = colorSchemeIonAdd(colorScheme, newIon, select)
% colorScheme = colorSchemeIonAdd(colorScheme, newIon)
%
% INPUT
% colorScheme:  The existing colorScheme
% newIon:       Name of the newIon, if the ion already exists, the function will
%               end and the command window shows "ion already exists in
%               colorScheme"
% selection:    'select' if you want to choose the color
%               if nothing is parsed, the color is randomly generated
%               'create' a new color by picking the one that is the farthest 
%               from all existing colors
% 
% OUTPUT
% colorScheme:  the new colorScheme with the new ion added at the end of
%               the table
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-N�rnberg

%% number of random colors to be generated for farthest point sampling
NCOLS = 10000; 

%% Check for ion name
ionTable = ionConvertName(newIon); % Convert the name in the table format
ionName = ionConvertName(ionTable); % make the correct name

length = size(colorScheme.ion);
ionAlreadyExist = false;
%% Check if ion already exists in colorScheme

for j = 1:length(1,1)
        if colorScheme.ion(j,1) == ionName(1,:)
            disp ('ion already exists in colorScheme')
            return % don't start the color searching loop 
        else 
            match = true; 
            ionAlreadyExist = false;
        end
end
        


%% check if color already exists for another ion, if the color already exists,
%   you need to choose a new color 
n = false;
while match == true
    % add a color to the new ion

   if ~exist('selection','var')
       selection = 'create';
   end
       
       if strcmp(selection,'select') && ionAlreadyExist == false 
                if n == false
                    color = uisetcolor();
                    n = true;
                else 
                    title = 'color already exists';
                    color = uisetcolor(title);
                end

       elseif strcmp(selection,'create')
                % generate new color from NCOLS number of randomly generated colors,
                % by picking the one that is the farthest from all existing colors
                cols = rand(NCOLS,3);
                dist = pdist2(colorScheme.color,cols,"euclidean");
                md = min(dist,[],1);
                [~, idx] = max(md);
                color = cols(idx,:);
           
        end
    
    % check if color already exists and change the match variable 
    for i = 1:length(1,1)
            if colorScheme.color(i,1) == color(1,1) && colorScheme.color(i,2) == color(1,2) && colorScheme.color(i,3) == color(1,3)
                match = true; % color already exist in colorScheme
                n = true; % countnumber for the second round of choosing a color 
                i = length(1,1);
            else
                match = false; % color does not exist in colorScheme
            end


    end
end



%% write ionName and color in the colorScheme
if ionAlreadyExist == false
    newIonRow = {ionName, color};
    colorScheme(end+1,:) = newIonRow;
end

end