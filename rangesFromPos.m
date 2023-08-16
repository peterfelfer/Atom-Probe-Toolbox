function [rangeTable colorScheme] = rangesFromPos(pos, colorScheme, isotopeTable)
%rangesFromPos extracts and creates the corresponding rangeTable to an allocated pos variable
%
% rangeTable = rangesFromPos(pos, colorScheme)
%
% INPUT
% pos           pos Variable with allocated ranges
% colorScheme   colorScheme of the toolbox/project
%
% OUTPUT
% rangeTable    Table variable with all the ranges in the pos variable
%
%


%% check if pos is in Raw format

pos = posUnDecompose(pos);

%% extract range file information

% find unique ion and chargeState information 
rangeTable = table();
rangeTable = addvars(rangeTable, pos.ion, pos.chargeState, 'NewVariableNames',{'ion', 'chargeState'});
rangeTable = rmmissing(unique(rangeTable));

% Create rangeTable with zeros

rangeTable.Properties.VariableNames(1) = "rangeName";

vec = zeros(height(rangeTable),1);
test = cell(height(rangeTable),1);

% rangeTable = addvars(rangeTable, zeros, zeros, zeros, zeros, [zeros zeros zeros], 'NewVariableNames',{'mcbegin', 'mcend', 'volume', 'ion', 'color'});
rangeTable = addvars(rangeTable, vec, vec, vec, test, [vec vec vec], 'NewVariableNames',{'mcbegin', 'mcend', 'volume', 'ion', 'color'});



%% find Data

for i = 1 :height(ionsUnique)
    % find mcbegin and mcend
    posIon = pos(pos.ion == rangeTable.rangeName(i), :);
    posIon = posIon(posIon.chargeState == rangeTable.chargeState(i), :);

    rangeTable.mcbegin(i) = min(posIon.mc);
    rangeTable.mcend(i) = max(posIon.mc);

    % find ion
    %rangeTable.ion(i) = ionConvertName(string(rangeTable.rangeName(i)));
    outputTable = ionConvertName(string(rangeTable.rangeName(i)));

    % find isotope

            % find ion in ionList and extract chargeState and ion
            % potentialIons = ionList(ionList.ion == rangeName(isotopes),:);
            % if isempty(potentialIons)
            %     ionName = rangeName(isotopes);
            %     disp(['ion ' ionName ' is not on ion list. pls expand list to include ion']);
            %     isIn = 0;
            % else
            %     isIn = potentialIons.mc > mcbegin(isotopes) & potentialIons.mc < mcend(isotopes);
            %     if any(isIn)
            %         potentialIons = potentialIons(isIn,:);
            %         ion{isotopes,:} = ionConvertName(string(potentialIons.ionIsotopic(1)));
            %         chargeStateSingleIon = nnz(char(potentialIons.ionIsotopic(1)) == '+');
            %         chargeState(isotopes) = chargeStateSingleIon;
            %     elseif isIn == 0 
            %         mcbeginNew = (mcbegin(isotopes)-0.2);
            %         mcendNew = (mcend(isotopes)+0.2);
            % 
            %     else
            %         isIn = 0;
            %     end
            % end
    

    

    rangeTable.ion{i} = outputTable;

    % find color
    
    if colorScheme.ion == rangeTable.rangeName(i)
        rangeTable.color(i,:) = colorScheme.color(colorScheme.ion == rangeTable.rangeName(i), :); 
    else
        colorScheme = colorSchemeIonAdd(colorScheme,string(rangeTable.rangeName(i)));
        rangeTable.color(i,:) = colorScheme.color(colorScheme.ion == rangeTable.rangeName(i), :); 
    end

   

end










end