
function ions = apl_ionListFromRangeTable(rangeTable, extraIons)
% apl_ionListFromRangeTable builds a unique ion list from rangeTable.
% Returns table with columns: ionName, ionLabel, chargeState, ionTable

if nargin < 2
    extraIons = {};
end

ionName = strings(0,1);
ionLabel = strings(0,1);
chargeState = [];
ionTableCell = {};

% From rangeTable
if istable(rangeTable) && all(ismember({'rangeName','chargeState'}, rangeTable.Properties.VariableNames))
    for i = 1:height(rangeTable)
        name = string(rangeTable.rangeName(i));
        if name == "" || name == "not assigned" || name == "background" || name == "unknown"
            continue;
        end
        cs = rangeTable.chargeState(i);
        if isnan(cs) || cs == 0
            continue;
        end
        if ismember('ion', rangeTable.Properties.VariableNames)
            ionEntry = rangeTable.ion{i};
            if istable(ionEntry) && height(ionEntry) > 0
                ionTbl = ionEntry;
            else
                [ionTbl, ~] = ionConvertName(char(name));
            end
        else
            [ionTbl, ~] = ionConvertName(char(name));
        end
        nm = string(ionConvertName(ionTbl));
        lbl = string(ionConvertName(ionTbl, cs));
        ionName(end+1,1) = nm; %#ok<AGROW>
        ionLabel(end+1,1) = lbl; %#ok<AGROW>
        chargeState(end+1,1) = cs; %#ok<AGROW>
        ionTableCell{end+1,1} = ionTbl; %#ok<AGROW>
    end
end

% From extraIons list
if ~isempty(extraIons)
    for i = 1:numel(extraIons)
        [ionTbl, cs] = ionConvertName(char(extraIons{i}));
        nm = string(ionConvertName(ionTbl));
        lbl = string(ionConvertName(ionTbl, cs));
        ionName(end+1,1) = nm; %#ok<AGROW>
        ionLabel(end+1,1) = lbl; %#ok<AGROW>
        chargeState(end+1,1) = cs; %#ok<AGROW>
        ionTableCell{end+1,1} = ionTbl; %#ok<AGROW>
    end
end

ions = table(ionName, ionLabel, chargeState, ionTableCell, ...
    'VariableNames', {'ionName','ionLabel','chargeState','ionTable'});

% Remove duplicates by ionLabel
[~, ia] = unique(ions.ionLabel, 'stable');
ions = ions(ia,:);
end
