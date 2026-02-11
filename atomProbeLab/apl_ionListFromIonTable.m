
function ions = apl_ionListFromIonTable(ionTable, extraIons)
% apl_ionListFromIonTable builds a unique ion list from ionTable
% ionTable from ionsExtractFromMassSpec (ionName, chargeState, ion, color, isTracer)

if nargin < 2
    extraIons = {};
end

ionName = strings(0,1);
ionLabel = strings(0,1);
chargeState = [];
ionTableCell = {};

if istable(ionTable)
    for i = 1:height(ionTable)
        cs = ionTable.chargeState(i);
        if isnan(cs) || cs == 0
            continue;
        end
        if ismember('ion', ionTable.Properties.VariableNames)
            ionEntry = ionTable.ion{i};
            if istable(ionEntry)
                ionTbl = ionEntry;
            else
                % ionEntry may be categorical array of elements
                ionTbl = table(categorical(ionEntry), NaN(numel(ionEntry),1), 'VariableNames', {'element','isotope'});
            end
        else
            [ionTbl, ~] = ionConvertName(char(ionTable.ionName(i)));
        end
        nm = string(ionConvertName(ionTbl));
        lbl = string(ionConvertName(ionTbl, cs));
        ionName(end+1,1) = nm; %#ok<AGROW>
        ionLabel(end+1,1) = lbl; %#ok<AGROW>
        chargeState(end+1,1) = cs; %#ok<AGROW>
        ionTableCell{end+1,1} = ionTbl; %#ok<AGROW>
    end
end

% extra ions
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

[~, ia] = unique(ions.ionLabel, 'stable');
ions = ions(ia,:);
end
