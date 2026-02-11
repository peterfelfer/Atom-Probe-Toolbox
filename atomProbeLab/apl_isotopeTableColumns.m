
function cols = apl_isotopeTableColumns(isotopeTable)
% apl_isotopeTableColumns detects column names in isotope table
%
% Returns struct with fields: element, isotope, mass, abundance

cols = struct('element','', 'isotope','', 'mass','', 'abundance','');

names = isotopeTable.Properties.VariableNames;

cols.element = pickCol(names, {"element","symbol","Element","Symbol"});
cols.isotope = pickCol(names, {"isotope","massnumber","mass_number","A","MassNumber"});
cols.mass = pickCol(names, {"weight","mass","atomicmass","Mass","AtomicMass"});
cols.abundance = pickCol(names, {"abundance","naturalabundance","relativeabundance","Abundance"});

if isempty(cols.element) || isempty(cols.isotope) || isempty(cols.mass) || isempty(cols.abundance)
    error('apl_isotopeTableColumns:missingColumns', 'Isotope table must contain element, isotope, mass/weight, and abundance columns.');
end

end

function col = pickCol(names, candidates)
col = '';
for i = 1:numel(candidates)
    idx = find(strcmpi(names, candidates{i}), 1, 'first');
    if ~isempty(idx)
        col = names{idx};
        return;
    end
end
end
