function w = ionWeight(ion, isotopeTable, chargeState)
% ionWeight calculates the weight of an ion, based on the provided isotope table
%
% w = ionWeight(ion, isotopeTable);
% w = ionWeight(ion, isotopeTable, chargeState);
%
% INPUTS
% ion:          the definition of the ion as a table with ion.element 
%               categorical vector of chemical element)and 
%               ion.isotope (int vector of isotope number)
%
% isotopeTable: table of all isotopes from APT Toolbox database
%
% chargeState:  int of charge state of the ion, optional
%
% OUTPUTS
% w: weight of ion in amu (w/o chargeState) or Da

%% extract the individual ions and get weight

w = 0;
for i = 1:height(ion)             % sums up the weight in case of complex/molecular ions
    w = w + isotopeTable.weight(isotopeTable.element == ion.element(i) ...     
        & isotopeTable.isotope == ion.isotope(i));
end

%% divide by chargestate if applicable
if exist('chargeState','var')
    w = w/chargeState;
end
