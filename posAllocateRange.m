function pos = posAllocateRange(pos,rng,options)
% posAllocateRange takes a pos and a range variable and
% allocates ion hits to it.
%
% INPUT
% pos:      table of reconstructed atom positions with ionIdx, x, y, z, and m/c
%
% rng:      table with extracted ranges from the mass spec
%           input possible as rangesExtractFromMassSpec(spec)
%
% options:  'decompose': the complex ions will be split up
%           'raw': complex ions will not be split up
%
% OUTPUT
% pos:      first case, with 'raw' as options input used: table of reconstructed atom
%           positions with ionIdx, x, y, z, m/c, and additional ion and
%           chargeState fields
%           second case, with 'decompose' as options input used: table of
%           reconstructed atom positions with ionIdx, x, y, z, m/c, and
%           additional ion, chargeState, atom, isotope, and ionComplexity fields

% find the range the ion is in
in = (pos.mc > rng.mcbegin') & (pos.mc < rng.mcend');
unranged = ~sum(in,2);

in = [in unranged];

[rngIdx, ~] = find(in');

if ~(length(rngIdx) == height(pos)) % if ions are allocated to more than one range
    error('overlapping ranges');
end

% extract ion names and append for unranged
ionNames = rng.rangeName;
ionNames = [ionNames; categorical(string(missing))];

% extract charge state
chargeStates = rng.chargeState;
chargeStates = [chargeStates; NaN];

% allocation by range
if strcmp(options,'raw')
    % allocate ion name to pos
    pos.ion = ionNames(rngIdx);
    
    % allocate charge state to pos variable
    pos.chargeState = chargeStates(rngIdx);
    
    
    
    % allocation with decomposition
elseif strcmp(options,'decompose')
    numIon = height(pos);
    for r = 1:height(rng)
        rngComplexity(r,:) = height(rng.ion{r});
        rngElements(r,1:rngComplexity(r)) = rng.ion{r}.element';
        rngIsotopes(r,1:rngComplexity(r)) = rng.ion{r}.isotope';
    end
    rngComplexity(end+1,:) = 1; % for unranged ions
    rngElements(end+1,1) = categorical(string(missing));
    rngIsotopes(end+1,1) = categorical(string(missing));
    
    % create vector with indices mapping into new table
    ionComplexity = rngComplexity(rngIdx); % complexity of each hit
    mapVec = (repelem(1:height(pos),ionComplexity))'; % map vector to decomposed pos variable
    numAtomDecomp = length(mapVec);
    
    % create vectors with ion name, chargeState, element, isotope,
    % ion complexity
    
    % find index of atom in ion table
    firstIdx = cumsum(ionComplexity); % index of first occurrence of ion
    atomIdx = (1:numAtomDecomp)' + ionComplexity(mapVec) - firstIdx(mapVec); % index of atom in ion table
    
    % index into rngElements and rngIsotopes linearly
    linIdx = sub2ind(size(rngElements),rngIdx(mapVec),atomIdx);
    
    
    pos = pos(mapVec,:);
    pos.ion = ionNames(rngIdx(mapVec));
    pos.chargeState = chargeStates(rngIdx(mapVec));
    pos.atom = rngElements(linIdx);
    pos.isotope = rngIsotopes(linIdx);
    
    % give all unranged and unknown ions an ion complexity of NaN
    pos.ionComplexity = ionComplexity(mapVec);
    pos.ionComplexity(isundefined(pos.ion) | isnan(pos.chargeState)) = NaN;
    
    % return notification if options input is not 'decomposed'
elseif ~strcmp(options,'decompose')
    
    error('Invalid options input. Only allowed option is: ''decompose''.');
    
end





