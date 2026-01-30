function conc = posCalculateConcentrationSimple(pos, detEff, excludeList, volumeName, varargin)
% posCalculateConcentrationSimple calculates the concentration of a
% categorical list of atoms or ions
%
% conc = posCalculateConcentrationSimple(pos, detEff, excludeList, volumeName);
% conc = posCalculateConcentrationSimple(pos, detEff, excludeList)
% conc = posCalculateConcentrationSimple(pos, detEff)
% conc = posCalculateConcentrationSimple(..., 'mode', mode)
%
% INPUT
% pos:          decomposed pos file that contains ion and charge state of
%               the individual atoms
%
% detEff:       detector efficiency of the atom probe, can be parsed as
%               or as a fraction (for a LEAP 4000X HR it is 0.37)
%
% excludeList:  cell array that contains as character the individual
%               ions that shall not be considered for the concentration
%               calculation, unranged atoms appear as 'unranged', if not
%               parsed, no atoms will be excluded
%
% volumeName:   name of the volume, parsed as character array, if not parsed,
%               the volume will be named after pos
%
% Name-Value Options:
% 'mode'        Concentration calculation mode:
%               'ionic' - each unique ion species counted separately
%                         (56Fe++ and 56Fe+ are different)
%               'isotopic' - groups by isotopic composition, ignores charge
%                         (56Fe++ and 56Fe+ counted together as 56Fe)
%               'atomic' - groups by element, ignores isotope and charge
%                         (all Fe isotopes counted together)
%               Default: 'atomic' if 'atom' column exists, 'ionic' otherwise
%
% OUTPUT
% conc:         is a table that contains the count, concentration, and
%               variance for each atom/ion that is not on the excludeList.
%               statistical deviation calculated after Danoix et al.,
%               https://doi.org/10.1016/j.ultramic.2007.02.005
%               variance(conc) = conc*(1-conc)/numAtomsDetected * (1 - detEff)
%
% EXAMPLES
% % Ionic concentration (each charge state separate)
% conc = posCalculateConcentrationSimple(pos, 0.37, {}, '', 'mode', 'ionic');
%
% % Isotopic concentration (charge states grouped, isotopes separate)
% conc = posCalculateConcentrationSimple(pos, 0.37, {}, '', 'mode', 'isotopic');
%
% % Atomic concentration (all isotopes and charge states grouped)
% conc = posCalculateConcentrationSimple(pos, 0.37, {}, '', 'mode', 'atomic');
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

%% Parse name-value options
mode = '';  % empty means auto-detect
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        if strcmpi(varargin{k}, 'mode')
            mode = lower(string(varargin{k+1}));
        end
    end
end

%% detector efficiency
if detEff > 1
    detEff = detEff/100;
end

%% Determine concentration mode and get appropriate data
% Priority: explicit mode > auto-detect from columns
hasAtomColumn = any(ismember(pos.Properties.VariableNames, 'atom'));
hasIonColumn = any(ismember(pos.Properties.VariableNames, 'ion'));

if isempty(mode) || mode == ""
    % Auto-detect: use atom column if available, otherwise ion
    if hasAtomColumn
        type = 'atomic';
        columnType = 'atom';
        atoms = pos.atom;
    else
        type = 'ionic';
        columnType = 'ion';
        atoms = pos.ion;
    end
else
    % Explicit mode specified
    switch mode
        case "ionic"
            if ~hasIonColumn
                error('posCalculateConcentrationSimple:missingColumn', ...
                    'Ionic mode requires ''ion'' column in pos table.');
            end
            type = 'ionic';
            columnType = 'ion';
            atoms = pos.ion;

        case "isotopic"
            if ~hasIonColumn
                error('posCalculateConcentrationSimple:missingColumn', ...
                    'Isotopic mode requires ''ion'' column in pos table.');
            end
            type = 'isotopic';
            columnType = 'ion';
            % Convert ionic names to isotopic (strip charge state)
            atoms = ionConvertMode(pos.ion, 'isotopic');

        case "atomic"
            if hasAtomColumn
                % Use existing atom column
                type = 'atomic';
                columnType = 'atom';
                atoms = pos.atom;
            elseif hasIonColumn
                % Convert from ion column
                type = 'atomic';
                columnType = 'ion';
                atoms = ionConvertMode(pos.ion, 'atomic');
            else
                error('posCalculateConcentrationSimple:missingColumn', ...
                    'Atomic mode requires ''atom'' or ''ion'' column in pos table.');
            end

        otherwise
            error('posCalculateConcentrationSimple:invalidMode', ...
                'Mode must be ''ionic'', ''isotopic'', or ''atomic''.');
    end
end

if ~exist('volumeName', 'var') || isempty(volumeName)
    volumeName = inputname(1);
end

% need to assign, otherwise not counted in countcats
atoms(isundefined(atoms)) = 'unranged';

%% check for excluded types
cats = categories(atoms);
if exist('excludeList','var')
    isExcluded = ismember(cats,excludeList);
else
    isExcluded = false(size(cats));
end
isExcluded = isExcluded';

%% calculate concentrations for not excluded variables
counts = countcats(atoms);
if iscolumn(counts)
    counts = counts';
end
counts(2,:) = counts./sum(counts(~isExcluded));
counts(2,isExcluded) = 0;
counts(3,:) = counts(2,:).*(1- counts(2,:))./counts(1,:) * (1-detEff);



%% creating output table
conc = array2table(counts,'VariableNames',cats');
conc.Properties.VariableDescriptions = repmat({columnType},size(cats'));

conc = [table(categorical({volumeName;volumeName;volumeName}),'VariableNames',{'volume'})...
    table([0;0;0],'VariableNames',{'distance'} )...
    table(categorical({type;type;type}),'VariableNames',{'type'}),...
    table(categorical({'count'; 'concentration';'variance'}),'VariableNames',{'format'}),...
    conc];

