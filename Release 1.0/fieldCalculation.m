% function fieldstrength_calculation = fieldCalculation(element, chargestate, CSRparameters, pos)
% calculates the chargeStateratio of an element and his chargeStates in the
% pos, as well as the fieldstrength after the paper of Levi Tegg et. al,
% 'Estimation of the Electric Field in Atom Probe Tomography Experiments Using Charge State Ratios'
% Published 06 June 2024
%
%
% INPUT
% 
% element:      character with the name of the element of interesst
%
% CSRparameters:file with specific parameters depending on the Element and
% charge states
%
% chargeStates: double withe the chargeStates of interest
%
%
% pos:          decomposed pos file that contains ion and charge state of
%               the individual atoms
%
%function chargeStateRatioTable: necessary to calculate chargesateratio,
%               will be called by fieldstrength_calculation
%
%
% OUTPUT
% fieldstrength_calculation:    
%               table with the ratio to
%               the chargeStates and the calculated fieldstrengths
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nürnberg, written by Niklas
% Zimmermann
function fieldstrength_calculation = fieldCalculation(element, chargestate, CSRparameters, pos)
    % Initialize an empty table to store the results
    fieldstrength_calculation = table();

    % Display the element and charge state
    disp('Element:');
    disp(element);
    disp('Charge state:');
    disp(chargestate);
    used_element = element;
    used_chargestate = chargestate;

    % Convert chargestate from string to a numeric vector, if needed
    if ischar(chargestate) || isstring(chargestate)
        chargeStates = str2num(chargestate); % Konvertiert '[2,3]' zu [2, 3]
    else
        chargeStates = chargestate; % Wenn es bereits numerisch ist, keine Konvertierung
    end

   % Calculate the chargeStateRatioTable with the converted chargeStates
    chargeStateRatioTable = chargeStateRatioCalculate(element, pos, chargeStates);

    % Access values in row 1, columns 3 and 4 of chargeStateRatioTable
    x = chargeStateRatioTable{1, 3};
    y = chargeStateRatioTable{1, 4};
    x = x{1, 1}; % Extract value from the cell
    y = y{1, 1}; % Extract value from the cell

   % Calculate the CSR value
    CSR = y / x;
    fprintf('CSR-Wert: %.4f\n', CSR);

    % Ensure CSRparameters has at least 4 columns
    if size(CSRparameters, 2) < 4
        error('CSRparameters muss mindestens 4 Spalten haben.');
    end

   % Access CSRparameters columns
    element_check = CSRparameters{:, 1}; % First column
    chargestate_check = CSRparameters{:, 2}; % First column
    a_values = CSRparameters{:, 3}; % Third column (a values)
    b_values = CSRparameters{:, 4}; % Fourth column (b values)

    % Compare and calculate using CSRparameters
    found = false;
    row = 1; % Row counter for the table
    for i = 1:length(element_check)
        if element_check(i) == used_element && chargestate_check(i) == chargestate
            % Read values for a and b
            a = a_values(i);
            b = b_values(i);

            if isnan(a) || isnan(b)
                warning('Invalid numeric values found in row %d. Skipping calculation.', i);
                continue;
            end

            % Calculate Fieldstrength
            Fieldstrength = a * (1 - (b / ((CSR^0.3) + b + 0.256)));
            fprintf('Fieldstrength: %.4f\n', Fieldstrength);
            found = true;

            % Add result to the table
            fieldstrength_calculation(row, :) = {element, chargestate, CSR, Fieldstrength};
            row = row + 1; % Nächste Zeile der Tabelle
        end
    end

     % Check for matches
    if ~found
        disp('No matching values found.');
        fieldstrength_calculation = table(); % Return an empty table if no matches found
    else
        % Set column names
        fieldstrength_calculation.Properties.VariableNames = {'Element', 'ChargeState', 'CSR_Wert', 'Fieldstrength'};
    end
end