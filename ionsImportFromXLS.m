function [ions, chargeStates] = ionsImportFromXLS(fileName)
% ionsImportFromXLS Imports data of ions names and charge states from
% a Microsoft Excel spreadsheet file named fileName.
%
% [ions, chargeStates] = ionsImportFromXLS(fileName)
%
% INPUT       
% fileName:        file name and path of excel spreadsheet with 1st column
%                  being the ion name, subsequent columns for charge states.
%                  Since the charge states are read as boolean values, the
%                  entry for an existing charge state must be 'true'.
%                  Nonexisitng charge states need not be denoted as 'false'.
%
% OUTPUT
% ions:            string vector of all ion names in the Excel file
%
% chargeStates:    NxM logical array with 'true' for any charge state that has
%                  been observed in the literature. N is the same length as
%                  ions, M is 6 (may be expanded if chargeStates >6 are
%                  observed experimentally)

%% Input handling

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = 1;
opts.DataRange = "A:G";

% Specify column names and types
opts.VariableNames = ["ion", "+", "++", "+++", "++++", "+++++", "++++++"];
opts.VariableTypes = ["string", "logical", "logical", "logical", "logical", "logical", "logical"];

% Specify variable properties
opts = setvaropts(opts, "ion", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "ion", "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["+", "++", "+++", "++++", "+++++", "++++++"], "FillValue", false);
opts = setvaropts(opts, ["+", "++", "+++", "++++", "+++++", "++++++"], "TreatAsMissing", '');

% Import the data
ionTable = readtable(fileName, opts, "UseExcel", false);
ionTable = ionTable(2:end,:);

ions = ionTable.ion;

chargeStates = table2array(ionTable(:,2:end));

