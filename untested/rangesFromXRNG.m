function ranges = rangesFromXRNG(filename, isotopeTable)
% RANGESFROMXRNG Parse an xrng (XML range) file into a range table.
%
% ranges = rangesFromXRNG(filename, isotopeTable)
%
% INPUT
%   filename     - path to an .xrng file
%   isotopeTable - isotope table (from isotopeTable_naturalAbundances.mat)
%
% OUTPUT
%   ranges - table with columns:
%       mcbegin      - range start (Da)
%       mcend        - range end (Da)
%       ion          - ion name string (e.g. 'Fe', 'Ti O')
%       chargeState  - charge state
%       volume       - atomic volume (nm^3)
%
% (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

% Read and clean XML — some xrng files have malformed tags (e.g. unclosed
% <item> elements). Fix by removing incomplete range entries before parsing.
xmlText = fileread(filename);

% Remove <isotopiccombinations> blocks — they contain nested <item> tags
% that confuse the range parser, and we don't need them.
xmlText = regexprep(xmlText, '<isotopiccombinations>.*?</isotopiccombinations>', '', 'dotall');

% Remove any <item> that is opened but not closed before </ranges>
xmlText = regexprep(xmlText, '<item>\s*</ranges>', '</ranges>');

% Fix duplicate <item><item> (missing closing tag between entries)
xmlText = regexprep(xmlText, '</item>\s*<item>\s*<item>', '</item>\n      <item>');

% Parse the cleaned XML from string
tmpFile = [tempname '.xrng'];
fid = fopen(tmpFile, 'w');
fwrite(fid, xmlText);
fclose(fid);
try
    doc = xmlread(tmpFile);
catch ME
    delete(tmpFile);
    rethrow(ME);
end
delete(tmpFile);

rangesNode = doc.getElementsByTagName('ranges').item(0);
rangeItems = rangesNode.getChildNodes();

mcbegin = [];
mcend = [];
ionNames = {};
chargeStates = [];
volumes = [];

for i = 0:rangeItems.getLength()-1
    node = rangeItems.item(i);
    if ~strcmp(char(node.getNodeName()), 'item'); continue; end

    % Get mcbegin, mcend
    mcBNode = node.getElementsByTagName('masstochargebegin');
    mcENode = node.getElementsByTagName('masstochargeend');
    if mcBNode.getLength() == 0 || mcENode.getLength() == 0; continue; end

    mcB = str2double(char(mcBNode.item(0).getTextContent()));
    mcE = str2double(char(mcENode.item(0).getTextContent()));

    % Get ion info from first <ions><item>
    ionsNode = node.getElementsByTagName('ions');
    if ionsNode.getLength() == 0; continue; end

    % Find the first <item> directly under <ions>
    ionItemNodes = ionsNode.item(0).getChildNodes();
    ionNode = [];
    for j = 0:ionItemNodes.getLength()-1
        if strcmp(char(ionItemNodes.item(j).getNodeName()), 'item')
            ionNode = ionItemNodes.item(j);
            break;
        end
    end
    if isempty(ionNode); continue; end

    csNode = ionNode.getElementsByTagName('chargestate');
    cs = str2double(char(csNode.item(0).getTextContent()));

    volNode = ionNode.getElementsByTagName('volume');
    vol = str2double(char(volNode.item(0).getTextContent()));

    % Get atomic numbers from <atoms><item><atomicnumber>
    atomsNode = ionNode.getElementsByTagName('atoms');
    if atomsNode.getLength() == 0; continue; end
    atomItemNodes = atomsNode.item(0).getChildNodes();
    atomicNums = [];
    for a = 0:atomItemNodes.getLength()-1
        aNode = atomItemNodes.item(a);
        if strcmp(char(aNode.getNodeName()), 'item')
            anNode = aNode.getElementsByTagName('atomicnumber');
            if anNode.getLength() > 0
                atomicNums(end+1) = str2double(char(anNode.item(0).getTextContent()));
            end
        end
    end

    if isempty(atomicNums); continue; end

    % Convert atomic numbers to ion name using ionConvertName
    ionName = ionConvertName(atomicNums(:), NaN, 'plain', isotopeTable);

    mcbegin(end+1,1) = mcB;
    mcend(end+1,1) = mcE;
    ionNames{end+1,1} = ionName;
    chargeStates(end+1,1) = cs;
    volumes(end+1,1) = vol;
end

ion = ionNames;
chargeState = chargeStates;
volume = volumes;
ranges = table(mcbegin, mcend, ion, chargeState, volume);

end
