function pos = vestaXYZimport(fileName)

% Select an xyz file
[filename, pathname] = uigetfile('*.xyz', 'Select an XYZ file');
filepath = fullfile(pathname, filename);

% Check if user canceled the file selection
if isequal(filename,0)
    disp('User selected Cancel');
    return;
end

% Open and read the file
fid = fopen(filepath, 'r');

% Read the number of entries (first line) and name (second line)
nEntries = str2double(fgetl(fid));
name = fgetl(fid);

% Pre-allocate the cell array to store data
data = cell(nEntries, 4); % 4 columns: Atom name, x, y, z

% Loop through the file and read atom names and coordinates
for i = 1:nEntries
    line = fgetl(fid);
    parts = strsplit(line);
    data{i, 1} = parts{1}; % Atom name
    data{i, 2} = str2double(parts{2}); % x-coordinate
    data{i, 3} = str2double(parts{3}); % y-coordinate
    data{i, 4} = str2double(parts{4}); % z-coordinate
end

fclose(fid);

% Convert cell array to table
pos = cell2table(data, 'VariableNames', {'atom', 'x', 'y', 'z'});
pos.atom = categorical(pos.atom);
end
