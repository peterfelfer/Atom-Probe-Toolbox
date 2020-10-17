function voxelsToRAW(volume, gv, fileName)
%NEEDS DOCUMENTATION
%saves a volume value matrix to an 8bit raw file for volumetric analysis
%and display in external programs such as Blender


if ~exist('fileName','var')
    [file, path] = uiputfile('*.RAW','Select RAW file');
    fileName = [path file];
    
else
    fileName = [fileName '.RAW'];
    
end




%% Info about data in file to be put in file name
%number of voxels
sz = size(volume);

%extent of volume
if exist('gv','var')
    for i=1:3
        %inc(i) = gv{i}(2) - gv{i}(1);
        d(i) = gv{i}(end) - gv{i}(1);% + inc(i);
    end
    
end

%minimum and maximum values
mn = min(volume(:));
mx = max(volume(:));

%embedding data in filename
[pathstr,name,ext] = fileparts(fileName);

name = [name num2str(sz(1)) 'x' num2str(sz(2)) 'x' num2str(sz(3))];
name = [name '_' num2str(d(1)) 'nm' num2str(d(2)) 'nm' num2str(d(3)) 'nm'];
name = [name '_' num2str(mn,3) '-' num2str(mx,3)];

fileName = [path name ext];









%% data preparation
%set NaNs to 0
volume(isnan(volume)) = 0;

%clamping of values
volume = (volume - mn) / (mx - mn) * 255;

%conversion to 8bit integers
volume = uint8(volume);




%% writing of file
% opens the file
fid = fopen(fileName,'W','b');

disp(['File opened: ' name ext]);

% reshape the matrix to get one row for RAW file
mat=reshape(volume,numel(volume),1);

% write the matrix as floats
fwrite(fid,mat,'uint8');
fclose(fid);
disp('RAW file saved');


end