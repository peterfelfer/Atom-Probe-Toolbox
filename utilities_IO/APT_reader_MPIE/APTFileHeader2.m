% paraprobe-transcoder
% script to reader Cameca/AMETEK *.APT files
% Markus K\"uhbach, m.kuehbach@mpie.de, 2020/03/19

classdef APTFileHeader2
    properties
        %the variables and members (C/C++) the class object holds
        healthy;
        %aptfn;
        cSignature;
        iHeaderSize;
        iHeaderVersion;
        wcFilename;
        ftCreationTime;
        llIonCount;
    end
    methods
        %the methods fucntions (C/C++) the class object can execute
        function obj = APTFileHeader2(fid)
            %constructor for the class object/instance
            if nargin == 1
                obj.healthy = true;
                %obj.aptfn = fn;
                obj.cSignature = "";
                obj.iHeaderSize = 0;
                obj.iHeaderVersion = 0;
                obj.wcFilename = "";
                obj.ftCreationTime = 0;
                obj.llIonCount = 0;
                
                %parse out the header
                %fn = 'Z:/GITHUB/MPIE_APTFIM_TOOLBOX/paraprobe/code/84d83abc-537d-4795-854e-b65b53d5c0c9.apt';
                %##MK::explore memmapfile function for a less verbose way of mapping data
                %fid = fopen(obj.aptfn, 'rb'); %open fn file read-only binary
                %read four uint8 from beginning of the file (0), little-endian, 
                %convert to character array and then fuse the characters
                %offs = 0; %eventual skip in addition !! to already implicitly advanced file pointer
                %evidently fread advances the file pointer implicitly C-style
                obj.cSignature = strcat(char(fread(fid,[1,4],'uint8'))); %offs = offs + 4*1;
                if strcmp(obj.cSignature,'APT') == false
                    obj.healthy = false;
                    disp('File is not a valid APT file!'); 
                    return;
                end
                obj.iHeaderSize = fread(fid,[1,1],'int32'); %little endian by default
                obj.iHeaderVersion = fread(fid,[1,1],'int32');
                %##MK::implement version testing
                %##MK::UTF16 not yet implemented, instead currently we read the uint16 raw and interpret into UTF8 !!
                obj.wcFilename = strcat(char(fread(fid,[1,256],'uint16')));
                obj.ftCreationTime = fread(fid,[1,1],'uint64');
                obj.llIonCount = fread(fid,[1,1],'uint64');
                %fclose(fid);
                disp(['Reading *.APT file header was successful']);
                disp([num2str(obj.llIonCount) ' ions']);
            end
        end
        function print(obj)
            disp(['cSignature ' obj.cSignature]);
            disp(['iHeaderSize ' num2str(obj.iHeaderSize)]);
            disp(['iHeaderVersion ' num2str(obj.iHeaderVersion)]);
            disp(['wcFilename ' obj.wcFilename]);
            disp(['ftCreationTime ' num2str(obj.ftCreationTime)]);
            disp(['llIonCount ' num2str(obj.llIonCount)]);
        end
        %function offs = get_offset(obj)
        %    offs = 540;
        %end
    end
end

