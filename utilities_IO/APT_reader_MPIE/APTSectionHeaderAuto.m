% paraprobe-transcoder
% script to reader Cameca/AMETEK *.APT files
% Markus K\"uhbach, m.kuehbach@mpie.de, 2020/03/19

classdef APTSectionHeaderAuto
    properties
        %the variables and members (C/C++) the class object holds
        healthy;
        %aptfn;
        cSignature;
        iHeaderSize;
        iHeaderVersion;
        wcSectionType;
        iSectionVersion;
        eRelationshipType;
        eRecordType;
        eRecordDataType;
        iDataTypeSize;
        iRecordSize;
        wcDataUnit;
        llRecordCount;
        llByteCount;
        iElements;
    end
    methods
        %the methods fucntions (C/C++) the class object can execute
        function obj = APTSectionHeaderAuto(fid)
            %constructor for the class object/instance, will throw error
            %when section is not of desired sectionname
            if nargin == 1
                obj.healthy = true;
                %obj.aptfn = fn;
                obj.cSignature = "";
                obj.iHeaderSize = 0;
                obj.iHeaderVersion = 0;
                obj.wcSectionType = 'Failure';
                obj.iSectionVersion = 0;
                obj.eRelationshipType = 0;
                obj.eRecordType = 0;
                obj.eRecordDataType = 0;
                obj.iDataTypeSize = 0;
                obj.iRecordSize = 0;
                obj.wcDataUnit = 0;
                obj.llRecordCount = 0;
                obj.llByteCount = 0;
                
                %parse out the header
                %offset = 540;
                %sectiontype = 'Mass';
                %fn = 'Z:/GITHUB/MPIE_APTFIM_TOOLBOX/paraprobe/code/84d83abc-537d-4795-854e-b65b53d5c0c9.apt';
                
                %see additional exemplary comments in APTFileHeader
                %fid = fopen(obj.aptfn, 'rb');
                %start reading file at byte offset
                %offs = offset;
                %fseek(fid, offs, 'bof'); %start reading at byte position offset from the beginning of the file
                obj.cSignature = strcat(char(fread(fid,[1,4],'uint8'))); %no additional skip on top of implicit fp advancement
                if strcmp(obj.cSignature,'SEC') == false
                    obj.healthy = false;
                    disp('Section is not a valid SEC!');
                    obj.wcSectionType = 'Failure'; return;
                end
                obj.iHeaderSize = fread(fid,[1,1],'int32');
                obj.iHeaderVersion = fread(fid,[1,1],'int32');
                if obj.iHeaderVersion ~= 2
                    obj.healthy = false;
                    disp('Section header has a different version!');
                    obj.wcSectionType = 'Failure'; return;
                end
                %##MK::implement version testing
                %##MK::UTF16 not yet implemented, instead currently we read the uint16 raw and interpret into UTF8 !!
                tmp = char(fread(fid,[1,32],'uint16'));
                if sum(tmp > 255) > 0
                    obj.healthy = false;
                    disp('UTF16 conversion of wcSectionType did not work!');
                    obj.wcSectionType = 'Failure'; return;
                end       
                obj.wcSectionType = strcat(tmp);
                obj.iSectionVersion = fread(fid,[1, 1],'int32');
                obj.eRelationshipType = fread(fid,[1, 1],'uint32');
                if obj.eRelationshipType ~= 1
                    obj.healthy = false;
                    disp('eRelationshipType ~= ONE_TO_ONE but that is the only currently support!');
                    obj.wcSectionType = 'Failure'; return;
                end                    
                obj.eRecordType = fread(fid,[1, 1],'uint32');
                if obj.eRecordType ~= 1
                    obj.healthy = false;
                    disp('eRecordType ~= FIXED_SIZE but that is the only currently support!');
                    obj.wcSectionType = 'Failure'; return;
                end                
                obj.eRecordDataType = fread(fid,[1, 1],'uint32');
                obj.iDataTypeSize = fread(fid,[1, 1],'int32');
                obj.iRecordSize = fread(fid,[1,1], 'int32');
                obj.wcDataUnit = strcat(char(fread(fid,[1,16],'uint16')));
                obj.llRecordCount = fread(fid,[1,1], 'uint64');
                obj.llByteCount = fread(fid,[1,1], 'uint64');
                obj.iElements = obj.iRecordSize / (obj.iDataTypeSize/8);
                %fclose(fid);
                disp(['Reading *.APT section __' obj.wcSectionType '__ successful']);                
            end
        end      
        function print(obj)
            disp(['cSignature ' obj.cSignature]);
            disp(['iHeaderSize ' num2str(obj.iHeaderSize)]);
            disp(['iHeaderVersion ' num2str(obj.iHeaderVersion)]);
            disp(['wcSectionType ' obj.wcSectionType]);
            disp(['iSectionVersion ' num2str(obj.iSectionVersion)]);
            disp(['eRelationshipType ' num2str(obj.eRelationshipType)]);
            disp(['eRecordType ' num2str(obj.eRecordType)]);
            disp(['eRecordDataType ' num2str(obj.eRecordDataType)]);
            disp(['iDataTypeSize ' num2str(obj.iDataTypeSize)]);
            disp(['iRecordSize ' num2str(obj.iRecordSize)]);
            disp(['wcDataUnit ' obj.wcDataUnit]);
            disp(['llRecordCount ' num2str(obj.llRecordCount)]);
            disp(['llByteCount ' num2str(obj.llByteCount)]);
            disp(['iElements ' num2str(obj.iElements)]);
        end
        %function offs = get_offset(obj)
        %    %##MK::possibly different for position
        %    offs = 148;
        %end
    end
end

