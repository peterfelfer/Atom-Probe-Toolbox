% paraprobe-transcoder
% script to reader Cameca/AMETEK *.APT files
% Markus K\"uhbach, m.kuehbach@mpie.de, 2020/03/19

classdef PARAPROBE_Transcoder2
    properties
        %the variables and members (C/C++) the class object holds
        % headers, sections, and other metadata
        healthy;
        aptfn;
        aptbranches;
        header;
        idtfyd_sections;
        tipbox;
        % heavy data
        tof;
        pulse;
        freq;
        tElapsed;
        erate;
        tstage;
        TargetErate;
        TargetFlux;
        pulseDelta;
        Pres;
        VAnodeMon;
        Temp;
        AmbTemp;
        FractureGuard;
        Vref;
        Noise;
        Uniformity;
        xstage;
        ystage;
        zstage;
        z;
        tofc;
        Mass;
        tofb;
        xs;
        ys;
        zs;
        rTip;
        zApex;
        zSphereCorr;
        XDet_mm;
        YDet_mm;
        Multiplicity;
        Vap;
        DetectorCoordinates;
        Position;
        % we create row matrices here because they are more performant than column matrices!
        known_sections;
    end
    methods
        %the methods fucntions (C/C++) the class object can execute
        function obj = PARAPROBE_Transcoder2(fn)
            %constructor for the class object/instance
            if nargin == 1
                obj.healthy = true;
                obj.aptfn = fn;
                %*.apt format as of 2020/03/20
                %FileHead, SectionHead('Mass'), Mass data,
                %SectionHead('Position'), xlo,xhi,ylo,yhi,zlo,zhi, pos data
                % define all *.APT branches possible 
                % as of 2020/03/20 and APSuite6/IVAS4
                obj.aptbranches = APTBranchesDict;
                
                %open the file for read
                %fn = 'Z:/GITHUB/MPIE_APTFIM_TOOLBOX/paraprobe/code/84d83abc-537d-4795-854e-b65b53d5c0c9.apt';
                fid = fopen(obj.aptfn, 'rb'); %open file aptfn read-only binary
                %convert to character array and then fuse the characters
                %like in C, file pointer is implicitly advanced with each fread
                %evidently fread advances the file pointer implicitly C-style
                %read four uint8 from beginning of the file (0), little-endian, 
                              
                % read the header of the file
                obj.header = APTFileHeader2(fid);
                obj.header.print();
                if obj.header.healthy == false
                    return;
                end
                
                %get file read offset skip = obj.header.get_offset();
                                
                %auto-detect all sections in the file search basically for
                %all known branches and take the existing ones
                for i = 1:1:length(obj.aptbranches.dict.keyword)
                    i
                    obj.idtfyd_sections{i} = APTSectionHeaderAuto(fid); %keep on parsing
                    if obj.idtfyd_sections{i}.healthy == false
                        disp(['Section ' num2str(i) ' failed!']);
                        return;
                    end
                    % ###check if not a duplicate
                    %obj.idtfyd_sections{i}.print();
                    
                    % load specific data based on what the section encodes
                    ni = obj.header.llIonCount; %number of ions assume FIXED_SIZE, ONE_TO_ONE mapping of the data 
                    nj = obj.idtfyd_sections{i}.iElements; %number of members in value tuple per ion
                    disp(['ni ' num2str(ni) ' nj ' num2str(nj)]);
                    obj.idtfyd_sections{i}
                    switch obj.idtfyd_sections{i}.wcSectionType
                        case 'Failure'
                            disp(['Autodetect the ' num2str(i) ' APT header failed!']); 
                            break;
                        case 'tof'
                            disp('tof');
                            % ##MK::it seems that the fread always casts to
                            % double anyway in an *.APT we have single
                            % precision floats will be upcasted to double
                            % and downcasted to single again, possible
                            % error in the order of 1.0e-5 > eps > 1.0e-6
                            obj.tof = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]); %in-place deinterleaving
                            %obj.tof = fread(fid,[1, n],'float32');
                            if ~isequal( size(obj.tof), [nj, ni])
                                disp('Reading tof failed!'); break;
                            end
                        case 'pulse'
                            disp('pulse');
                            obj.pulse = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.pulse), [nj, ni])
                                disp('Reading pulse failed!'); break;
                            end
                        case 'freq'
                            disp('freq');
                            obj.freq = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.freq), [nj, ni])
                                disp('Reading freq failed!'); break;
                            end
                        case 'tElapsed'
                            disp('tElapsed');
                            obj.tElapsed = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.tElapsed), [nj, ni])
                                disp('Reading tElapsed failed!'); break;
                            end
                        case 'erate'
                            disp('erate');
                            obj.erate = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.erate), [nj, ni])
                                disp('Reading erate failed!'); break;
                            end
                        case 'tstage'
                            disp('tstage');
                            obj.tstage = reshape(uint16(fread(fid,[1, ni*nj],'uint16')),[nj, ni]);
                            if ~isequal( size(obj.tstage), [nj, ni])
                                disp('Reading tstage failed!'); break;
                            end
                        case 'TargetErate'
                            disp('TargetErate');
                            obj.TargetErate = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.TargetErate), [nj, ni])
                                disp('Reading TargetErate failed!'); break;
                            end
                        case 'TargetFlux'
                            disp('TargetFlux');
                            obj.TargetFlux = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.TargetFlux), [nj, ni])
                                disp('Reading TargetFlux failed!'); break;
                            end
                        case 'pulseDelta'
                            disp('pulseDelta');
                            obj.pulseDelta = reshape(int16(fread(fid,[1, ni*nj],'int16')),[nj, ni]);
                            if ~isequal( size(obj.pulseDelta), [nj, ni])
                                disp('Reading pulseDelta failed!'); break;
                            end
                        case 'Pres'
                            disp('Pres');
                            obj.Pres = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Pres), [nj, ni])
                                disp('Reading Pres failed!'); break;
                            end
                        case 'VAnodeMon'
                            disp('VAnodeMon');
                            obj.VAnodeMon = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.VAnodeMon), [nj, ni])
                                disp('Reading VAnodeMon failed!'); break;
                            end
                        case 'Temp'
                            disp('Temp');
                            obj.Temp = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Temp), [nj, ni])
                                disp('Reading Temp failed!'); break;
                            end
                        case 'AmbTemp'
                            disp('AmbTemp');
                            obj.AmbTemp = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.AmbTemp), [nj, ni])
                                disp('Reading AmbTemp failed!'); break;
                            end
                        case 'FractureGuard'
                            disp('FractureGuard');
                            obj.FractureGuard = reshape(uint16(fread(fid,[1, ni*nj],'uint16')),[nj, ni]);
                            if ~isequal( size(obj.FractureGuard), [nj, ni])
                                disp('Reading FractureGuard failed!'); break;
                            end
                        case 'Voltage'
                            disp('Vref');
                            obj.Vref = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Vref), [nj, ni])
                                disp('Reading Vref failed!'); break;
                            end
                        case 'Noise'
                            disp('Noise');
                            obj.Noise = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Noise), [nj, ni])
                                disp('Reading Noise failed!'); break;
                            end
                        case 'Uniformity'
                            disp('Uniformity');
                            obj.Uniformity = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Uniformity), [nj, ni])
                                disp('Reading Uniformity failed!'); break;
                            end
                        case 'xstage'
                            disp('xstage');
                            obj.xstage = reshape(int32(fread(fid,[1, ni*nj],'int32')),[nj, ni]);
                            if ~isequal( size(obj.xstage), [nj, ni])
                                disp('Reading xstage failed!'); break;
                            end
                        case 'ystage'
                            disp('ystage');
                            obj.ystage = reshape(int32(fread(fid,[1, ni*nj],'int32')),[nj, ni]);
                            if ~isequal( size(obj.ystage), [nj, ni])
                                disp('Reading ystage failed!'); break;
                            end
                        case 'zstage'
                            disp('zstage');
                            obj.zstage = reshape(int32(fread(fid,[1, ni*nj],'int32')),[nj, ni]);
                            if ~isequal( size(obj.zstage), [nj, ni])
                                disp('Reading zstage failed!'); break;
                            end
                        case 'z'
                            disp('z');
                            obj.z = reshape(uint64(fread(fid,[1, ni*nj],'uint64')),[nj, ni]);
                            if ~isequal( size(obj.z), [nj, ni])
                                disp('Reading z failed!'); break;
                            end
                        case 'tofc'
                            disp('tofc');
                            obj.tofc = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.tofc), [nj, ni])
                                disp('Reading tofc failed!'); break;
                            end
                        case 'Mass'
                            disp('Mass');
                            obj.Mass = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Mass), [nj, ni])
                                disp('Reading Mass failed!'); break;
                            end
                        case 'tofb'
                            disp('tofb');
                            obj.tofb = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.tofb), [nj, ni])
                                disp('Reading tofb failed!'); break;
                            end
                        case 'xs'
                            disp('xs');
                            obj.xs = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.xs), [nj, ni])
                                disp('Reading xs failed!'); break;
                            end 
                        case 'ys'
                            disp('ys');
                            obj.ys = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.ys), [nj, ni])
                                disp('Reading ys failed!'); break;
                            end 
                        case 'zs'
                            disp('zs');
                            obj.zs = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.zs), [nj, ni])
                                disp('Reading zs failed!'); break;
                            end
                        case 'rTip'
                            disp('rTip');
                            obj.rTip = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.rTip), [nj, ni])
                                disp('Reading rTip failed!'); break;
                            end
                        case 'zApex'
                            disp('zApex');
                            obj.zApex = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.zApex), [nj, ni])
                                disp('Reading zApex failed!'); break;
                            end
                        case 'zSphereCorr'
                            disp('zSphereCorr');
                            obj.zSphereCorr = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.zSphereCorr), [nj, ni])
                                disp('Reading zSphereCorr failed!'); break;
                            end
                        case 'XDet_mm'
                            disp('XDet_mm');
                            obj.XDet_mm = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.XDet_mm), [nj, ni])
                                disp('Reading XDet_mm failed!'); break;
                            end
                        case 'YDet_mm'
                            disp('YDet_mm');
                            obj.YDet_mm = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.YDet_mm), [nj, ni])
                                disp('Reading YDet_mm failed!'); break;
                            end
                        case 'Multiplicity'
                            disp('Multiplicity');
                            obj.Multiplicity = reshape(int32(fread(fid,[1, ni*nj],'int32')),[nj, ni]);
                            if ~isequal( size(obj.Multiplicity), [nj, ni])
                                disp('Reading Multiplicity failed!'); break;
                            end
                        case 'Vap'
                            disp('Vap');
                            obj.Vap = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Vap), [nj, ni])
                                disp('Reading Vap failed!'); break;
                            end
                        case 'Detector Coordinates'
                            disp('Detector Coordinates');
                            obj.DetectorCoordinates = reshape(int32(fread(fid,[1, ni*nj],'int32')),[nj, ni]);
                            if ~isequal( size(obj.DetectorCoordinates), [nj, ni])
                                disp('Reading DetectorCoordinates failed!'); break;
                            end
                        case 'Position'
                            disp('Position');
                            % first read the preceeding header with bounds
                            obj.tipbox = reshape(single(fread(fid,[1, 6],'float32')),[2, 3]); % ##check
                            if ~isequal( size(obj.tipbox), [2, 3])
                                disp('Reading Position box bounds failed!'); break;
                            end
                            obj.Position = reshape(single(fread(fid,[1, ni*nj],'float32')),[nj, ni]);
                            if ~isequal( size(obj.Position), [nj, ni])
                                disp('Reading Position failed!'); break;
                            end                     
                        otherwise
                            disp('Stop reading sections');
                            break; %break out of the for loop
                    end
                end 
                disp('Done, closing the *.APT file');
                fclose(fid);
                
%                 disp(['skip ' num2str(skip)]);
%                 obj.mass_sect = APTSectionHeader(fn,'Mass', skip);
%                 obj.mass_sect.print();
%                 if obj.healthy == false
%                     return;
%                 end
%                 
%                 skip = skip + obj.mass_sect.get_offset();
%                 disp(['skip ' num2str(skip)]);
%                 fid = fopen(obj.aptfn, 'rb'); %##MK::optimize later for smaller number of file management
%                 fseek(fid, skip,'bof');
%                 obj.mq = fread(fid,[1, obj.header.llIonCount*1],'float32',0,'l');
%                 fclose(fid);
%                 
%                 skip = skip + obj.header.llIonCount*1*4;
%                 disp(['skip ' num2str(skip)]);
%                 obj.xyz_sect = APTSectionHeader(fn,'Position', skip);
%                 obj.xyz_sect.print();
%                 if obj.healthy == false;
%                     return;
%                 end
%                 % after the position section there are the recon bounds
%                 skip = skip + obj.xyz_sect.get_offset();
%                 disp(['skip ' num2str(skip)]);
%                 fid = fopen(obj.aptfn, 'rb');
%                 fseek(fid,skip,'bof');
%                 %readinfg with in-place deinterleacing
%                 %resulting in a [xmi,ymi,zmi; xmx, ymx, zmx] matrix ! i.e.
%                 %x,y,z rows instead of columns!
%                 bounds = reshape(fread(fid,[1,6],'float32',0,'l'),[2,3]);
%                 disp('bounds read');
%                 disp(bounds);
%                 
%                 %no skip now instead, move on an read position data               
%                 %skip = skip + 6*4;
%                 %disp(['skip ' num2str(skip)]);
%                 %load the position data with in-place deinterleacing
%                 % BE CAREFUL x,y,z are now as rows instead of columns because
%                 % Matlab inherits Fortran-style indexing, i.e.
%                 % [3,llIonCount] is faster than using [llIonCount,3] matrices ...
%                 obj.xyz = reshape(fread(fid,[1, obj.header.llIonCount*3],'float32', 0,'l'),[3, obj.header.llIonCount]);
%                 disp(['xyz read']);
%                 disp(['xyz shape']);
%                 disp(size(obj.xyz));                
%                 fclose(fid);
            end
        end
    end
end

