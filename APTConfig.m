classdef APTConfig < handle
    % APTCONFIG Global configuration management for Atom Probe Toolbox
    %
    % Singleton class providing centralized configuration management.
    % User settings persist across sessions via a JSON config file.
    %
    % USAGE:
    %   cfg = APTConfig.getInstance();
    %   cfg.reconstruction.detectionEfficiency = 0.52;
    %   cfg.save();
    %
    %   % Or use static methods:
    %   APTConfig.set('reconstruction.detectionEfficiency', 0.52);
    %   eta = APTConfig.get('reconstruction.detectionEfficiency');
    %
    % CONFIGURATION CATEGORIES:
    %   reconstruction - Reconstruction parameters
    %   analysis       - Analysis defaults
    %   visualization  - Plotting and display settings
    %   io             - Input/output settings
    %   performance    - Parallel computing and memory settings
    %
    % (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

    properties (Access = private)
        configFilePath
    end

    properties
        % Reconstruction parameters
        reconstruction = struct(...
            'detectionEfficiency', 0.37, ...  % Default for LEAP 4000X HR
            'flightPathLength', 110, ...      % mm
            't0', 55.7, ...                   % ns, field-free region time
            'kFactor', 3.3, ...               % Image compression factor
            'ICF', 1.65, ...                  % Image compression factor
            'fieldFactor', 0.18, ...          % V/nm
            'atomicVolume', 0.012, ...        % nm^3 (default for Fe)
            'method', 'Geiser' ...            % Reconstruction algorithm
        )

        % Analysis parameters
        analysis = struct(...
            'defaultBinWidth', 0.1, ...       % Da for mass spectrum
            'concentrationKernelSize', 100, ... % atoms for local concentration
            'clusterMinPts', 10, ...          % Minimum atoms for cluster
            'clusterEpsilon', 0.5, ...        % nm, DBSCAN neighborhood
            'voronoiMaxVolume', 2.0, ...      % nm^3, max Voronoi cell volume
            'proxigramBinWidth', 0.2, ...     % nm
            'spatialResolution', 0.3 ...      % nm, estimated lateral resolution
        )

        % Visualization settings
        visualization = struct(...
            'defaultColormap', 'parula', ...
            'atomMarkerSize', 2, ...
            'atomMarkerAlpha', 0.5, ...
            'meshFaceAlpha', 0.5, ...
            'meshEdgeAlpha', 0.1, ...
            'figurePosition', [100 100 800 600], ...
            'exportDPI', 300, ...
            'defaultFontSize', 12, ...
            'axesEqual', true, ...
            'showAxesBox', true ...
        )

        % I/O settings
        io = struct(...
            'defaultDataPath', '', ...
            'defaultExportPath', '', ...
            'hdf5Compression', 'deflate', ...
            'hdf5CompressionLevel', 6, ...
            'autoBackup', true, ...
            'backupInterval', 300 ...         % seconds
        )

        % Performance settings
        performance = struct(...
            'useParallel', true, ...          % Use parallel computing if available
            'maxWorkers', 0, ...              % 0 = use all available
            'useGPU', true, ...               % Use GPU if available
            'gpuFallbackToCPU', true, ...     % Fall back to CPU if GPU fails
            'memoryLimit', 0.8, ...           % Fraction of available memory to use
            'chunkSize', 1000000, ...         % Atoms per chunk for large datasets
            'showProgress', true, ...         % Show progress bars
            'verboseOutput', true ...         % Print status messages
        )

        % Instrument presets
        instrumentPresets = struct(...
            'LEAP4000XHR', struct(...
                'detectionEfficiency', 0.37, ...
                'flightPathLength', 110, ...
                't0', 55.7), ...
            'LEAP5000XR', struct(...
                'detectionEfficiency', 0.52, ...
                'flightPathLength', 120, ...
                't0', 50.0), ...
            'LEAP5000XS', struct(...
                'detectionEfficiency', 0.80, ...
                'flightPathLength', 120, ...
                't0', 50.0), ...
            'EIKOS', struct(...
                'detectionEfficiency', 0.50, ...
                'flightPathLength', 100, ...
                't0', 45.0) ...
        )

        % Material presets (atomic volumes in nm^3)
        materialPresets = struct(...
            'Fe', 0.01178, ...
            'Al', 0.01660, ...
            'Cu', 0.01182, ...
            'Ni', 0.01094, ...
            'Ti', 0.01766, ...
            'W',  0.01583, ...
            'Si', 0.02000, ...
            'Mg', 0.02305 ...
        )
    end

    methods (Access = private)
        function obj = APTConfig()
            % Private constructor for singleton
            toolboxPath = fileparts(mfilename('fullpath'));
            obj.configFilePath = fullfile(toolboxPath, '.aptconfig.json');
            obj.load();
        end
    end

    methods (Static)
        function obj = getInstance()
            % Get the singleton instance
            persistent instance
            if isempty(instance) || ~isvalid(instance)
                instance = APTConfig();
            end
            obj = instance;
        end

        function value = get(key)
            % Static method to get a configuration value
            % Example: APTConfig.get('reconstruction.detectionEfficiency')
            cfg = APTConfig.getInstance();
            value = cfg.getValue(key);
        end

        function set(key, value)
            % Static method to set a configuration value
            % Example: APTConfig.set('reconstruction.detectionEfficiency', 0.52)
            cfg = APTConfig.getInstance();
            cfg.setValue(key, value);
        end

        function reset()
            % Reset configuration to defaults
            cfg = APTConfig.getInstance();
            cfg.resetToDefaults();
        end

        function show()
            % Display current configuration
            cfg = APTConfig.getInstance();
            cfg.display();
        end
    end

    methods
        function value = getValue(obj, key)
            % Get a configuration value using dot notation
            parts = strsplit(key, '.');
            value = obj;
            for i = 1:length(parts)
                if isstruct(value)
                    value = value.(parts{i});
                else
                    value = value.(parts{i});
                end
            end
        end

        function setValue(obj, key, value)
            % Set a configuration value using dot notation
            parts = strsplit(key, '.');
            if length(parts) == 1
                obj.(parts{1}) = value;
            elseif length(parts) == 2
                obj.(parts{1}).(parts{2}) = value;
            else
                error('APTConfig:nestedTooDeep', 'Configuration nesting too deep');
            end
        end

        function save(obj)
            % Save configuration to file
            configData = struct();
            configData.reconstruction = obj.reconstruction;
            configData.analysis = obj.analysis;
            configData.visualization = obj.visualization;
            configData.io = obj.io;
            configData.performance = obj.performance;

            try
                jsonStr = jsonencode(configData);
                fid = fopen(obj.configFilePath, 'w');
                if fid == -1
                    warning('APTConfig:saveError', 'Could not save configuration file');
                    return;
                end
                fprintf(fid, '%s', jsonStr);
                fclose(fid);
            catch ME
                warning('APTConfig:saveError', 'Error saving configuration: %s', ME.message);
            end
        end

        function load(obj)
            % Load configuration from file
            if ~isfile(obj.configFilePath)
                return;  % Use defaults
            end

            try
                fid = fopen(obj.configFilePath, 'r');
                if fid == -1
                    return;
                end
                jsonStr = fread(fid, '*char')';
                fclose(fid);

                configData = jsondecode(jsonStr);

                % Merge with existing (to preserve new fields in code)
                obj.mergeStruct('reconstruction', configData);
                obj.mergeStruct('analysis', configData);
                obj.mergeStruct('visualization', configData);
                obj.mergeStruct('io', configData);
                obj.mergeStruct('performance', configData);
            catch ME
                warning('APTConfig:loadError', 'Error loading configuration: %s', ME.message);
            end
        end

        function mergeStruct(obj, fieldName, configData)
            % Merge loaded config with existing defaults
            if isfield(configData, fieldName)
                loadedFields = fieldnames(configData.(fieldName));
                for i = 1:length(loadedFields)
                    fn = loadedFields{i};
                    if isfield(obj.(fieldName), fn)
                        obj.(fieldName).(fn) = configData.(fieldName).(fn);
                    end
                end
            end
        end

        function resetToDefaults(obj)
            % Reset all settings to defaults
            defaultConfig = APTConfig();
            obj.reconstruction = defaultConfig.reconstruction;
            obj.analysis = defaultConfig.analysis;
            obj.visualization = defaultConfig.visualization;
            obj.io = defaultConfig.io;
            obj.performance = defaultConfig.performance;

            % Delete config file
            if isfile(obj.configFilePath)
                delete(obj.configFilePath);
            end
        end

        function applyInstrumentPreset(obj, instrumentName)
            % Apply instrument-specific settings
            % Example: cfg.applyInstrumentPreset('LEAP5000XR')

            if ~isfield(obj.instrumentPresets, instrumentName)
                validNames = fieldnames(obj.instrumentPresets);
                error('APTConfig:unknownInstrument', ...
                    'Unknown instrument: %s. Valid options: %s', ...
                    instrumentName, strjoin(validNames, ', '));
            end

            preset = obj.instrumentPresets.(instrumentName);
            fields = fieldnames(preset);
            for i = 1:length(fields)
                obj.reconstruction.(fields{i}) = preset.(fields{i});
            end

            fprintf('Applied instrument preset: %s\n', instrumentName);
        end

        function applyMaterialPreset(obj, materialName)
            % Apply material-specific settings (atomic volume)
            % Example: cfg.applyMaterialPreset('Al')

            if ~isfield(obj.materialPresets, materialName)
                validNames = fieldnames(obj.materialPresets);
                error('APTConfig:unknownMaterial', ...
                    'Unknown material: %s. Valid options: %s', ...
                    materialName, strjoin(validNames, ', '));
            end

            obj.reconstruction.atomicVolume = obj.materialPresets.(materialName);
            fprintf('Applied material preset: %s (atomic volume = %.5f nm^3)\n', ...
                materialName, obj.reconstruction.atomicVolume);
        end

        function display(obj)
            % Display current configuration
            fprintf('\n');
            fprintf('==========================================\n');
            fprintf('  Atom Probe Toolbox Configuration\n');
            fprintf('==========================================\n\n');

            categories = {'reconstruction', 'analysis', 'visualization', 'io', 'performance'};

            for c = 1:length(categories)
                cat = categories{c};
                fprintf('%s:\n', upper(cat));
                fprintf('------------------------------------------\n');

                fields = fieldnames(obj.(cat));
                for f = 1:length(fields)
                    fn = fields{f};
                    val = obj.(cat).(fn);
                    if isnumeric(val)
                        fprintf('  %-25s: %g\n', fn, val);
                    elseif islogical(val)
                        if val
                            fprintf('  %-25s: true\n', fn);
                        else
                            fprintf('  %-25s: false\n', fn);
                        end
                    elseif ischar(val)
                        fprintf('  %-25s: %s\n', fn, val);
                    end
                end
                fprintf('\n');
            end

            fprintf('Config file: %s\n', obj.configFilePath);
            fprintf('==========================================\n\n');
        end

        function tf = useGPU(obj)
            % Check if GPU should be used (available and enabled)
            tf = false;
            if ~obj.performance.useGPU
                return;
            end
            try
                if exist('gpuDeviceCount', 'file') && gpuDeviceCount > 0
                    tf = true;
                end
            catch
                tf = false;
            end
        end

        function tf = useParallel(obj)
            % Check if parallel computing should be used
            tf = false;
            if ~obj.performance.useParallel
                return;
            end
            if ~isempty(ver('parallel'))
                tf = true;
            end
        end

        function n = getNumWorkers(obj)
            % Get number of parallel workers to use
            n = 0;
            if ~obj.useParallel()
                return;
            end

            if obj.performance.maxWorkers > 0
                n = obj.performance.maxWorkers;
            else
                % Use all available
                c = parcluster('local');
                n = c.NumWorkers;
            end
        end
    end
end
