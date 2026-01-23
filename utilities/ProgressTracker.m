classdef ProgressTracker < handle
    % PROGRESSTRACKER Unified progress tracking for long-running operations
    %
    % Provides both command-line and graphical progress feedback with
    % automatic time estimation.
    %
    % USAGE:
    %   % Simple iteration tracking
    %   prog = ProgressTracker(nIterations, 'Processing atoms');
    %   for i = 1:nIterations
    %       % ... do work ...
    %       prog.update(i);
    %   end
    %   prog.finish();
    %
    %   % With sub-operations
    %   prog = ProgressTracker(nFiles, 'Processing files', 'showGUI', true);
    %   for i = 1:nFiles
    %       prog.update(i, sprintf('File %d: %s', i, filenames{i}));
    %   end
    %   prog.finish();
    %
    %   % Manual percentage updates
    %   prog = ProgressTracker(100, 'Computing', 'mode', 'percentage');
    %   prog.update(25);  % 25% complete
    %
    % OPTIONS:
    %   'showGUI'     - Show graphical waitbar (default: false)
    %   'mode'        - 'iteration' or 'percentage' (default: 'iteration')
    %   'minInterval' - Minimum seconds between display updates (default: 0.5)
    %   'showETA'     - Show estimated time remaining (default: true)
    %
    % (c) by Prof. Peter Felfer Group @FAU Erlangen-Nurnberg

    properties (Access = private)
        totalSteps
        currentStep
        taskName
        startTime
        lastUpdateTime
        minInterval
        showGUI
        showETA
        waitbarHandle
        mode
        isFinished
        lastMessage
    end

    methods
        function obj = ProgressTracker(totalSteps, taskName, options)
            arguments
                totalSteps (1,1) double {mustBePositive}
                taskName (1,:) char = 'Processing'
                options.showGUI (1,1) logical = false
                options.mode (1,:) char {mustBeMember(options.mode, {'iteration', 'percentage'})} = 'iteration'
                options.minInterval (1,1) double = 0.5
                options.showETA (1,1) logical = true
            end

            obj.totalSteps = totalSteps;
            obj.taskName = taskName;
            obj.currentStep = 0;
            obj.startTime = tic;
            obj.lastUpdateTime = 0;
            obj.minInterval = options.minInterval;
            obj.showGUI = options.showGUI;
            obj.showETA = options.showETA;
            obj.mode = options.mode;
            obj.isFinished = false;
            obj.lastMessage = '';

            if obj.showGUI
                obj.waitbarHandle = waitbar(0, obj.taskName, ...
                    'Name', 'Progress', ...
                    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
                setappdata(obj.waitbarHandle, 'canceling', 0);
            else
                % Print initial message
                fprintf('%s: 0%%', obj.taskName);
            end
        end

        function cancelled = update(obj, step, message)
            % Update progress
            % Returns true if user cancelled (GUI mode only)

            arguments
                obj
                step (1,1) double
                message (1,:) char = ''
            end

            cancelled = false;
            obj.currentStep = step;

            % Check for cancellation in GUI mode
            if obj.showGUI && isvalid(obj.waitbarHandle)
                if getappdata(obj.waitbarHandle, 'canceling')
                    cancelled = true;
                    obj.finish('Cancelled');
                    return;
                end
            end

            % Throttle updates
            elapsed = toc(obj.startTime);
            if elapsed - obj.lastUpdateTime < obj.minInterval && step < obj.totalSteps
                return;
            end
            obj.lastUpdateTime = elapsed;

            % Calculate progress
            if strcmp(obj.mode, 'percentage')
                progress = step / 100;
                pctComplete = step;
            else
                progress = step / obj.totalSteps;
                pctComplete = progress * 100;
            end

            % Estimate time remaining
            etaStr = '';
            if obj.showETA && progress > 0.01
                totalTime = elapsed / progress;
                remainingTime = totalTime - elapsed;
                etaStr = obj.formatTime(remainingTime);
            end

            % Build status message
            if isempty(message)
                if strcmp(obj.mode, 'percentage')
                    statusMsg = sprintf('%s: %.0f%%', obj.taskName, pctComplete);
                else
                    statusMsg = sprintf('%s: %d/%d (%.0f%%)', ...
                        obj.taskName, step, obj.totalSteps, pctComplete);
                end
            else
                statusMsg = sprintf('%s: %s', obj.taskName, message);
            end

            if ~isempty(etaStr)
                statusMsg = sprintf('%s - ETA: %s', statusMsg, etaStr);
            end

            % Update display
            if obj.showGUI && isvalid(obj.waitbarHandle)
                waitbar(progress, obj.waitbarHandle, statusMsg);
            else
                % Clear previous line and print new status
                fprintf(repmat('\b', 1, length(obj.lastMessage)));
                fprintf('%s', statusMsg);
                obj.lastMessage = statusMsg;
            end
        end

        function finish(obj, message)
            % Mark progress as complete

            arguments
                obj
                message (1,:) char = 'Complete'
            end

            if obj.isFinished
                return;
            end
            obj.isFinished = true;

            elapsed = toc(obj.startTime);
            elapsedStr = obj.formatTime(elapsed);

            if obj.showGUI
                if isvalid(obj.waitbarHandle)
                    delete(obj.waitbarHandle);
                end
            else
                fprintf(repmat('\b', 1, length(obj.lastMessage)));
                fprintf('%s: %s (%s)\n', obj.taskName, message, elapsedStr);
            end
        end

        function delete(obj)
            % Destructor - ensure waitbar is closed
            if obj.showGUI && ~isempty(obj.waitbarHandle) && isvalid(obj.waitbarHandle)
                delete(obj.waitbarHandle);
            end
            if ~obj.isFinished
                fprintf('\n');  % Ensure newline
            end
        end
    end

    methods (Static)
        function str = formatTime(seconds)
            % Format seconds into human-readable string
            if seconds < 60
                str = sprintf('%.0fs', seconds);
            elseif seconds < 3600
                mins = floor(seconds / 60);
                secs = mod(seconds, 60);
                str = sprintf('%dm %02.0fs', mins, secs);
            else
                hours = floor(seconds / 3600);
                mins = floor(mod(seconds, 3600) / 60);
                str = sprintf('%dh %02dm', hours, mins);
            end
        end

        function demo()
            % Demonstrate ProgressTracker usage
            fprintf('Demo 1: Command-line progress\n');
            n = 50;
            prog = ProgressTracker(n, 'Processing items');
            for i = 1:n
                pause(0.05);  % Simulate work
                prog.update(i);
            end
            prog.finish();

            fprintf('\nDemo 2: With custom messages\n');
            files = {'data1.pos', 'data2.pos', 'data3.pos', 'data4.pos', 'data5.pos'};
            prog = ProgressTracker(length(files), 'Loading files');
            for i = 1:length(files)
                pause(0.3);  % Simulate work
                prog.update(i, files{i});
            end
            prog.finish();

            fprintf('\nDemo complete!\n');
        end
    end
end
