% SIGLENT SSA3021X: Absorption Dip Detection

% Parameters
analyzerIP = '192.168.0.50';
port = 5025;
timeout = 120;

scanStartFreq  = 10e6;        % determined by step 1 rough scan (find comb profile)
scanStopFreq   = 20e6;        % determined by step 1 rough scan (find comb profile)
segmentSpan    = 400e3;       % 400 kHz segments
targetAvg      = 50;          % 500 averages
analogRBW      = 300;         % 300 Hz RBW
dipThresholdDB = 2.5;         % 2.5 dB threshold
expectedWidthHz = 4000;       % Expected physical dip width (4000 Hz)

try
    fprintf('STATUS: Checking TCP connection to %s...\n', analyzerIP);
    if ~exist('ssa', 'var')
        fprintf('STATUS: Establishing new connection...\n');
        ssa = tcpclient(analyzerIP, port, 'Timeout', timeout);
        flush(ssa); % Clean pipe on first boot
        writeline(ssa, '*IDN?');
        fprintf('STATUS: Connected -> %s\n', readline(ssa));
    else
        try
            % Remove any trapped data from a previous aborted run
            flush(ssa);
            
            writeline(ssa, '*IDN?');
            readline(ssa);
            fprintf('STATUS: Reusing existing TCP connection.\n');
        catch
            fprintf('STATUS: Old connection died. Rebuilding socket...\n');
            clear ssa;
            pause(2.0);
            ssa = tcpclient(analyzerIP, port, 'Timeout', timeout);
            flush(ssa);
            writeline(ssa, '*IDN?');
            fprintf('STATUS: Connected -> %s\n', readline(ssa));
        end
    end
    
    % Instrument setup
    fprintf('STATUS: Configuring global instrument parameters...\n');
    writeline(ssa, ':INSTrument SANalyzer');
    writeline(ssa, ':DISPlay:WINDow:TRACe:Y:SCALe:RLEVel -90');
    writeline(ssa, ':DISPlay:WINDow:TRACe:Y:SCALe:PDIVision 5');
    writeline(ssa, ':SENSe:BWIDth:RESolution ' + string(analogRBW));
    
    writeline(ssa, ':TRACe1:MODE AVERage');
    writeline(ssa, ':AVERage:TRACe1:STATe ON');
    writeline(ssa, ':AVERage:TRACe1:COUNt ' + string(targetAvg));
    writeline(ssa, ':DETector:TRACe1 AVERage');
    
    % Disk streaming setup
    filename = createUniqueFile('rf_absorption_scan', '.txt');
    fileID = fopen(filename, 'a');
    if fileID == -1, error('Failed to open file for writing.'); end

    dipFilename = createUniqueFile('dip_locations', '.txt');
    dipFileID = fopen(dipFilename, 'a');
    if dipFileID == -1, error('Failed to open dip log file.'); end
    fprintf(dipFileID, 'RF_Frequency_MHz, Dip_Depth_dB\n');

    totalSegments = round((scanStopFreq - scanStartFreq) / segmentSpan);
    
    fprintf('\n======================================================\n');
    fprintf('STARTING LIVE SEGMENTED DIP DETECTION\n');
    fprintf('Total Segments: %d\n', totalSegments);
    fprintf('Streaming raw data to:  %s\n', filename);
    fprintf('======================================================\n');
    
    % Setup plot window
    figure(1); clf; hold on; grid on;
    xlabel('Frequency (MHz)'); ylabel('Amplitude (dBm)');
    title('Live Segmented Dip Detection');
    xlim([scanStartFreq, scanStopFreq] / 1e6);
    
    currentSegNum = 1;
    
    % Segmented Survey Loop
    for segStart = scanStartFreq : segmentSpan : (scanStopFreq - segmentSpan)
        segStop = segStart + segmentSpan;
        
        fprintf('\n--- Processing Segment %d of %d [%.3f to %.3f MHz] ---\n', ...
            currentSegNum, totalSegments, segStart/1e6, segStop/1e6);
        
        writeline(ssa, ':SENSe:FREQuency:STARt ' + string(segStart));
        writeline(ssa, ':SENSe:FREQuency:STOp ' + string(segStop));
        
        writeline(ssa, ':INITiate:CONTinuous ON');
        writeline(ssa, ':AVERage:CLEar');
        writeline(ssa, ':INITiate:IMMediate');
        
        % Polling Loop
        currentAvg = 0;
        icheck = 0;
        lastAvg = -1;
        stuckCounter = 0;
        
        fprintf('STATUS: Hardware sweeping. Waiting for %d averages...\n', targetAvg);
        
        while currentAvg < targetAvg || icheck < 3
            pause(2.0);
            writeline(ssa, ':SENSe:AVERage:TRACe1?'); 
            response = readline(ssa);
            
            if ~isempty(response)
                currentAvg = str2double(response);
                if isnan(currentAvg), currentAvg = 0; end
            end
            
            % Check if instrument frozen
            if currentAvg == lastAvg
                stuckCounter = stuckCounter + 1;
            else
                stuckCounter = 0;
            end
            
            lastAvg = currentAvg;
            
            % If stuck for ~40 seconds (20 checks * 2 sec)
            if stuckCounter > 20
                fprintf('\n[!!!] WATCHDOG TRIGGERED [!!!]\n');
                fprintf('Instrument hung at %d averages. Forcing data download to save segment...\n', currentAvg);
                break;
            end
            
            icheck = icheck + 1;
        end
        
        % Freeze and download
        writeline(ssa, ':INITiate:CONTinuous OFF');
        pause(0.5);
        
        fprintf('STATUS: Downloading trace data matrix...\n');
        writeline(ssa, 'TRACe:DATA? 1');
        traceData = str2double(split(readline(ssa), ','));
        traceData = traceData(:).'; % Force to row vector
        
        freqAxis = linspace(segStart, segStop, length(traceData));
        
        % Stream to disk
        M = [freqAxis(:), traceData(:)];
        fprintf(fileID, '%.6f\t%.6f\n', M.');
        
        % Dip Detection Algorithm
        fprintf('STATUS: Running dip detection algorithm...\n');
        
        smoothedData = movmedian(traceData, 3);
        digitalSpacingHz = segmentSpan / (length(traceData) - 1);
        minDipPoints = floor((expectedWidthHz / digitalSpacingHz) * 0.5);
        
        baselineWindowPoints = round((expectedWidthHz / digitalSpacingHz) * 10);
        if baselineWindowPoints > floor(length(traceData) / 2)
            baselineWindowPoints = floor(length(traceData) / 2);
        end
        if baselineWindowPoints < 10, baselineWindowPoints = 10; end
        
        localBaseline = movmedian(smoothedData, baselineWindowPoints);
        isBelowThreshold = smoothedData < (localBaseline - dipThresholdDB);
        
        paddedArray = [0, isBelowThreshold, 0];
        edges = diff(paddedArray);
        startIndices = find(edges == 1);
        endIndices = find(edges == -1) - 1;
        dipWidths = endIndices - startIndices + 1;
        
        validDipIDs = find(dipWidths >= minDipPoints);
        
        allValidIndices = [];
        if isempty(validDipIDs)
            fprintf('  -> RESULT: Segment clear. No dips > %.1f dB found.\n', dipThresholdDB);
        else
            fprintf('  -> RESULT: >>> ALERT! %d dip(s) detected! <<<\n', length(validDipIDs));
            for k = 1:length(validDipIDs)
                thisDip = validDipIDs(k);
                dipIndices = startIndices(thisDip) : endIndices(thisDip);
                allValidIndices = [allValidIndices, dipIndices];
                
                centerOfDip = freqAxis(round(median(dipIndices))) / 1e6;

                dipDepth = max(localBaseline(dipIndices) - smoothedData(dipIndices));
                
                fprintf('      - Dip found roughly at: %.4f MHz (Depth: %.2f dB)\n', centerOfDip, dipDepth);
                fprintf(dipFileID, '%.6f, %.2f\n', centerOfDip, dipDepth);
            end
        end
        
        % Plotting
        fprintf('STATUS: Updating live plot...\n');
        plot(freqAxis / 1e6, traceData, 'b', 'LineWidth', 1);
        if ~isempty(allValidIndices)
            plot(freqAxis(allValidIndices) / 1e6, traceData(allValidIndices), 'r.', 'MarkerSize', 10);
        end
        drawnow;
        
        currentSegNum = currentSegNum + 1;
        
    end % End of the Segmented Survey Loop
    
    fclose(fileID);
    fclose(dipFileID);
    fprintf('\n======================================================\n');
    fprintf('SURVEY COMPLETE!\n');
    fprintf('All requested segments analyzed and safely secured in:\n%s\n', filename);
    fprintf('======================================================\n');
    
catch ME
    fprintf('\n======================================================\n');
    fprintf('!!! CRITICAL SYSTEM CRASH DETECTED !!!\n');
    fprintf('Error Log: %s\n', ME.message);
    fprintf('Attempting to safely close file handles and sockets...\n');
    
    if exist('fileID', 'var') && fileID ~= -1
        fclose(fileID);
        fprintf('File handle secured.\n');
    end

    if exist('dipFileID', 'var') && dipFileID ~= -1
        fclose(dipFileID);
        fprintf('Dip log handle secured.\n');
    end

    if exist('ssa', 'var')
        clear ssa;
        fprintf('TCP socket severed.\n');
    end
    fprintf('======================================================\n');
end