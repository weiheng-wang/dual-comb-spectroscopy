% SIGLENT SSA3021X: High-Resolution Scan of Dips

% Parameters
analyzerIP = '192.168.0.50';
port = 5025;
timeout = 120; 

dipLogFilename = 'dip_locations.txt'; % CHANGE THIS TO EXACT FILE NAME
segmentSpan    = 4e3;                 % 4 kHz segment width
targetAvg      = 100;                 % 100 averages
analogRBW      = 3;                   % 3 Hz RBW

try
    fprintf('STATUS: Loading dip coordinates from %s...\n', dipLogFilename);
    % readmatrix skips the text header row
    dipData = readmatrix(dipLogFilename);
    
    if isempty(dipData)
        error('The dip log file is empty or could not be read.');
    end
    
    targetFreqs_MHz = dipData(:, 1);
    numDips = length(targetFreqs_MHz);
    fprintf('STATUS: Successfully loaded %d target dips.\n\n', numDips);
    
    fprintf('STATUS: Checking TCP connection...\n');
    if ~exist('ssa', 'var')
        ssa = tcpclient(analyzerIP, port, 'Timeout', timeout);
        flush(ssa); 
        writeline(ssa, '*IDN?');
        fprintf('STATUS: Connected -> %s\n', readline(ssa));
    else
        try
            flush(ssa); 
            writeline(ssa, '*IDN?');
            readline(ssa);
            fprintf('STATUS: Reusing existing TCP connection.\n');
        catch
            clear ssa; pause(2.0);
            ssa = tcpclient(analyzerIP, port, 'Timeout', timeout);
            flush(ssa);
            writeline(ssa, '*IDN?');
            fprintf('STATUS: Connected -> %s\n', readline(ssa));
        end
    end
    
    % Instrument Setup
    fprintf('STATUS: Configuring sweep parameters...\n');
    writeline(ssa, ':INSTrument SANalyzer');
    writeline(ssa, ':SENSe:BWIDth:RESolution ' + string(analogRBW));
    
    writeline(ssa, ':TRACe1:MODE AVERage');
    writeline(ssa, ':AVERage:TRACe1:STATe ON');
    writeline(ssa, ':AVERage:TRACe1:COUNt ' + string(targetAvg));
    writeline(ssa, ':DETector:TRACe1 AVERage');
    
    % Setup plot window
    figure(1); clf; grid on;
    xlabel('Frequency (MHz)'); ylabel('Amplitude (dBm)');
    title('Live High-Resolution Scan');
    
    fprintf('\n======================================================\n');
    fprintf('STARTING DIP EXTRACTION\n');
    fprintf('======================================================\n');
    
    % Scanning Loop
    for i = 1:numDips
        centerFreq_Hz = targetFreqs_MHz(i) * 1e6;
        segStart = centerFreq_Hz - (segmentSpan / 2);
        segStop  = centerFreq_Hz + (segmentSpan / 2);
        
        fprintf('\n--- Extracting Dip %d of %d [Center: %.4f MHz] ---\n', i, numDips, targetFreqs_MHz(i));
        
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
        
        fprintf('STATUS: Sweeping at 3 Hz RBW. Waiting for %d averages...\n', targetAvg);
        
        while currentAvg < targetAvg || icheck < 3
            pause(3.0); % Pause longer because 3Hz sweeps take time
            writeline(ssa, ':SENSe:AVERage:TRACe1?'); 
            response = readline(ssa);
            
            if ~isempty(response)
                currentAvg = str2double(response);
                if isnan(currentAvg), currentAvg = 0; end
            end
            
            if currentAvg == lastAvg
                stuckCounter = stuckCounter + 1;
            else
                stuckCounter = 0;
            end
            
            lastAvg = currentAvg;
            
            % If stuck for ~90 seconds (30 checks * 3 sec)
            if stuckCounter > 30
                fprintf('\n[!!!] WATCHDOG TRIGGERED [!!!]\n');
                fprintf('Instrument hung. Forcing data download...\n');
                break;
            end
            icheck = icheck + 1;
        end
        
        % Freeze and download
        writeline(ssa, ':INITiate:CONTinuous OFF');
        pause(0.5);
        
        fprintf('STATUS: Downloading trace...\n');
        writeline(ssa, 'TRACe:DATA? 1');
        traceData = str2double(split(readline(ssa), ','));
        traceData = traceData(:).';
        
        freqAxis = linspace(segStart, segStop, length(traceData));
        
        % Plotting
        plot(freqAxis / 1e6, traceData, 'b', 'LineWidth', 1.5);
        xlim([segStart, segStop] / 1e6);
        title(sprintf('Dip %d: Center %.4f MHz', i, targetFreqs_MHz(i)));
        drawnow;
        
        % Creates a file like: dip_01_at_38.502MHz.txt
        % CHANGE FILE NAME IF USING SAME DIP LOCATIONS TO AVOID OVERWRITING
        % Ex: scanning a gas cell and reference using the same dip locations
        outFilename = sprintf('dip_%02d_at_%.4fMHz.txt', i, targetFreqs_MHz(i));
        fileID = fopen(outFilename, 'w');
        fprintf(fileID, 'RF_Frequency_MHz, Amplitude_dBm\n');
        
        M = [freqAxis(:) / 1e6, traceData(:)];
        fprintf(fileID, '%.6f, %.6f\n', M.');
        fclose(fileID);
        
        fprintf('STATUS: Successfully saved to %s\n', outFilename);
        
    end
    
    fprintf('\n======================================================\n');
    fprintf('ALL EXTRACTIONS COMPLETE!\n');
    fprintf('======================================================\n');
    
catch ME
    fprintf('\n======================================================\n');
    fprintf('!!! CRITICAL SYSTEM CRASH DETECTED !!!\n');
    fprintf('Error Log: %s\n', ME.message);
    if exist('fileID', 'var') && fileID ~= -1
        fclose(fileID);
    end
    if exist('ssa', 'var')
        clear ssa;
    end
    fprintf('======================================================\n');
end