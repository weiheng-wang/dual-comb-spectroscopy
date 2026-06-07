% SIGLENT SSA3021X: Automated Comb Profile Sweep

% Parameters
analyzerIP = '192.168.0.50';
port = 5025;
timeout = 60;

broadStartFreq = 5e6;
broadStopFreq  = 25e6;
analogRBW       = 1000;        % 1 kHz
snr_margin     = 26;           % dB drop from absolute peak for cutoff

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
    
    % Instrument Setup
    fprintf('Configuring Sweep...\n');
    writeline(ssa, ':INSTrument SANalyzer');
    writeline(ssa, ':SENSe:FREQuency:STARt ' + string(broadStartFreq));
    writeline(ssa, ':SENSe:FREQuency:STOp ' + string(broadStopFreq));
    writeline(ssa, ':SENSe:BWIDth:RESolution ' + string(analogRBW));
    writeline(ssa, ':DISPlay:WINDow:TRACe:Y:SCALe:RLEVel -40');
    
    writeline(ssa, ':TRACe1:MODE AVERage');
    writeline(ssa, ':AVERage:TRACe1:STATe ON');
    writeline(ssa, ':AVERage:TRACe1:COUNt 1'); % Only need 1 full sweep
    writeline(ssa, ':DETector:TRACe1 AVERage');
    
    writeline(ssa, ':INITiate:CONTinuous ON');
    writeline(ssa, ':AVERage:CLEar');
    writeline(ssa, ':INITiate:IMMediate');
    
    % Polling Loop
    currentAvg = 0;
    icheck = 0;
    lastAvg = -1;
    stuckCounter = 0;
    
    fprintf('Sweeping... waiting for hardware confirmation...\n');
    
    % Waits until 1 sweep is finished
    while currentAvg < 1 || icheck < 3
        pause(1.0); % Poll every 1 second
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
        
        % If stuck for ~120 seconds (120 checks * 1 sec)
        if stuckCounter > 120
            fprintf('\n[!!!] WATCHDOG TRIGGERED [!!!]\n');
            fprintf('Instrument hung. Forcing data download...\n');
            break;
        end
        
        icheck = icheck + 1;
    end
    
    % Stop continuous sweep
    writeline(ssa, ':INITiate:CONTinuous OFF');
    pause(0.5);
    
    % Download data
    fprintf('Transferring data...\n');
    writeline(ssa, 'TRACe:DATA? 1');
    traceData = str2double(split(readline(ssa), ','));
    traceData = traceData(:).'; % Force to row vector
    
    numPoints = length(traceData);
    rf_freq_Hz = linspace(broadStartFreq, broadStopFreq, numPoints);
    rf_freq_MHz = rf_freq_Hz / 1e6;
    mag_dB = traceData;
    
    % Comb Profile Algorithm
    max_power_dB = max(mag_dB);
    dynamic_threshold = max_power_dB - snr_margin;
    valid_mask = mag_dB >= dynamic_threshold;
    
    freq_window_Hz = rf_freq_Hz(valid_mask);
    freq_window_MHz = rf_freq_MHz(valid_mask);
    mag_window_dB = mag_dB(valid_mask);
    
    if isempty(freq_window_Hz)
        error('No data found above %.1f dB. Is the laser firing?', dynamic_threshold);
    end
    
    recommendedStart = min(freq_window_Hz);
    recommendedStop  = max(freq_window_Hz);
    
    % Plotting
    figure('Name', 'Comb Profile', 'Position', [100, 100, 900, 500]);
    clf; hold on;
    
    stem(rf_freq_MHz, mag_dB, 'Marker', 'none', 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'BaseValue', -120);
    stem(freq_window_MHz, mag_window_dB, 'Marker', 'none', 'Color', '#0072BD', 'LineWidth', 1.5, 'BaseValue', -120);
    
    yline(dynamic_threshold, 'r--', 'LineWidth', 2, ...
        'Label', sprintf('Dynamic Cutoff: %.1f dB', dynamic_threshold), ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
        
    title('Automated Comb Profile Extraction');
    xlabel('RF Frequency (MHz)');
    ylabel('Magnitude (dB)');
    legend('Discarded Noise Floor', 'Valid Comb Profile', 'Location', 'northeast');
    grid on;
    ylim([-115, max_power_dB + 5]); 
    xlim([min(rf_freq_MHz), max(rf_freq_MHz)]);
    hold off;
    
    % Output
    fprintf('\n======================================================\n');
    fprintf('SCAN COMPLETE.\n');
    fprintf('Valid Comb Envelope Identified.\n');
    fprintf('======================================================\n');
    
    fprintf('scanStartFreq  = %.0f;      %% %.3f MHz\n', recommendedStart, recommendedStart/1e6);
    fprintf('scanStopFreq   = %.0f;      %% %.3f MHz\n', recommendedStop, recommendedStop/1e6);
    fprintf('\n======================================================\n');

catch ME
    fprintf('\n!!! CRASH DETECTED !!!\nError: %s\n', ME.message);
    if exist('ssa', 'var'), clear ssa; end
end