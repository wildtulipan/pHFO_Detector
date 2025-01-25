function [FastRipples] = FastRippleDetector(LFP, Fs, BandPass, stdsRMS, stdsPeaks)
    % FastRippleDetector Detect High-Frequency Oscillations (HFOs) in Local Field Potential signals
    %
    % Adapted from: 
    % Staba, R. J., Wilson, C. L., Bragin, A., Fried, I., & Engel Jr, J. (2002). 
    % Quantitative analysis of high-frequency oscillations (80–500 Hz) recorded in human epileptic hippocampus and entorhinal cortex. 
    % Journal of neurophysiology, 88(4), 1743-1752.
    %
    % Inputs:
    %   LFP         - Input local field potential signal
    %   Fs          - Sampling frequency
    %   BandPass    - Bandpass filter coefficients
    %   stdsRMS     - RMS threshold multiplier (default: 5)
    %   stdsPeaks   - Peak threshold multiplier (default: 3)
    %
    % Outputs:
    %   HFOsStaba   - Indices of detected HFO events

    % Set default parameters
    if nargin < 4 || isempty(stdsRMS), stdsRMS = 5; end
    if nargin < 5 || isempty(stdsPeaks), stdsPeaks = 3; end

    % Detection parameters
    minLength = round(0.006 * Fs);     % Minimum event duration (6 ms)
    maxDistance = round(0.010 * Fs);   % Maximum gap between event segments (10 ms)
    minRequiredPeaks = 6;              % Minimum peaks per event

    % Apply band-pass filter
    LFP_BP = filtfilt(BandPass, LFP);

    % Calculate moving RMS
    RMSpoints = round(0.003 * Fs);
    movingRMS = movmean(abs(LFP_BP), RMSpoints);

    % Set thresholds
    thresholdRMS = mean(movingRMS) + (std(movingRMS) * stdsRMS);
    peakthreshold = mean(abs(LFP_BP)) + (std(abs(LFP_BP)) * stdsPeaks);

    % Find RMS threshold-crossing indices
    overthresholdIdx = find(movingRMS >= thresholdRMS);
    
    % Return empty if no events detected
    if isempty(overthresholdIdx)
        FastRipples = [];
        return;
    end

    % Identify consecutive segments
    consecutiveValues = [];
    startValue = overthresholdIdx(1);
    endValue = overthresholdIdx(1);

    for i = 2:length(overthresholdIdx)
        if overthresholdIdx(i) == overthresholdIdx(i-1) + 1
            endValue = overthresholdIdx(i);
        else
            if endValue - startValue + 1 >= minLength
                consecutiveValues = [consecutiveValues; startValue, endValue];
            end
            
            startValue = overthresholdIdx(i);
            endValue = overthresholdIdx(i);
        end
    end

    % Add last segment if valid
    if endValue - startValue + 1 >= minLength
        consecutiveValues = [consecutiveValues; startValue, endValue];
    end

    % Combine nearby segments
    while true
        combined = false;
        for i = 1:size(consecutiveValues, 1) - 1
            distance = consecutiveValues(i+1, 1) - consecutiveValues(i, 2) - 1;
            if distance <= maxDistance
                consecutiveValues(i, 2) = consecutiveValues(i+1, 2);
                consecutiveValues(i+1, :) = [];
                combined = true;
                break;
            end
        end
        if ~combined
            break;
        end
    end

    % Validate events by peak characteristics
    FastRipples = [];
    for i = 1:size(consecutiveValues, 1)
        startIdx = consecutiveValues(i, 1);
        endIdx = consecutiveValues(i, 2);
        
        signalChunk = LFP(startIdx:endIdx);
        rectifiedSignal = abs(signalChunk);
        
        [peaks, ~] = findpeaks(rectifiedSignal, 'MinPeakHeight', peakthreshold);
        
        if length(peaks) >= minRequiredPeaks
            FastRipples = [FastRipples; consecutiveValues(i, :)];
        end
    end
end