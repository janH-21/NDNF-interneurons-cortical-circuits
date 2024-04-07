function [spikeDelay_ms, f] = spikeDelay(allData, cellNo, tStepStart, si_us, plotBool)
% INPUTs
%   allData         -   allData structure from barrage master script
%   cellNo          -   ID# of cell to analyze
%   tStepStart      -   start time of current injection [ms]
%   si_us           -   sampling interval [µs]
%   plotBool        -   boolan. plot sanity check plot?
% OUTPUTs
%   spikeDelay_ms   -   spike delay [ms]

    % find suprathreshold sweep
    sweepInds = characterizationPlot(allData, cellNo);
    supraThresholdSweep = sweepInds(3);
    
    % calculate delay to first spike
    firstTime = allData(cellNo).analysis.spikeTimes_ms.cell{supraThresholdSweep}(1);
    firstTime = (firstTime - (supraThresholdSweep-1)*4000); % substract previous sweeps converted to ms
    spikeDelay_ms = firstTime - tStepStart;
    
    % sanity check plot
    if plotBool == true
        f = figure;
        plot((1:length(allData(cellNo).recordingData.rawData_mV(:,sweepInds(3))))*si_us/1000,...
             allData(cellNo).recordingData.rawData_mV(:,sweepInds(3)),...
             'Color', 'black'); hold on
        line([firstTime firstTime], [ylim], 'Color', 'red', 'LineStyle', '-.');
                title('spike delay sanity check'); 
                xlabel('time [ms]'); 
                ylabel('V_m [mV]'); 
                set(gca, 'FontSize', 12)
    end
    
end