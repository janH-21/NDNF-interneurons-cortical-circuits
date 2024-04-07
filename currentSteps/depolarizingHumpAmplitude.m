function [humpVal, humpAmpl, humpPercentAmpl, humpTs, f] = depolarizingHumpAmplitude(allData, cellNo, baseTWind, maxTWind, si_us, plotBool)
% INPUTs
%   allData             -   allData structure from barrage master script
%   cellNo              -   ID# of cell to analyze
%   baseTWind           -   time window in which baseline during current 
%                           injection should be calculated [ms ms]
%   maxTWind            -   time window in which peak should be detected [ms ms]
%   si_us               -   sampling interval [µs]
%   plotBool            -   boolan. plot sanity check plot?
% OUTPUTs
%   depolarizingHump    -   amplitude of depolarizing hump [ms]. CAVE: will
%                           alwas calculate an amplitude and can thus 
%                           return noise if no hump is present

        % find subthreshold sweep
        sweepInds = characterizationPlot(allData, cellNo);
        subThreshSweep = allData(cellNo).recordingData.rawData_mV(:,sweepInds(2));

        % calculate hump amplitude
        baseTWind = baseTWind./si_us.*1000;
        maxTWind = maxTWind./si_us.*1000 - baseTWind(1) + 1;
        step = detrend(subThreshSweep(baseTWind(1):baseTWind(2))) - abs(mean(subThreshSweep(baseTWind(1):baseTWind(2))));  % = detrend()
        stepBaseline = median(step);
        sweepBaseline = median(subThreshSweep(1:1000*1000/20));
        [humpVal, humpInd] = max(step(maxTWind(1):maxTWind(2)));
        humpTs = (humpInd-1)*20/1000 + 1062;
        humpAmpl = humpVal - stepBaseline;
        humpPercentAmpl = (humpVal - sweepBaseline)/(stepBaseline - sweepBaseline) * 100 - 100;
        
        % sanity check plot
        if plotBool == true
            f = figure;
            plot((1:length(subThreshSweep))*20/1000, subThreshSweep, 'Color', 'black'); hold on; grid minor
            plot(1062:0.02:2062, step, 'Color', 'green'); ylim([-100 50])
            line([xlim], [stepBaseline stepBaseline], 'LineStyle', '-.', 'Color', 'blue', 'LineWidth', 1)
            line([xlim], [sweepBaseline sweepBaseline], 'LineStyle', '-.', 'Color', 'blue', 'LineWidth', 1)
            line(([maxTWind(2) maxTWind(2)]+baseTWind(1)-1)/1000*20, [ylim], 'LineStyle', '-', 'Color', 'magenta', 'LineWidth', 1)
            plot(humpTs, humpVal, '*', 'MarkerSize', 5, 'Color', 'red')
            title('depolarizing hump sanity check')
            xlabel('time [ms]'); 
            ylabel('V_m [mV]'); 
            set(gca, 'FontSize', 12)
            hold off
        end  
end