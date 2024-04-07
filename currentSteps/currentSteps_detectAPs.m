function APstats = barrage_detectAPs(d, si, t, n, cutOffmVpms, tLim, tSpike, tSpikeMax, plotXLim)
    APstats = struct;
    % create empty figure
    %no plot!	f = figure('units','normalized','outerposition',[0 0 1 1]);
    % create struct to hold data
    for jter = 1:size(n,2)
        [jter size(n,2)]
        % extract data
        realData = d(:, n(jter));
        [firstD, secondD] = alculateDerivative(t, realData);
        
        % prepare figure
        %no plot!	a = axis; hold on; 
        %no plot!	xlim(plotXLim); ylim([-100, 60]); 
        %no plot!	xlabel('time [ms]'); ylabel('V_{m} [mV]');
        %no plot!	set(gca,'fontsize',20)

        % plot data
        %no plot!	plot(t,realData, 'LineWidth', 1.5, 'Color', 'black')

        % find and plot max value of real data
        maxReal = max(realData);
        %no plot!	line([0 max(t)],[maxReal maxReal], 'LineStyle', '--');
             
        % find AP thresholds
        spkCount = 0; % count APs
        APthreshTs = [];   % find AP threshold time stamps
        APthresVal = [];   % find corresponding potentials
        APmaxVal = []; % collect AP max vals
        APmaxTs = []; % collect corresponting ts
        APmaxSlopeVal = []; % collect maximum spike slopes
        APmaxSlopeTs = []; % collect time stamp of maximum spike slopes
        AHPminVal = []; % collect AP max vals
        AHPminTs = []; % collect corresponting ts
        APhalfVal = [];
        APhalfCrossT1 = [];
        APhalfCrossT2 = [];
        APhalfWidth = []; % collect spike half-widths
        AHPhalfCrossT1 = [];
        AHPhalfCrossT2 = [];
        AHPhalfWidth = [];
        AHPhalfVal = [];
        APind = 0;
        tmpMaxInd = 0;
        for iter = 2:length(firstD) % iterate first derivative
            if firstD(iter) > cutOffmVpms && firstD(iter-1) <= cutOffmVpms % primary detection threshold in dVm/dt
                if ~(iter > length(firstD)-(tSpike/si)) % avoid end of recording artifacts
                    [peakData, ~] = findpeaks(realData(iter:iter+75)); 
                    [~, valleyIndcs] = findpeaks(-realData(iter:iter+75));
                    if ~isempty(peakData) && (isempty(valleyIndcs) || valleyIndcs(1) > 20) 
                    % SANITY CHECKS TO DIFFERENTIATE SPIKES FROM NOISE
                    % ~isempty(peakData)            --> if any peak at all --> catch errors 
                    % isempty(valleyIndcs)          --> if any valley at all --> good if no valleys OR                   
                    % valleyIndcs(1) > 10           --> if valley more than a few microseconds in the future
                    if peakData(1) > realData(iter)+15  % cancel out noise
                    % if max(realData(iter:iter+75)) > realData(iter)+15   % cancel out noise
                        if iter > tmpMaxInd+50 % require 1ms refractory period 
                            % get AP time & plot it
                            spkCount = spkCount + 1; % increase spike count
                            APind = iter;
                            APthreshTs(spkCount) = iter*si/1000; % save time stamp
                            APthresVal(spkCount) = realData(iter); % save potential
                            
%                             if APthreshTs(spkCount) > tLim(1) && APthreshTs(spkCount) <= tLim(2)
%                                 line([APthreshTs(spkCount) APthreshTs(spkCount)],[-100 60], 'Color', 'green', 'LineWidth', 1, 'LineStyle', '-.')
%                             else
%                                 line([APthreshTs(spkCount) APthreshTs(spkCount)],[-100 60], 'Color', 'cyan', 'LineWidth', 1, 'LineStyle', '-.')
%                             end
                            % extract time window after threshold crossing
                            APts = [iter:iter+(tSpike/si)-1];
                            APtsMax = [iter:iter+(tSpikeMax/si)-1];
                            % find local maximum                
                            %[maxVal, maxInd] = max(realData(APtsMax));
                            [maxVal, maxInd] = findpeaks(realData(APts));
                            APmaxVal(spkCount) = maxVal(1);
                            APmaxTs(spkCount) = (iter + maxInd(1) - 1)*20/1000;
                            plot(APmaxTs(spkCount), APmaxVal(spkCount), 'r*', 'Color', 'red', 'MarkerSize', 5)
                            tmpMaxInd = iter+maxInd(1)-1;
                            
                            % find maximal spike slope
                            [APmaxSlopeVal(spkCount) APmaxSlopeInd] = max(firstD(iter:iter+(0.5/20*1000)));
                            APmaxSlopeTs(spkCount) = (iter + APmaxSlopeInd - 1) * 20/1000;
%                             line([APmaxSlopeTs(spkCount) APmaxSlopeTs(spkCount)], [ylim], 'LineStyle', '-.', 'Color', 'magenta', 'LineWidth', 1)
                            
                            % find local minimum
% %                             tmpAHPminTs = APmaxTs(spkCount);
% %                             tmpAHPminVal = APmaxVal(spkCount);
% %                             for kter = APts
% %                                 % find AHP min
% %                                 if kter > APmaxTs & realData(kter) < tmpAHPminVal
% %                                     tmpAHPminTs = kter;
% %                                     tmpAHPminVal = realData(kter);
% %                                 end
% %                             end
% %                             AHPminTs(spkCount) = tmpAHPminTs*20/1000;
% %                             AHPminVal(spkCount) = tmpAHPminVal;
                            [AHPvals, AHPlocs] = findpeaks(abs(realData),...
                                                           'MinPeakDistance', 5,... 
                                                           'MinPeakWidth', 5,...
                                                           'MinPeakProminence', 0.5,...
                                                           'Annotate', 'extents');
                            peaksLargerMaxTs = find(AHPlocs > APmaxTs(spkCount)*1000/20);  
                            if ~isempty(peaksLargerMaxTs) % exception for borderline barrage peaks at end of sweep
                                AHPminInd = AHPlocs(peaksLargerMaxTs(1));
                                AHPminTs(spkCount) = AHPminInd*20/1000;
                                AHPminVal(spkCount) = -AHPvals(peaksLargerMaxTs(1));
%no plot!                                         plot(AHPminTs(spkCount), AHPminVal(spkCount), 'r*', 'Color', 'blue', 'MarkerSize', 5)
                            else
                                AHPminTs(spkCount) = NaN;
                                AHPminVal(spkCount) = NaN;
                            end
                            % find AP half-width
                            APhalfVal(spkCount) = (APthresVal(spkCount) + APmaxVal(spkCount))/2;
                            [maxVals, maxIndcs] = findpeaks(-abs(realData(APts)-APhalfVal(spkCount)));
                            if maxIndcs(1) < maxIndcs(2)
                                APhalfCrossInd1 = maxIndcs(1);
                                APhalfCrossInd2 = maxIndcs(2);
                                APhalfCrossVal1 = maxIndcs(1);
                                APhalfCrossVal2 = maxIndcs(2);
                            else
                                APhalfCrossInd1 = maxIndcs(2);
                                APhalfCrossInd2 = maxIndcs(1);
                                APhalfCrossVal1 = maxIndcs(2);
                                APhalfCrossVal2 = maxIndcs(1);
                            end
                            APhalfWidth(spkCount) = (APhalfCrossInd2 - APhalfCrossInd1)*20/1000;
                            APhalfCrossT1(spkCount) = (iter + APhalfCrossInd1 - 1)*20/1000;
                            APhalfCrossT2(spkCount) = (iter + APhalfCrossInd2 - 1)*20/1000;
%                             line([APhalfCrossT1(spkCount) APhalfCrossT2(spkCount)],...
%                                  [APhalfVal(spkCount) APhalfVal(spkCount)],...
%                                  'LineStyle', '-.', 'Color', 'cyan', 'LineWidth', 1)
                            
                            % find AHP half width
                            if ((jter == 1) || (APstats.nSpikesDuringStep{jter-1} <= 30)) && ~isempty(peaksLargerMaxTs)
                                AHPhalfVal(spkCount) = (APthresVal(spkCount) + AHPminVal(spkCount))/2;
                                [maxVals, maxIndcs] =  findpeaks(-abs(realData-AHPhalfVal(spkCount)),...
                                                                 'MinPeakProminence', 1,...
                                                                 'MinPeakDistance', 5,...
                                                                 'Annotate', 'extents');
                                peaksAfterSpike = find(maxIndcs > iter);
                                if length(peaksAfterSpike)>1 && ~(maxIndcs(peaksAfterSpike(2))*20/1000 > 2062)   
                                    AHPhalfCrossT1(spkCount) = (maxIndcs(peaksAfterSpike(1)))*20/1000;
                                    AHPhalfCrossT2(spkCount) = (maxIndcs(peaksAfterSpike(2)))*20/1000;
                                    AHPhalfWidth(spkCount) = AHPhalfCrossT2(spkCount) - AHPhalfCrossT1(spkCount);
%                                     line([AHPhalfCrossT1(spkCount) AHPhalfCrossT2(spkCount)],...
%                                          [AHPhalfVal(spkCount) AHPhalfVal(spkCount)],...
%                                          'LineStyle', '-.', 'Color', 'cyan', 'LineWidth', 1)
                                end
                            end % else AHP half width parameters remain empty, as detection is not reliable anymore with high frequencies & additionally biologically not plausible
                        end
                    end
                    end
                end
            end
        end
        
        % save data
        APstats.nSpikesDuringStep{jter} = sum((APthreshTs>tLim(1) & APthreshTs<tLim(2)));
        APstats.nSpikesTotal{jter} = spkCount;
        APstats.threshTs{jter} = APthreshTs;
        APstats.threshVal{jter} = APthresVal;
        APstats.APmaxVal{jter} = APmaxVal;
        APstats.APmaxTs{jter} = APmaxTs;
        APstats.APamplitude{jter} = APmaxVal - APthresVal;
        APstats.APmaxSlopeVal{jter} = APmaxSlopeVal;
        APstats.APmaxSlopeTs{jter} = APmaxSlopeTs;
        APstats.AHPminVal{jter} = AHPminVal;
        APstats.AHPminTs{jter} = AHPminTs;
        APstats.AHPamplitude{jter} = APthresVal - AHPminVal; 
        APstats.peakToTrough{jter} = AHPminTs - APmaxTs;
        APstats.APhalfWidthTs1{jter} = APhalfCrossT1;
        APstats.APhalfWidthTs2{jter} = APhalfCrossT2;
        APstats.APhalfWidth{jter} = APhalfWidth;
        APstats.APhalfVal{jter} = APhalfVal;
        APstats.AHPhalfWidthTs1{jter} = AHPhalfCrossT1;
        APstats.AHPhalfWidthTs2{jter} = AHPhalfCrossT2;
        APstats.AHPhalfWidth{jter} = AHPhalfWidth;
        APstats.AHPhalfVal{jter} = AHPhalfVal;
        
        % add figure title
        title1 = strjoin({'Sweep #', num2str(n(jter)), ':'},'');
        title2 = strjoin({num2str(spkCount), ' action potentials detected'},'');
        title({title1, title2})
        
        % give inspection time, then close window
        %no plot!	pause(0.1);
        %no plot!	delete(gca);
    end
    close gcf;
end