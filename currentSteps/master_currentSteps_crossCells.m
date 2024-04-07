%% prepare workspace
clear all; close all; clc; 
cd('/main/')
addpath('./analysis/currentSteps/')

%% --------------------- USER EDITABLE BLOCK --------------------------- %%
% set working directory
path = './analysis/currentSteps_NDNF/';
% path = './analysis/currentSteps_NDNF_synBlockers/';
% path = './analysis/currentSteps_FSIN/';

regexPattern = '_barrageDataSingleCell.mat';
% regexPattern = '_synBlockers_barrageDataSingleCell.mat';

% NOTE: 
% 1. Change if conditions in lines 36++ if working with non-NDNF cells.
% In particular, some pyramidal cells and some non-FS INs don't actually 
% reach FRs > 19Hz, such that this condition as well as downstream
% follow-ups need to be changed.
% 2. Not all parameters extracted will be used in the subsequent analysis,
% and non-used parameters were not consistently checked for bugs.

%---------------------END USER EDITABLE BLOCK-----------------------------%


%% -------------------- read & combine data ---------------------------- %%
% make list of all matching files in directory
lDir = getFilesFromDir(path, regexPattern, 'names');

% prepare vars to hold data
loadData = struct();
allData = struct();
nSweeps = zeros(1,length(lDir));

% add data to allData struct
counter = 0;
for iter = 1:length(lDir)
    loadData(iter).data = load(strcat(path, lDir(iter)));
    
    if (loadData(iter).data.saveData.cellType == "NDNF+" ||...
        loadData(iter).data.saveData.cellType == "NDNF+/NPY+") &&...
        max([loadData(iter).data.saveData.analysis.FR_Hz_step]) > 19
%        loadData(iter).data.saveData.cellType == "L23_FS_IN") 

        counter = counter + 1;
        allData(counter).recordingFileName = loadData(iter).data.saveData.recordingFileName;
        allData(counter).cellType = loadData(iter).data.saveData.cellType;
        allData(counter).recordingData = loadData(iter).data.saveData.recordingData;
        allData(counter).analysis = loadData(iter).data.saveData.analysis;
        nSweeps(counter) = length(allData(counter).recordingData.Isteps_pA);

        if isfield(loadData(iter).data.saveData, 'paradigm')
            allData(counter).paradigm = loadData(iter).data.saveData.paradigm;
        else
            allData(counter).paradigm.paradigmName = "basic";
        end
    else
        {loadData(iter).data.saveData.cellType,...
         num2str(max([loadData(iter).data.saveData.analysis.FR_Hz_step])),...
         loadData(iter).data.saveData.recordingFileName}
    end
    
end

% extract meta info - [pA], NOTE: THIS ASSUMES SAME STEP SIZE IN ALL RECORDINGS
[nSweepsMax, maxSweepInd] = max(nSweeps);
nCells = counter;
stepCenter = mean(allData(1).analysis.stepTimes_ms(1,:));
stepSize = allData(1).recordingData.stepSize_pA; 
sweepDuration = allData(1).recordingData.sweepDuration;
interSweepDuration = allData(1).recordingData.interSweepTime;
totalDuration = sweepDuration + interSweepDuration;
si_us = allData(counter).recordingData.si_us;


%% ------------------------- calculate stats --------------------------- %%
for iter = 1:length(allData)
    
    % create time stamps
    ts = (1:length(allData(iter).recordingData.rawData_mV(:, 1)))*20/1000;

    % find first sweep with AP & subthreshold sweep
    sweepInds = characterizationPlot(allData, iter);
    allData(iter).analysis.supraThresholdSweep = sweepInds(3);
    
    % input resistance
    % extract data
    subThreshSweepNos = 1:allData(iter).analysis.supraThresholdSweep-1;
    passiveSweeps = find(cell2mat(allData(iter).analysis.APstats.nSpikesDuringStep(subThreshSweepNos)) == 0); % find all passive sweeps
    if ~isempty(find(diff(passiveSweeps) > 1, 1))
        passiveSweeps = passiveSweeps(1:find(diff(passiveSweeps) > 1, 1));
    end
    Vm_base = mean(allData(iter).recordingData.rawData_mV(1:1000*1000/si_us, passiveSweeps), 1);
    Vm_step = mean(allData(iter).recordingData.rawData_mV(1300*1000/si_us:2000*1000/si_us, passiveSweeps), 1);
    Vm_diff = (Vm_step - Vm_base) / 1000; % /1000 to convert mV to V
    I_step = allData(iter).recordingData.Isteps_pA(passiveSweeps) / 10e12 ; % /10e12 to convert pA to A
    tmp(iter) = length(Vm_diff);
    % calculate linear fit
    RinFitCoeff = polyfit(I_step, Vm_diff, 1);
    xFit = linspace(min(I_step), max(I_step));
    yFit = polyval(RinFitCoeff , xFit);
    % plot results
% % %     fRin = figure; 
% % %         plot(I_step, Vm_diff, '*'); hold on; grid on
% % %         plot(xFit, yFit, 'Color', 'black')
% % %         title('Input resistance sanity check'); legend({'data', 'linear fit'})
% % %         xlabel('current injected [A]'); ylabel('voltage deflection [V]'); 
% % %         set(gca, 'FontSize', 12)
    % store results
    R_in = RinFitCoeff(1)/10e6; % convert Ohm to MOhm
    allData(iter).analysis.Rin_MOhm = R_in;
        
    % calculate membrane time constant 
    % NOTE: filter parameters and exponential fit window & start parameters 
    % might need to be changed depending on the cell type
    % exponential fit time window [ms]
    tRCfit = [1062.6 1112.6] .*1000./20; 
    % lowpass butterworth parameters
    SOS = [1,2,1,1,-1.99376360357872,0.993862527198525;...
           1,2,1,1,-1.98391284186570,0.984011276724947;...
           1,1,0,1,-0.990087914364395,0];
    G = [2.47309049513286e-05; 2.46087148124447e-05; 0.00495604281780237;1];
    % find sweep with -40pA current injection
    minus40pASweepNo = find(allData(iter).recordingData.Isteps_pA == 0)-4; 
    minus40pAData = allData(iter).recordingData.rawData_mV(tRCfit(1):tRCfit(2), minus40pASweepNo);
    minus40pAData = minus40pAData - min(minus40pAData); % substract baseline for better fit
    minus40pAfiltData = filtfilt(SOS, G, minus40pAData); % filter data for better fit
    minus40pAfitTs = (1:length(minus40pAfiltData))/1000*20; % ts for fit [ms]
    % single exponential fit
    fitModel = 'a * (1-exp(-x*b)) + c';
    RCfit = fit(minus40pAfitTs', minus40pAData, fitModel, 'StartPoint', [-10 0.1 10]);
    fitCoeffs = coeffvalues(RCfit);
    allData(iter).analysis.membraneConstant_tau_ms = ...
        1/fitCoeffs(2); %[ms]
% % %     fRCfit = figure; plot(RCfit, minus40pAfitTs, minus40pAfiltData); grid on
% % %         title('membrane constant exponential fit sanity check'); 
% % %         xlabel('time [ms]'); 
% % %         ylabel('voltage deflection [mV]'); 
% % %         set(gca, 'FontSize', 12)

                                              
    % delay to first spike in first suprathreshold sweep
    [allData(iter).analysis.firstSpikeDelay_ms, fDelay] = spikeDelay(allData,...
                                                                     iter,... 
                                                                     1062,...
                                                                     si_us,...
                                                                     true);
                                                       
    % depolarizing hump amplitude in last subthreshold sweep
    [allData(iter).analysis.depolarizingHumpVal_mV,...
     allData(iter).analysis.depolarizingHumpAmpl_mV,... 
     allData(iter).analysis.depolarizingHumpAmpl_percent,...
     allData(iter).analysis.depolarizingHumpTs_ms,...
     fHump] = depolarizingHumpAmplitude(allData,...
                                        iter,...
                                        [1062 2062],...
                                        [1062 1212],...
                                        si_us,...
                                        true);
                                                                       
    % find Vm_rest
    allData(iter).analysis.VmRest_mV = ...
        mean(allData(iter).recordingData.rawData_mV(1:1000*1000/si_us, 1));
    
    % calculate delayed rise due to slow potassium channels
    % extract data
    fitTWinStart = 1500;
    fitTWinEnd = 2000;
	subThreshSweep = allData(iter).recordingData.rawData_mV(:, allData(iter).analysis.supraThresholdSweep-1);
    KvFitTWin = smooth(subThreshSweep(fitTWinStart*1000/20:fitTWinEnd*1000/20));
    KvFitTs = ts(fitTWinStart*1000/20:fitTWinEnd*1000/20);
    % calculate fit
    KvFitCoeff = polyfit(KvFitTs, KvFitTWin, 1);
    xFit = linspace(min(KvFitTs), max(KvFitTs), 1000);
    yFit = polyval(KvFitCoeff , xFit);
    % plot results
% % %     fKvFit = figure;
% % %         plot(ts, subThreshSweep); hold on; grid minor
% % %         ylim([-100 50]); ylabel('V_m [mV]'); xlabel('t [ms]')
% % %         line([fitTWinStart fitTWinStart],[ylim],'LineStyle', '--', 'Color', 'black')
% % %         line([fitTWinEnd fitTWinEnd],[ylim],'LineStyle', '--', 'Color', 'black')
% % %         plot(KvFitTs, KvFitTWin, 'Color', 'cyan')
% % %         plot(xFit, yFit, 'r-', 'LineWidth', 1);
    % store data
    allData(iter).analysis.subThreshKv_mVperms = KvFitCoeff(1);
    
    % find first sweep with FR >= 20 Hz                                   
    Hz20SweepInd = ...
        find(cell2mat(allData(iter).analysis.APstats.nSpikesDuringStep) >= 20, 1);
% % %     f20Hz = figure;
% % %         plot((1:length(allData(iter).recordingData.rawData_mV(:, Hz20SweepInd)))*20/1000,...
% % %               allData(iter).recordingData.rawData_mV(:, Hz20SweepInd),...
% % %               'Color', 'black'); hold on
% % %         plot(allData(iter).analysis.APstats.threshTs{Hz20SweepInd},...
% % %              allData(iter).analysis.APstats.threshVal{Hz20SweepInd},...
% % %              '*', 'Color', 'red')
% % %         plot(allData(iter).analysis.APstats.APmaxTs{Hz20SweepInd},...
% % %             allData(iter).analysis.APstats.APmaxVal{Hz20SweepInd},...
% % %             '*', 'Color', 'green')
% % %         plot(allData(iter).analysis.APstats.AHPminTs{Hz20SweepInd},...
% % %             allData(iter).analysis.APstats.AHPminVal{Hz20SweepInd},...
% % %             '*', 'Color' , 'blue')
% % %         title('first sweep w/ FR >= 20Hz'); 
% % %         xlabel('time [ms]'); 
% % %         ylabel('V_m [mV]'); 
% % %         set(gca, 'FontSize', 12)
    
    % calculate ISI variance & spike frequency adatation
    ISIvarSpikeTs = ...
        allData(iter).analysis.spikeTimes_ms.cell{Hz20SweepInd}((allData(iter).analysis.spikeTimes_ms.cell{Hz20SweepInd} - (Hz20SweepInd-1) * 4000) > 1062 &...
                                                                (allData(iter).analysis.spikeTimes_ms.cell{Hz20SweepInd} - (Hz20SweepInd-1) * 4000) <= 2062);
    ISIvarSpikeDiff = diff(ISIvarSpikeTs);                                                    
    allData(iter).analysis.sweep20Hz.ISIvar = ...
        var(ISIvarSpikeDiff);    
    allData(iter).analysis.sweep20Hz.spikeFreqAdapt = ...
        ISIvarSpikeDiff(9)/ISIvarSpikeDiff(2);        
    
    
    % analyze sag
    negativeSweepIDs = find(allData(iter).recordingData.Isteps_pA < 0); % find all sweeps with pA < 0
    fSag = figure('units','normalized','outerposition',[0 0 0.5 1]); sagPlotCount = 1;
    base75 = NaN;
    base75cellInd = NaN;
    base75SweepInd = NaN;
    base80 = NaN;
    base80cellInd = NaN;
    base80SweepInd = NaN;
    for jter = 1:length(negativeSweepIDs)   
        % extract data
        filtDats = filtfilt(SOS,G,allData(iter).recordingData.rawData_mV(:,jter));
        plot(ts, filtDats, 'LineWidth', 2, 'Color', [0,0,0]+0.66*jter/length(negativeSweepIDs) ); hold on;
        [minVal, minInd] =  min(filtDats(1062*1000/20:1362*1000/20)); % maybe go to 500ms window instead?
        minTs = minInd * 20 / 1000;
        stepBaseline = median(filtDats(1360*1000/20:2060*1000/20));
        sweepBaseline = median(filtDats(1:1060*1000/20));       
        allData(iter).analysis.sagVal_mV{jter}  = minVal;
        allData(iter).analysis.sagAmpl_mv{jter}  = abs(stepBaseline - minVal);
        allData(iter).analysis.sagAmpl_percent{jter}  = (abs(sweepBaseline - minVal)/abs(sweepBaseline - stepBaseline)*100)-100;
        allData(iter).analysis.sagStepBaseline_mV{jter}  = stepBaseline;
        if stepBaseline <= -75
            if isnan(base75)
            	base75 = stepBaseline;
                base75cellInd = iter;
                base75SweepInd = jter;
            elseif stepBaseline > base75
                base75 = stepBaseline;
                base75cellInd = iter;
                base75SweepInd = jter;
            end
        end
        if stepBaseline <= -80
            if isnan(base80)
            	base80 = stepBaseline;
                base80cellInd = iter;
                base80SweepInd = jter;
            elseif stepBaseline > base80
                base80 = stepBaseline;
                base80cellInd = iter;
                base80SweepInd = jter;
            end
        end
        
    end
    % plot cosmetics
    ylim([-120 60]); grid minor; 
    ylabel('V_m [mV]'); xlabel('time [ms]'); title('sag sanity check')
    set(gca, 'FontSize', 18, 'LineWidth', 1) 
    % save data
    allData(iter).analysis.sag.base75 = base75;
    allData(iter).analysis.sag.base75cellInd = base75cellInd;
    allData(iter).analysis.sag.base75SweepInd = base75SweepInd;
    allData(iter).analysis.sag.base80 = base80;
    allData(iter).analysis.sag.base80cellInd = base80cellInd;
    allData(iter).analysis.sag.base80SweepInd = base80SweepInd;   
    % end parameter extraction
    allData(iter).recordingFileName
    iter
%     pause
    close all 
    
    % calculate M-value in first sweep with AP   
    threshSweepApTs = allData(iter).analysis.APstats.threshTs{:, allData(iter).analysis.supraThresholdSweep};
    APind = threshSweepApTs*1000./20;
	thresholdSweep = allData(iter).recordingData.rawData_mV(:, allData(iter).analysis.supraThresholdSweep);
    MVal = zeros(size(threshSweepApTs));
    diffSweep = diff(thresholdSweep).*50;
    for APiter = 1:length(threshSweepApTs)
        diffAPdats = diffSweep(APind(APiter)-5:APind(APiter)+75);
        filtAPdats = smooth(diffAPdats, 1);
        vals = [0, 0, 0];
        locs = [NaN, NaN, NaN];
        [vals(1), locs(1)] = max(filtAPdats);
        try % if there is another peak before absolute max, measure
            [vals(2), locs(2)] = findpeaks(filtAPdats(1:locs(1)),...
                                           'SortStr', 'descend', 'NPeaks', 1); % find first peak = spikelet component in spike/spikelet hybrid    
            [vals(3), tmp] = min(filtAPdats(locs(2):locs(1)));                 % find valley betweenn peaks
            locs(3) = locs(2) + tmp - 1;                                       % calculate location of valley
            deltaDiff = vals(2) - vals(3);                                     % calculate difference = "M value" = measure of shoulder prominence
        catch
            deltaDiff = 0;
        end
        if false
            fControlShoulder = figure;
            set(gcf, 'Position', get(0, 'Screensize'));
            plot(diffAPdats, 'Color', [158, 157, 153]./255); hold on
            plot(filtAPdats, 'Color', 'black')
            plot(locs(1), vals(1), '*', 'Color', 'red')
            plot(locs(2), vals(2), '*', 'Color', 'green')
            plot(locs(3), vals(3), '*', 'Color', 'blue')
            legend('1st derivative',...
                   'smoothed data',...
                   'primary peak',...
                   'secondary peak',...
                   'min. \DeltaV/\Deltat')
        end 
        MVal(APiter) = deltaDiff;
    end
    allData(iter).analysis.supraThresholdMVal = MVal;
 
    
end


%% ------------------------ inter-step firing -------------------------- %%
for iter = 1:nCells
    % preallocate collector variable
    interStepFR = NaN(length(allData(iter).analysis.spikeTimes_ms.cell)-1,1); 
    
    % loop sweeps
	for jter = 1:length(allData(iter).analysis.spikeTimes_ms.cell)-1
        % calculate time window
        tStart = (jter-1) * totalDuration + 2065; % actually, 2062 - 2065 used as safety margin
        tStop = jter * totalDuration + 1059; % actually, 1062 - 1059 used as safety margin 
        % find number of spikes in time window and calculate FR
        interStepSpkCount = sum(allData(iter).analysis.spikeTimes_ms.mat > tStart &...
                                allData(iter).analysis.spikeTimes_ms.mat < tStop);
       
        interStepFR(jter) = interStepSpkCount / (((4000-2070)+1060)/1000); % /t[sec]; written explicitly for understandability
    end
    % collect data and add to allData
    allData(iter).analysis.interStepFR = interStepFR;
end


% ---------------- calculate barrage firing duration and FR --------------- 
% loop cells
maxInterStepFR = [];
for iter = 1:length(allData)
    maxInterStepFR = max(allData(iter).analysis.interStepFR);
    barrageFR = [];
    barrageDuration = [];
    nSpikes = [];
    % exlude cells without significant inter-step FR
    if maxInterStepFR == 0
        allData(iter).analysis.maxBarrageDuration_ms = nan;
        allData(iter).analysis.avgBarrageDuration_ms = nan;
        allData(iter).analysis.firstBarrageDuration_ms = nan;
        allData(iter).analysis.barrageDuration_ms = zeros(size(allData(iter).analysis.interStepFR));
        allData(iter).analysis.maxBarrageFR_Hz = nan;
        allData(iter).analysis.avgBarrageFR_Hz = nan;
        allData(iter).analysis.firstBarrageFR_Hz = nan;
        allData(iter).analysis.barrageFR_Hz = zeros(size(allData(iter).analysis.interStepFR));
        allData(iter).analysis.nBarrageSpikes = nan;
        allData(iter).analysis.nBarrageSteps = nan;
    else
        % loop inter-step intervals
        for jter = 1:length(allData(iter).analysis.spikeTimes_ms.cell)-1
            tStart = (jter-1) * totalDuration + 2062;
            tStop = jter * totalDuration + 1062;
            % extract spike times in current inter-step interval
            interStepAPs = allData(iter).analysis.spikeTimes_ms.mat(...
                                allData(iter).analysis.spikeTimes_ms.mat > tStart &...
                                allData(iter).analysis.spikeTimes_ms.mat < tStop);
            % if no spikes in current interval, assign zeros
            if isempty(interStepAPs)
                barrageDuration(jter) = 0;
                barrageFR(jter) = 0;
                nSpikes(jter) = 0;
            % else calculate those stats
            else
                % get number of spikes in current interval
                nSpikes(jter) = length(interStepAPs);
                % measure time from end of step to last inter-step AP [ms]
                barrageDuration(jter) = max(interStepAPs) - ((jter-1)*4000 + 2062); 
                 % calculate amount of APs over that time period [Hz]
                barrageFR(jter) = nSpikes(jter)/(barrageDuration(jter)/1000);
            end
        end
        
        % find maximum barrage FR & duration
        maxBarrageFR = max(barrageFR);
        maxBarrageDuration = max(barrageDuration);

        % find all intervals with spikes  
        boolInd = find(nSpikes >= 1);
        % find first intervals and extract stats
        firstBarrageFR = barrageFR(boolInd(1));
        firstBarrageDuration = barrageDuration(boolInd(1));
        % calculate average stats over all intervals (with more than 5 spikes)
        avgBarrageDuration = mean( barrageDuration(boolInd));
        avgBarrageFR = mean( barrageFR(boolInd));
        % number of steps wih pF
        pf_steps = nSpikes > 0; 
        n_pF_steps= sum(pf_steps);
        
        % add data to allData
        allData(iter).analysis.maxBarrageDuration_ms = maxBarrageDuration;
        allData(iter).analysis.avgBarrageDuration_ms = avgBarrageDuration;
        allData(iter).analysis.firstBarrageDuration_ms = firstBarrageDuration;
        allData(iter).analysis.barrageDuration_ms = barrageDuration;
        allData(iter).analysis.maxBarrageFR_Hz = maxBarrageFR;
        allData(iter).analysis.avgBarrageFR_Hz = avgBarrageFR;
        allData(iter).analysis.firstBarrageFR_Hz = firstBarrageFR;
        allData(iter).analysis.barrageFR_Hz = barrageFR;
        allData(iter).analysis.nBarrageSpikes = nSpikes;
        supraSweep = allData(iter).analysis.supraThresholdSweep;
        allData(iter).analysis.nEctopicSupraTotal = sum(nSpikes(supraSweep:end));
        allData(iter).analysis.nEctopicSupraMean = mean(nSpikes(supraSweep:end));
        allData(iter).analysis.nEctopicSupraMax = max(nSpikes(supraSweep:end));
        allData(iter).analysis.nBarrageSteps = n_pF_steps;
    end
end
    


%% -------------------------- RESULTS TABLE ---------------------------- %%
% ----------------------------------------------------------------------- %
% intrinsic properties
firstSpikeDelay_ms = NaN(size(allData))';
AP_amplitude_mV = NaN(size(allData))';
AP_threshold_mV = NaN(size(allData))';
AP_halfWidth_ms = NaN(size(allData))';
AP_maxSpikeSlope_mVpms = NaN(size(allData))';
peakToTrough_ms = NaN(size(allData))';
AHP_amplitude_mV = NaN(size(allData))';
AHP_halfWidth_ms = NaN(size(allData))';
ISI_adapt = NaN(size(allData))';
ISI_var = NaN(size(allData))';
Vm_rest = NaN(size(allData))';
Rinput_MOhm = NaN(size(allData))';
tau_ms = NaN(size(allData))';
depolHump_mV = NaN(size(allData))';
depolHump_percent = NaN(size(allData))';
barrageFiring_maxFR_Hz = NaN(size(allData))';
barrageFiring_maxT_ms = NaN(size(allData))';
sagAmpl_percent80 = NaN(size(allData))';
sagAmpl_percent75 = NaN(size(allData))';
sagBaseMin_mV = NaN(size(allData))';
subThreshKv_mVperms = NaN(size(allData))';
reboundSpikes = NaN(size(allData))';
maxFR = NaN(size(allData))';

% meta info
cellNo = NaN(size(allData))';
cellType = NaN(size(allData))';
suitableBFanalysis = NaN(size(allData))';
isBF = NaN(size(allData))';
lineID = NaN(size(allData))';
mouseID = NaN(size(allData))';

% barrage data
maxBarrageDuration_ms  = NaN(size(allData))';
avgBarrageDuration_ms = NaN(size(allData))';
firstBarrageDuration_ms = NaN(size(allData))';
barrageDuration_ms = NaN(size(allData))';
maxBarrageFR_Hz = NaN(size(allData))';
avgBarrageFR_Hz = NaN(size(allData))';
firstBarrageFR_Hz = NaN(size(allData))';
barrageFR_Hz = NaN(size(allData))';
nBarrageSpikes = NaN(size(allData))';
nEctopicSupraTotal = NaN(size(allData))';
nEctopicSupraMean = NaN(size(allData))';
nEctopicSupraMax = NaN(size(allData))';
nBarrageSteps  = NaN(size(allData))';

% experiment
experimentID = {};

% extract
for iter = 1:length(allData)
    
    % cellNo for later identification
    cellNo(iter) = iter;
    experimentID{iter} = convertCharsToStrings(allData(iter).recordingFileName);
    
    % extract previously calculated intrinsic properties
    Hz20SweepInd = find(cell2mat(allData(iter).analysis.APstats.nSpikesDuringStep) >= 20, 1);
    sweepInds = characterizationPlot(allData, iter);
    supraThresholdSweep = sweepInds(3);  
    firstSpikeDelay_ms(iter) = allData(iter).analysis.firstSpikeDelay_ms;
    AP_amplitude_mV(iter) = allData(iter).analysis.APstats.APamplitude{supraThresholdSweep}(1);
    [iter allData(iter).analysis.APstats.APamplitude{supraThresholdSweep}(1)];
    AP_threshold_mV(iter) = allData(iter).analysis.APstats.threshVal{supraThresholdSweep}(1);
    AP_halfWidth_ms(iter) = allData(iter).analysis.APstats.APhalfWidth{supraThresholdSweep}(1);
    AP_maxSpikeSlope_mVpms(iter) = allData(iter).analysis.APstats.APmaxSlopeVal{supraThresholdSweep}(1);
    peakToTrough_ms(iter) = allData(iter).analysis.APstats.peakToTrough{supraThresholdSweep}(1);
    AHP_amplitude_mV(iter) = allData(iter).analysis.APstats.AHPamplitude{supraThresholdSweep}(1);
    AHP_halfWidth_ms(iter) = allData(iter).analysis.APstats.AHPhalfWidth{supraThresholdSweep}(1);
    ISI_adapt(iter) = allData(iter).analysis.sweep20Hz.spikeFreqAdapt;
    ISI_var(iter) = allData(iter).analysis.sweep20Hz.ISIvar;
    Vm_rest(iter) = allData(iter).analysis.VmRest_mV;
    Rinput_MOhm(iter) = allData(iter).analysis.Rin_MOhm;
    tau_ms(iter) = allData(iter).analysis.membraneConstant_tau_ms;
    depolHump_mV(iter) = allData(iter).analysis.depolarizingHumpAmpl_mV;
    depolHump_percent(iter) = allData(iter).analysis.depolarizingHumpAmpl_percent;
    subThreshKv_mVperms(iter) = allData(iter).analysis.subThreshKv_mVperms;
    reboundSpikes(iter) = allData(iter).analysis.APstats.nSpikesTotal{1};
    
    % extract sag data
    if ~isnan(allData(iter).analysis.sag.base75SweepInd)
        sagAmpl_percent75(iter) = ...
            allData(iter).analysis.sagAmpl_percent{allData(iter).analysis.sag.base75SweepInd};   
    end 
    if ~isnan(allData(iter).analysis.sag.base80SweepInd)
        sagAmpl_percent80(iter) = ...
            allData(iter).analysis.sagAmpl_percent{allData(iter).analysis.sag.base80SweepInd}; 
    end
    sagBaseMin_mV(iter) = min([allData(iter).analysis.sagStepBaseline_mV{:}]);
   
    
    % extract BF data
    if allData(iter).analysis.suitableForBFanalysis
        suitableBFanalysis(iter) = 1;
        if ~isnan(allData(iter).analysis.maxBarrageFR_Hz)
            if allData(iter).analysis.nEctopicSupraTotal > 5; isBF(iter) = 1; % cutoff to exclude cells which fire once in a while without it being BF, see histogram of ectopic spikes
            else; isBF(iter) = 0; end
            maxBarrageDuration_ms(iter)  = allData(iter).analysis.maxBarrageDuration_ms;
            avgBarrageDuration_ms(iter) = allData(iter).analysis.avgBarrageDuration_ms;
            firstBarrageDuration_ms(iter) = allData(iter).analysis.firstBarrageDuration_ms;
            maxBarrageFR_Hz(iter) = allData(iter).analysis.maxBarrageFR_Hz;
            avgBarrageFR_Hz(iter) = allData(iter).analysis.avgBarrageFR_Hz;
            firstBarrageFR_Hz(iter) = allData(iter).analysis.firstBarrageFR_Hz;
            nEctopicSupraTotal(iter) = allData(iter).analysis.nEctopicSupraTotal;
            nEctopicSupraMean(iter) = allData(iter).analysis.nEctopicSupraMean;
            nEctopicSupraMax(iter) = allData(iter).analysis.nEctopicSupraMax;
            nBarrageSteps(iter) = allData(iter).analysis.nBarrageSteps;
        else
            isBF(iter) = 0;
            maxBarrageDuration_ms(iter)  = 0;
            avgBarrageDuration_ms(iter) = 0;
            firstBarrageDuration_ms(iter) = 0;
            maxBarrageFR_Hz(iter) = 0;
            avgBarrageFR_Hz(iter) = 0;
            firstBarrageFR_Hz(iter) = 0;
            nEctopicSupraTotal(iter) = 0;
            nEctopicSupraMean(iter) = 0;
            nEctopicSupraMax(iter) = 0;
            nBarrageSteps(iter) = 0;
        end
    else
        suitableBFanalysis(iter) = 0;
    end
    
    % extract maxFR: only valid if max FR is actually reached, which is only known if depoli block, which is criterion for BF analysis also
    if suitableBFanalysis(iter) == 1       
        maxFR(iter) = max(allData(iter).analysis.FR_Hz_step);
    end
    
    % check cell type
    if allData(iter).cellType == "NDNF+"; cellType(iter) = 1;
    elseif allData(iter).cellType == "NDNF+/NPY+"; cellType(iter) = 2;
    elseif allData(iter).cellType == "L23_FS_IN"; cellType(iter) = 3;
    end 

end


% ----------------- INTRINSIC PROPERTIES FOR CLUSTERING -------------------
allNames = {'cellNo',...
            'cellType',...
            'useBFanalysis',...
            'isBF',...
            'firstSpikeDelay_ms',...
            'AP_amplitude_mV',...
            'AP_threshold_mV',...
            'AP_halfWidth_ms',...
            'AP_maxSpikeSlope_mVpms',...
            'AHP_amplitude_mV',...
            'AHP_halfWidth_ms',...
            'ISI_adapt',...
            'ISI_var',...
            'Vm_rest',...
            'Rinput_MOhm',...
            'tau_ms',...
            'depolHump_percent',...
            'sagAmpl_75_precent',...
            'maxFR_Hz_evoked',...
            'maxBarrageDuration_ms',...
            'avgBarrageDuration_ms',...
            'firstBarrageDuration_ms',...
            'maxBarrageFR_Hz',...
            'avgBarrageFR_Hz',...
            'firstBarrageFR_Hz',...
            'nEctopicSupraTotal',...
            'nEctopicSupraMean',...
            'nEctopicSupraMax',...
            'nBarrageSteps'};

metaInfo = [cellNo,...
            cellType,...
            suitableBFanalysis,...
            isBF];
        
intrProp = [firstSpikeDelay_ms,...
            AP_amplitude_mV,...
            AP_threshold_mV,...
            AP_halfWidth_ms,...
            AP_maxSpikeSlope_mVpms,...
            AHP_amplitude_mV,...
            AHP_halfWidth_ms,...
            ISI_adapt,...
            ISI_var,...
            Vm_rest,...
            Rinput_MOhm,...
            tau_ms,...
            depolHump_percent,...
            sagAmpl_percent75,...
            maxFR];
           
BFparams = [maxBarrageDuration_ms,...
            avgBarrageDuration_ms,...
            firstBarrageDuration_ms,...
            maxBarrageFR_Hz,...
            avgBarrageFR_Hz,...
            firstBarrageFR_Hz,...
            nEctopicSupraTotal,...
            nEctopicSupraMean,...
            nEctopicSupraMax,...
            nBarrageSteps];

resultTable_currentSteps = array2table([metaInfo, intrProp, BFparams], "VariableNames", allNames);


%% ------------------ INTRINSIC PROPERTIES FOR CLUSTERING ------------------
colNames = {'firstSpikeDelay_ms',...
            'AP_amplitude_mV',...
            'AP_threshold_mV',...
            'AP_halfWidth_ms',...
            'AP_maxSpikeSlope_mVpms',...
            'AHP_amplitude_mV',...
            'AHP_halfWidth_ms',...
            'ISI_adapt',...
            'ISI_var',...
            'Vm_rest',...
            'Rinput_MOhm',...
            'tau_ms',...
            'depolHump_percent',...
            'sagAmpl_75_precent'};
        
collectData = [firstSpikeDelay_ms....
               AP_amplitude_mV,...
               AP_threshold_mV,...
               AP_halfWidth_ms,...
               AP_maxSpikeSlope_mVpms,...
               AHP_amplitude_mV,...
               AHP_halfWidth_ms,...
               ISI_adapt,...
               ISI_var,...
               Vm_rest,...
               Rinput_MOhm,...
               tau_ms,...
               depolHump_percent,...
               sagAmpl_percent75]; % use this if need the actual values for all cells  
           
% clean Data
cleanDataInd = ~isnan(sagAmpl_percent75);% &...
               %Vm_rest < -60 &...
               %AP_threshold_mV < -30;

% all data
barrageFiring_maxT_ms_all = barrageFiring_maxT_ms;
collectData_all = collectData;
cellType_all = cellType;

% cleaned data for intrinsic property analysis
collectData = collectData(cleanDataInd, :);
barrageFiring_maxT_ms = barrageFiring_maxT_ms(cleanDataInd);
cellType = cellType(cleanDataInd);

% clustering data
normCollectData = normalize(collectData, 1, 'range'); % use this if all cells should be included in method, e.g. PCA

%% clean up
% clearvars -except allData allNames resultTable_currentSteps colNames collectData cleanDataInd collectData_all normCollectData experimentID
