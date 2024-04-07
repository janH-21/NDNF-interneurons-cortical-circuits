%% ---------------------- user editable block -------------------------- %%
% working directory
cd '/main/'
close all; clc; clear all;


%% ------------------------------ main --------------------------------- %%
% ------------------------ choose and load file ------------------------- % 
[file, path] = uigetfile({'*.abf'});                                       % NOTE:
[d,si,h]=abfload(fullfile(path,file));                                     % d - (datapoint, channel, sweep)
                                                                           % si - sampling interval in µs                                                                   % h - meta
% format data
nSweeps = size(d,3);
d = d(:,:,1:10);                                                           % only consider first 10 sweeps
dats = reshape(d(:,1,:),[],1);
sweep = reshape(d(:,2,:),[],1);
t = (1:length(dats)).*si/1000;

% ----------------------------- analysis -------------------------------- %
%  detect spikes
[APts, APtsVals, intT, intD, diffInt] = ...
    detectSpikes(dats, si, 10, 10, 3, 15, 100, false);

% measure AP stats
MVal = NaN(size(APts));                                                    
APampl =  NaN(size(APts));                                                 
APmax = NaN(size(APts));                                                   
diffDats = diff(dats).*1000./si;                                          
spkMaxSlope =  NaN(size(APts));
for iter = 1:length(APts)                                                  % iterate AP time stamps
    % M-value
    diffSpike = diffDats(APts(iter)*1000/si - 0.2*1000/si:...
                         APts(iter)*1000/si + 2*1000/si);                  % extract 1st deriv. data around spike ts +- window
    [~, ~, MVal(iter)] = measureSpikeShoulder(diffSpike, 1, false);        
    spkMaxSlope(iter) = max(diffSpike);
    % AP amplitude
    APdats = dats(APts(iter)*1000/si:...                                   % as above, but with raw data
                  APts(iter)*1000/si + 1.5*1000/si);
    APmax(iter) = max(APdats);
    APampl(iter) = max(APdats) - APtsVals(iter);
end

% instantaneous firing rate
ifreq = 1000./diff(APts);

% find step times
[stepTimes, interStepTimes] = findStepTimes(sweep, si);

% find FRs during and between steps
[stepFRs, interStepFRs] = findStepFRs(stepTimes, interStepTimes, APts);

% analyse Vm between steps
interStepDuration = interStepTimes(:,2)-interStepTimes(:,1);
interStepBool = interStepDuration < 251 & interStepDuration > 249;
interSweepBool = interStepDuration < 751 & interStepDuration > 749;
Vm = NaN(size(interStepBool));
Vmin = NaN(size(interStepBool));
for iter = 2:size(interStepTimes,1)-1                                      % skip first before stimulation starts
    % inter step of inter sweep?
    % convert times to samples
    interStepSamples = interStepTimes(iter, :)*1000/si;
    % minimum in first 3rd of interval
    Vmin(iter) = min(dats(interStepSamples(1):...
                          interStepSamples(1)+floor(diff(interStepSamples)/3)));
    % baseline in last 3rd of interval
    Vm(iter) = mean(dats(interStepSamples(1)+floor(diff(interStepSamples)*2/3):...
                         interStepSamples(1)+floor(diff(interStepSamples)*0.95)));
end
AHP = Vm - Vmin;
interStep_Vm = Vm; interStep_Vm(~interStepBool) = NaN;
interStep_AHP = AHP; interStep_AHP(~interStepBool) = NaN;
interSweep_Vm = Vm; interSweep_Vm(~interSweepBool) = NaN;
interSweep_AHP = AHP; interSweep_AHP(~interSweepBool) = NaN;


% plot results 
fOverview = plot250msProtocol(t, dats, APts, MVal, APtsVals);
fStats = plotStats(MVal, APtsVals, APampl, APmax, 100, "pdf");
figure; plot(stepFRs, "o-"); hold on; plot(interStepFRs, "o-"); 
file
[nSweeps si stepFRs(1)]

%% save data
% If there is no file yet, run:
    % allData_250ms_2Hz_7APs = struct();
    % cellCount = 1;
% To continue an exsiting file, run:
    load('allData_250ms_2Hz_7APs_blockers.mat'); 
    cellCount = length(allData_250ms_2Hz_7APs)+1; 

% meta info
allData_250ms_2Hz_7APs(cellCount).meta.file = file;
allData_250ms_2Hz_7APs(cellCount).meta.si = si;
allData_250ms_2Hz_7APs(cellCount).meta.nSweeps = size(d,3);
% action potential stats
allData_250ms_2Hz_7APs(cellCount).analysis.APs.APts_ms = APts;
allData_250ms_2Hz_7APs(cellCount).analysis.APs.APtsVals_mV = APtsVals;
allData_250ms_2Hz_7APs(cellCount).analysis.APs.MVal_mVpms = MVal;
allData_250ms_2Hz_7APs(cellCount).analysis.APs.ifreq_Hz= ifreq;
allData_250ms_2Hz_7APs(cellCount).analysis.APs.APampl_mV = APampl;
allData_250ms_2Hz_7APs(cellCount).analysis.APs.APmax_mV = APmax;
% other stats
allData_250ms_2Hz_7APs(cellCount).analysis.step.stepTimes_ms = stepTimes;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interStepTimes_ms = interStepTimes;
allData_250ms_2Hz_7APs(cellCount).analysis.step.stepFRs_Hz = stepFRs;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interStepFRs_Hz = interStepFRs;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interSweepVm_mV = interSweep_Vm;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interSweepAHP_mV = interSweep_AHP;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interStepVm_mV = interStep_Vm;
allData_250ms_2Hz_7APs(cellCount).analysis.step.interStepAHP_mV = interStep_AHP;
% output
save('allData_250ms_2Hz_7APs.mat', 'allData_250ms_2Hz_7APs')


%% subroutines
function [APts, APtsVals, intT, intD, diffInt] = detectSpikes(rawD, si, vThres, diffThres, spikeWindow, slopeCrit, interpFact, controlPlot)
    % ---------------------------INPUTS------------------------------------
    % raw              - raw data
    % si               - sampling interval
    % diffThres        - dV/dt threshold in mV/ms for event detection
    % vThres           - amplitude threshold in mV for counting event as
    %                    spike
    % spikeWindow      - time in ms to consider after threshold crossing
    % slopeCrit        - average minimum Vm rise after detected onset crit-
    %                    to count onset as true onset, over 0.5ms, in mV/ms
    % interpFact       - factor by which to interpolate, 1: no
    %                    interpolation, 10: 10x interpolation etc.; change
    %                    according to sampling interval; WARNING: might
    %                    affect spike detection severly if out of range
    % controlPlot      - boolean, plot sanity check plots
    %
    % ---------------------------OUTPUTS-----------------------------------
    % APts             - action potential time stamps in interpolated data
    % intD             - interpolated data
    
    % prepare data
    t = (1:size(rawD)).*si/1000;                                           % time in ms
    intT = (1:(1/interpFact):size(rawD)).*si/1000;                         % time stamps for interpolation
    intD = interp1(t, rawD, intT);                                         % interpolate data for more precise measurements
    diffRaw = diff(rawD)./si*1000;                                         % dV/dt in mV/ms
    diffInt = interp1(t(2:end), diffRaw, intT(2:end));                     % interpolate dV/dt
    
    % detect spike onsets
    [pks, lcs] = findpeaks(-abs(diffThres-diffInt), 'MinPeakHeight', -1.5);  % find dV/dt threshold crossings
    lcBoolInds = ones(size(lcs));                                          % preallocate
    for iter = 1:length(lcs)                                               % check whether threshold crossings are really spike onsets
        if ~(iter == length(lcs))                                          % only do comparisons to subsequent threshold crossing if there is a subsequent one
            if (intT(lcs(iter+1)) - intT(lcs(iter))) > spikeWindow         % spike window criterion
                 lcBoolInds(iter) = 0; end
            if (intD(lcs(iter+1)) - intD(lcs(iter))) < vThres              % spike amplitude criterion
                lcBoolInds(iter) = 0; end
        end

        siFact = si/100;
        if lcs(iter)+100/siFact > length(intD)
            lcBoolInds(iter) = 0;
        else
            if intD(lcs(iter)) > mean(intD(lcs(iter):lcs(iter)+100/siFact))    % directionality of spike slope criterion -> this also deals with last threshold crossing in trace in most cases
                lcBoolInds(iter) = 0; end
            if mean(diffInt(lcs(iter):lcs(iter)+200/siFact)) < slopeCrit       % average voltage slope criterion -> sorts out events with slope to low for spike
                lcBoolInds(iter) = 0; end
        end
    end
    lcsTrue = lcs(logical(lcBoolInds));                                    % exclude ts which didn't pass
    pksTrue = pks(logical(lcBoolInds));                                    % same for APs
    APts = intT(lcsTrue+1);                                                % output ts
    APtsVals = intD(lcsTrue);                                              % output values
    
    % control figure
    if controlPlot
        col1 = [4, 160, 191]./255;
        col2 = [235, 173, 5]./255;
        colThres = [46, 247, 10]./255;
        colAPts = 'red';
        fControl = figure;
        subplot(3,1,1)
            plot(t, rawD, 'o-', 'MarkerSize', 1, 'Color', col1); hold on
            plot(intT, intD, 'o', 'MarkerSize', 1, 'Color', col2)
            plot(intT(lcs+1), intD(lcs), 'o', 'Color', colThres)
            plot(intT(lcsTrue+1), intD(lcsTrue), '*', 'Color', colAPts)
            grid on
            legend('raw data',...
                   'interpolated data',...
                   'points of threshold crossing',...
                   'true action potential onsets',...
                   'Location', 'southeast')
            title('Sweep Data');
            ylabel('V_m [mV]'); xlabel('Time [ms]')
        subplot(3,1,2)
            plot(t(2:end), diffRaw, 'o-', 'MarkerSize', 1, 'Color', col1); hold on
            plot(intT(2:end), diffInt, 'o', 'MarkerSize', 1, 'Color', col2)
            line([xlim], [diffThres diffThres], 'LineStyle', ':', 'Color', 'red')
            plot(intT(lcs+1), diffInt(lcs), 'o', 'Color', colThres)
            plot(intT(lcsTrue+1), diffInt(lcsTrue), 'x', 'Color', colAPts)
            grid on
            legend('\DeltaV/\Deltat',...
                   'interpolated \DeltaV/\Deltat',...
                   'dV/dt threshold',...
                   'points of threshold crossing',...
                   'true action potential onsets',...
                   'Location', 'southeast')
            title('1st derivative'); 
            ylabel('\DeltaV/\Deltat'); xlabel('Time [ms]')
        subplot(3,1,3)
            plot(intT(2:end), -abs(diffThres-diffInt), 'o-', 'MarkerSize', 1, 'Color', col1); hold on
            plot(intT(lcs+1), pks, 'o', 'Color', colThres)
            plot(intT(lcsTrue+1), pksTrue, 'x', 'Color', colAPts)
            grid on
            legend('-abs(\DeltaV/\Deltat threshold - interpolated vV/\Deltat)',...
                   'points of threshold crossing',...
                   'true action potential onsets',...
                   'Location', 'southeast')
            title('-|1st derivative threshold - 1st derivative|'); 
            ylabel('\DeltaV/\Deltat'); xlabel('Time [ms]')
    end
end


function [vals, locs, deltaDiff] = measureSpikeShoulder(diffAPdats, smoothFactor, controlPlot)
    vals = [0, 0, 0];                                                      % preallocate
    locs = [NaN, NaN, NaN];                                                % preallocate
    filtAPdats = smooth(diffAPdats, smoothFactor);                         % filter data
    [vals(1), locs(1)] = max(filtAPdats);                                  % find absolute maximum
    try % if there is another peak before absolute max, measure
        [vals(2), locs(2)] = findpeaks(filtAPdats(1:locs(1)),...
                                       'SortStr', 'descend', 'NPeaks', 1); % find first peak = spikelet component in spike/spikelet hybrid    
        [vals(3), tmp] = min(filtAPdats(locs(2):locs(1)));                 % find valley betweenn peaks
        locs(3) = locs(2) + tmp - 1;                                       % calculate location of valley
        deltaDiff = vals(2) - vals(3);                                     % calculate difference = "M value" = measure of shoulder prominence
    catch % else, there is no deltaDiff
        deltaDiff = 0;
    end
    
    if controlPlot
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
end


function f = plot250msProtocol(intT, intD, APts, MVal, APtsVals)
    f = figure;
    set(gcf, 'Position', get(0, 'Screensize'));
    
    subplot(3,1,1)
    plot(intT, intD, 'Color', 'black', 'LineWidth', 0.1)
    title('raw data'); ylabel('V_m [mV]'); xlabel('Time [ms]')
    set(gca, 'FontSize', 14); ylim([-85 60]); grid on
    
    subplot(3,1,2)
    plot(APts, MVal, '.')
    title('shoulder prominendce'); ylabel('max. deceleration [mV/ms]'); xlabel('time [ms]')
    set(gca, 'FontSize', 14); ylim([0 150]); grid on
    
    subplot(3,1,3)
    plot(APts, APtsVals, '.'); hold on
    title('AP threshold'); ylabel('\DeltaV/\Deltat > 10mV/ms [mV]'); xlabel('time [ms]')
    set(gca, 'FontSize', 14); ylim([-80 0]); grid on
end


function f = plotStats(MVal, APtsVals, APampl, APmax, nBins, histNorm)
    fSize = get(0, 'Screensize');
    f = figure('Position', [1 1 fSize(4) fSize(4)]);
    % M-Value
    subplot(5,5,1); 
        histogram(MVal, "Normalization", histNorm, "NumBins", nBins, "Orientation", "horizontal")
        set(gca, "XDir", "reverse", "FontSize", 14); ylabel("max. deceleration [mV/ms]"); ylim([0 max(MVal)*1.1]); grid on
        subplot(5,5,22); 
        histogram(MVal, "Normalization", histNorm, "NumBins", nBins, "Orientation", "vertical")
        set(gca, "YDir", "reverse", "FontSize", 14); xlabel("max. deceleration [mV/ms]"); xlim([0 max(MVal)*1.1]); grid on
    % AP threshold
    subplot(5,5,6); 
        histogram(APtsVals, "Normalization", histNorm, "NumBins", nBins, "Orientation", "horizontal")
        set(gca, "XDir", "reverse", "FontSize", 14); ylabel("AP threshold [mV]"); ylim([min(APtsVals)*1.1 max(APtsVals)*0.9]); grid on
        subplot(5,5,23); 
        histogram(APtsVals, "Normalization", histNorm, "NumBins", nBins, "Orientation", "vertical")
        set(gca, "YDir", "reverse", "FontSize", 14); xlabel("AP threshold [mV]"); xlim([min(APtsVals)*1.1 max(APtsVals)*0.9]); grid on
    % AP amplitude
    subplot(5,5,11); 
        histogram(APampl, "Normalization", histNorm, "NumBins", nBins, "Orientation", "horizontal")
        set(gca, "XDir", "reverse", "FontSize", 14); ylabel("AP amplitude [mV]"); ylim([min(APampl)*0.9 max(APampl)*1.1]); grid on
        subplot(5,5,24); 
        histogram(APampl, "Normalization", histNorm, "NumBins", nBins, "Orientation", "vertical")
        set(gca, "YDir", "reverse", "FontSize", 14); xlabel("AP amplitude [mV]"); xlim([min(APampl)*0.9 max(APampl)*1.1]); grid on
    % AP absolute max. value
    subplot(5,5,16); 
        histogram(APmax, "Normalization", histNorm, "NumBins", nBins, "Orientation", "horizontal")
        set(gca, "XDir", "reverse", "FontSize", 14); ylabel("AP abs(maxVal) [mV]"); ylim([min(APmax)-5 max(APmax)*1.1]); grid on
        subplot(5,5,25); 
        histogram(APmax, "Normalization", histNorm, "NumBins", nBins, "Orientation", "vertical")
        set(gca, "YDir", "reverse", "FontSize", 14); xlabel("AP abs(maxVal) [mV]"); xlim([min(APmax)-5 max(APmax)*1.1]); grid on  
    % column 1
    subplot(5,5,2)
        scatter(MVal, MVal, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([0 max(MVal)*1.1]); ylim([0 max(MVal)*1.1]); grid on
    subplot(5,5,7)
        scatter(MVal, APtsVals, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([0 max(MVal)*1.1]); ylim([min(APtsVals)*1.1 max(APtsVals)*0.9]); grid on
    subplot(5,5,12)
        scatter(MVal, APampl, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([0 max(MVal)*1.1]); ylim([min(APampl)*0.9 max(APampl)*1.1]); grid on
    subplot(5,5,17)
        scatter(MVal, APmax, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([0 max(MVal)*1.1]); ylim([min(APmax)-5 max(APmax)*1.1]); grid on
    % column 2
    subplot(5,5,8)
        scatter(APtsVals, APtsVals, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APtsVals)*1.1 max(APtsVals)*0.9]); ylim([min(APtsVals)*1.1 max(APtsVals)*0.9]); grid on
    subplot(5,5,13)
        scatter(APtsVals, APampl, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APtsVals)*1.1 max(APtsVals)*0.9]); ylim([min(APampl)*0.9 max(APampl)*1.1]); grid on
    subplot(5,5,18)
        scatter(APtsVals, APmax, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APtsVals)*1.1 max(APtsVals)*0.9]); ylim([min(APmax)-5 max(APmax)*1.1]); grid on
    % column 3
    subplot(5,5,14)
        scatter(APampl, APampl, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APampl)*0.9 max(APampl)*1.1]); ylim([min(APampl)*0.9 max(APampl)*1.1]); grid on
    subplot(5,5,19)
        scatter(APampl, APmax, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APampl)*0.9 max(APampl)*1.1]); ylim([min(APmax)-5 max(APmax)*1.1]); grid on
    % column 4
    subplot(5,5,20)
        scatter(APmax, APmax, '.'); set(gca, "FontSize", 14, "XTickLabel", [], "YTickLabel", [])
        xlim([min(APmax)-5 max(APmax)*1.1]); ylim([min(APmax)-5 max(APmax)*1.1]); grid on
end


function [stepTimes, interStepTimes] = findStepTimes(sweep, si)
    % find step onset via first derivative
    [~, onsets] = findpeaks(diff(sweep), "MinPeakDistance", 10, "MinPeakHeight", 0.1);
    [~, offsets] = findpeaks(-diff(sweep), "MinPeakDistance", 10, "MinPeakHeight", 0.1);
    % samples to time stamps
    onsets = onsets.*si./1000;
    offsets = offsets.*si./1000;
    % collect & return results
    stepTimes = [onsets, offsets];
    interStepTimes = [[0; offsets], [onsets; length(sweep)*si/1000]];
end

function [stepFRs, interStepFRs] = findStepFRs(stepTimes, interStepTimes, APts)
    % preallocate
    stepFRs = zeros(length(stepTimes),1);
    interStepFRs = zeros(length(interStepTimes),1);
    % step firing rates
    for iter = 1:length(stepTimes)
        spikeCount = sum(APts > stepTimes(iter,1) & APts < stepTimes(iter,2));
        stepDuration = (stepTimes(iter,2) - stepTimes(iter,1))/1000; % [secs]
        stepFRs(iter) = spikeCount/stepDuration; % [Hz]
    end
    % inter step firing rates
    for iter = 1:length(interStepTimes)
        spikeCount = sum(APts > interStepTimes(iter,1) & APts < interStepTimes(iter,2));
        stepDuration = (interStepTimes(iter,2) - interStepTimes(iter,1))/1000; % [secs]
        interStepFRs(iter) = spikeCount/stepDuration; % [Hz]
    end
end