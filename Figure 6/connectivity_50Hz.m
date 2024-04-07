%% Stimulation with 45% LED power over 4 seconds at 50Hz with 5ms pulse width
clear all; close all; clc; 
cd '/main/'
stim_on_ind = 135939;
stim_on_ts = stim_on_ind * 20 / 1000;
stim_off_ind = 335189;
stim_off_ts = stim_off_ind * 20 / 1000;

%% Calibration curve
path = './data/calibration_50Hz_NDNF_ChR2/';
lDir = dir(path);
AP_max_ts = {};
AP_max_val = {};
AP_max_ts_bin = [];
counter = 1;
for iter = 3:length(lDir)
    % load data
    file = lDir(iter).name;
    [d,si,h]=abfload(fullfile(path,file));
    dats = d(:,1,1);
    stim = d(:,2,1);
    stim = stim./max(stim) * 3 - 80;
    ts = (1:length(d)) .* si ./ 1000;
    % AP detection including spikelets
    [PKS,LOCS] = findpeaks(dats, "MinPeakDistance", 2*1000/si, "MinPeakProminence", 3);
    % control for spikelet amplitude
    plateau = median(dats(2718*1000/si:6704*1000/si));
    LOCS_ind_1 = (LOCS < 2718*1000/si | LOCS > 6704*1000/si) & PKS > -55;
    LOCS_ind_2 = LOCS > 2718*1000/si & LOCS < 6704*1000/si & PKS-plateau>10 & PKS > -55;
    LOCS_ind = LOCS_ind_1 | LOCS_ind_2;
    LOCS = LOCS(LOCS_ind);
    PKS = PKS(LOCS_ind);
    % control plot for AP detection
    figure;
    plot(ts, dats, "Color", "black"); hold on
    plot(ts, stim, "Color", "cyan")
    scatter(LOCS.*si./1000, PKS, "filled", "CData", [1 0 0])
    title(file, "Interpreter", "none");
    plotAesthetics(gca, 1, 14)
    xlabel("Time [ms]")
    ylabel("V_m [mV]")
    % collect data
    AP_max_ts{counter} = LOCS.*si./1000;
    AP_max_val{counter} = PKS;
    binned = histcounts(LOCS.*si./1000, "Normalization", "count", "BinEdges", [(2718-1000):500:(2718+4000+4000)]);
    AP_max_ts_bin(counter,:) = binned;
    counter = counter + 1;
end
AP_max_ts_bin_freq = AP_max_ts_bin*2;
[AP_max_ts_bin_freq_mean, ~, SEM, ~] = statistics(AP_max_ts_bin_freq, 1);

% plot NDNF IN population average firing rate
figure
plot((2718-1000):500:(2718+4000+3500), AP_max_ts_bin_freq_mean, "o", "Color", "black", "MarkerFaceColor", "black", "MarkerSize", 10); hold on
plot([(2718-1000):500:(2718+4000+3500);(2718-1000):500:(2718+4000+3500)],...
     [AP_max_ts_bin_freq_mean-SEM; AP_max_ts_bin_freq_mean+SEM], "LineWidth", 2, "Color", "black"); 
rectangle("Position", [2718, -3, 6703-2718, 2], "FaceColor", "cyan", "EdgeColor", "none")
xlim([1618 10318])
line(xlim, [0 0], "LineStyle", "--", "LineWidth", 2, "Color", [.5 .5 .5])
plotAesthetics(gca, 1, 14)
xlabel("Time [ms]"); ylabel("Firing Rate [Hz]")


%% Connectivity to PCs
path = './data/connectivity_NDNF_PC_50Hz';
lDir = dir(path);
traces = [];
stats = [];
for iter = 4:length(lDir)
    % load data
    file = lDir(iter).name;
    [d,si,h]=abfload(fullfile(path,file));
    d_raw = squeeze(d(:,1,:));
    % prepare data
    [filtDats1, filtDats2, filtDats3] = filtData(d_raw, 100, 500);
    d_avg = mean(filtDats1, 2);
    d_base = d_avg - mean(d_avg);
    traces(iter-3,:) = d_base;
    stim = d(:,2,1);
    ts = (1:size(d,1)) .* si ./ 1000;
    % analyze data
    [stats(iter-3,1) stats(iter-3,2) stats(iter-3,3) stats(iter-3,4) stats(iter-3,5) stats(iter-3,6)] = extractIPSCstats(ts, stim_on_ts, stim_off_ts, mean(filtDats2,2), si, 0);
end
stim = stim./max(stim).*3 - 10;

% plot IPSC traces
figure
% single traces
subplot(2,1,1)
    plot(ts, traces, "LineWidth", 1);
    hold on;
    plot(ts, stim, "Color", "cyan")
    line(xlim, [0 0], "LineStyle", "--", "LineWidth", 1, "Color", [.5 .5 .5])
    xlim([2500 9500]); ylim([-15 150])
    plotAesthetics(gca,1,14)
    xlabel("Time [ms]"); ylabel("Amplitude [pA]")
% average trace
subplot(2,1,2)
    plot(ts, mean(traces,1), "LineWidth", 2, "Color", "black");
    hold on;
    plot(ts, stim, "Color", "cyan")
    xlim([2500 9500]); ylim([-15 50])
	line(xlim, [0 0], "LineStyle", "--", "LineWidth", 1, "Color", [.5 .5 .5])
    plotAesthetics(gca,1,14)
    xlabel("Time [ms]"); ylabel("Amplitude [pA]")
    
% plot IPSC stats
figure
% max. amplitude
subplot(1,4,1) 
    bp1 = boxplot(stats(:,3))
        set(bp1, 'LineWidth', 2, "Color", "black");
    h = findobj('LineStyle','--'); 
        set(h, 'LineStyle','-');
    h = findobj(gca,'Tag','Box');
        patch(get(h,'XData'),get(h,'YData'), [.9 .9 .9],'FaceAlpha',.8);
    ylabel("Peak amplitude [pA]")
    plotAesthetics(gca, 1, 12)
    set(gca, "XTick", [])
% rise time
subplot(1,4,2)
    bp2 = boxplot(stats(:,4))
        set(bp2, 'LineWidth', 2, "Color", "black");
    h = findobj('LineStyle','--'); 
        set(h, 'LineStyle','-');
    h = findobj(gca,'Tag','Box');
        patch(get(h,'XData'),get(h,'YData'), [.9 .9 .9],'FaceAlpha',.8);
    ylabel("20-80 rise time [ms]")
    plotAesthetics(gca, 1, 12)
    set(gca, "XTick", [])
% decay time
subplot(1,4,3)
    bp3 = boxplot(stats(:,5))
        set(bp3, 'LineWidth', 2, "Color", "black");
    h = findobj('LineStyle','--'); 
        set(h, 'LineStyle','-');
    h = findobj(gca,'Tag','Box');
        patch(get(h,'XData'),get(h,'YData'), [.9 .9 .9],'FaceAlpha',.8);
    ylabel("80-20 decay time [ms]")
    plotAesthetics(gca, 1, 12)
    set(gca, "XTick", [])
% equilibrium amplitude
subplot(1,4,4)
    bp4 = boxplot(stats(:,6))
        set(bp4, 'LineWidth', 2, "Color", "black");
    h = findobj('LineStyle','--'); 
        set(h, 'LineStyle','-');
    h = findobj(gca,'Tag','Box');
        patch(get(h,'XData'),get(h,'YData'), [.9 .9 .9],'FaceAlpha',.8);
    ylabel("Equilibrium amplitude [pA]")
    plotAesthetics(gca, 1, 12)
    set(gca, "XTick", [])
    
    
% stats
median(stats(:,3)) % median max. amplitude
median(stats(:,4)) % median rise time
median(stats(:,5)) % median decay time
median(stats(:,6)) % median equi. amplitude
median(stats(:,6)./stats(:,3)) % equi./max. amplitude ratio


 %% subroutines
 % filter averaged traces (copied from "IPSCanalysis_singleTrials_v5.m")
 function [filtDats1, filtDats2, filtDats3] = filtData(dats, movAvgWin1, movAvgWin2)
    % 8-pole butterworth with Fc = 1000Hz
    G = [0.0038;
         0.0037;
         0.0036;
         0.0035;
         1.0000];
    SOS = [1.0000    2.0000    1.0000    1.0000   -1.9369    0.9523;
           1.0000    2.0000    1.0000    1.0000   -1.8551    0.8698;
           1.0000    2.0000    1.0000    1.0000   -1.7970    0.8112;
           1.0000    2.0000    1.0000    1.0000   -1.7670    0.7811];
    filtDats1 = filtfilt(SOS, G, dats);
    filtDats2 = movmean(filtDats1, movAvgWin1);
    filtDats3 = movmean(filtDats1, movAvgWin2);
 end

 % extract IPSC stats (based on "IPSCanalysis_singleTrials_v5.m")
 function [maxAmpl_ind, maxAmpl_ts, maxAmpl_val, rise2080, decay8020, equiAmpl] = extractIPSCstats(ts, LED_on, LED_off, dats_pA, si, controlPlot)
    % substractive normalization to 200ms pre-LED baseline
    ind_base = (LED_on-200)*1000/si:(LED_on)*1000/si;
    dats_pA = dats_pA - mean(dats_pA(ind_base));
    % max. amplitude amplitude & location
    [maxAmpl_val, maxAmpl_ind] = max(dats_pA(LED_on/si*1000:(LED_on+200)/si*1000));
    maxAmpl_ind = maxAmpl_ind + LED_on*1000/si - 1;
    maxAmpl_ts = maxAmpl_ind *si/1000;
    % rise & decay times
    dats_norm = dats_pA./maxAmpl_val;
    ind_20_rise = find(dats_norm(LED_on*1000/si:end) > 0.2, 1);
    ind_80_rise = find(dats_norm(LED_on*1000/si:end) > 0.8, 1);
    ind_20_fall = find(dats_norm(maxAmpl_ind:end) < 0.2, 1);
    ind_80_fall = find(dats_norm(maxAmpl_ind:end) < 0.8, 1);
    ind_20_rise = ind_20_rise - 1 + LED_on*1000/si;
    ind_80_rise = ind_80_rise - 1 + LED_on*1000/si;
    ind_20_fall = ind_20_fall - 1 + maxAmpl_ind;
    ind_80_fall = ind_80_fall - 1 + maxAmpl_ind;
    ts_20_rise = ind_20_rise * si/1000;
    ts_80_rise = ind_80_rise * si/1000;
    ts_20_fall = ind_20_fall * si/1000;
    ts_80_fall = ind_80_fall * si/1000;
    rise2080 = ts_80_rise - ts_20_rise;
    decay8020 = ts_20_fall - ts_80_fall;
    % equilibrium amplitude current
    ind_start = (LED_off - ((LED_off-LED_on)/2))*1000/si;
    ind_end = LED_off*1000/si;
    equiAmpl = mean(dats_pA(ind_start:ind_end));
    % analysis control plot
    if controlPlot
        figure;
        plot(ts, dats_pA); hold on
        plot(ts, dats_norm, "Color", "black")
        plot(maxAmpl_ts, maxAmpl_val, "o", "Color", "red")
        line([ts_20_rise ts_20_rise], ylim, "Color", "cyan")
        line([ts_80_rise ts_80_rise], ylim, "Color", "cyan")
        line([ts_80_fall ts_80_fall], ylim, "Color", "green")
        line([ts_20_fall ts_20_fall], ylim, "Color", "green")
        line(xlim, [0 0], "LineStyle", ":", "Color", "red")
        line(xlim, [1 1], "LineStyle", ":", "Color", "red")
        line(xlim, [.2 .2], "LineStyle", "--", "Color", "red")
        line(xlim, [.8 .8], "LineStyle", "--", "Color", "red")
        pause; close
    end
 end
 