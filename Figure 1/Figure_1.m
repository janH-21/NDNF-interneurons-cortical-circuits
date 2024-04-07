%% prepare workspace and load data
clear all; close all; clc;
% working directory
cd('/main/');
% aesthetis
load('./analysis/misc/aesthetics.mat')
load('./analysis/misc/IPSC_lowpassFilter.mat')
colCodes_con_ndnf = [[0.9882    0.6157    0.0118];
                     [0.1647    0.6353    0.7804];
                     [0.4039    0.7608    0.2275];
                     [0.8000    0.1000    0.1000]];
% LED calibration data
load('./analysis/Figure 1/LED_ndnf_cre_ChR2.mat');
load('./analysis/Figure 1/LED_ndnf_flpo_chrimsonR.mat');
% chosen LED powers in mm^2
dio_chr2 = 5.08;
fdio_chrimsonr = 1.82;
% LEDtimes [ms]
LEDtimes = [2.7188 3.7188 4.7188 5.7188]*1000;

%% load connectvity data
regexPattern = '_meanTrace.mat';
% NDNF > PV
dirPath = './analysis/connectivity/';
addpath(genpath(dirPath)) 
 ndnf_to_pv = readtable([dirPath 'Ndnf_PV_5_results.csv']);
 [traces_pv, ts_pv] = readTraces([dirPath 'NDNF_PV_meanTraces/'], regexPattern, 700000, 20);
 [traces_pv_filt, ~, ~] = filtData(traces_pv', 100, 500);
% NDNF > non-FS
 ndnf_to_nonfs = readtable([dirPath 'Ndnf_nonFS_5_results.csv']);
 [traces_nonfs, ts_nonfs] = readTraces([dirPath 'NDNF_nonFS_meanTraces/'], regexPattern, 700000, 20);
 [traces_nonfs_filt, ~, ~] = filtData(traces_nonfs', 100, 500);
% NDNF > VIP
 ndnf_to_vip = readtable([dirPath 'Ndnf_VIP_5_results.csv']);
 [traces_vip, ts_vip] = readTraces([dirPath 'NDNF_VIP_meanTraces/'], regexPattern, 700000, 20);
 [traces_vip_filt, ~, ~] = filtData(traces_vip', 100, 500);
% NDNF > SST
 ndnf_to_sst = readtable([dirPath 'Ndnf_SST_5_results.csv']);

% amplitude cutoff: recordings with smaller amplituide not included: SNR to
% high to reliable estimate parameters other than amplitude
ampl_cutoff = 10; % [pA]
% n cells connnected
% collect data
n_connect = [[sum(ndnf_to_vip.isConnected == 1),...
              sum(ndnf_to_pv.isConnected == 1),...                   
              sum(ndnf_to_nonfs.isConnected == 1),...
              sum(ndnf_to_sst.isConnected == 1)];...
             [sum(ndnf_to_vip.isConnected == 0),...
              sum(ndnf_to_pv.isConnected == 0),...
              sum(ndnf_to_nonfs.isConnected == 0),...
              sum(ndnf_to_sst.isConnected == 0)]]';
n_connect_percent = n_connect ./ sum(n_connect,2) * 100;

% figure
F1_ndnf_connectivity = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");

%% ------------------------ plot example cell -----------------------------
% 2019_11_15_NdnfVIP_0008.abf
subplot(6,5,2:3);
plot(NaN); hold on
    % plot LED linies
    for iter = 1:4
        rectangle('Position', [LEDtimes(iter), -10, 0.5, 85],....
            'FaceColor', 'cyan', 'EdgeColor', 'cyan'); hold on;
    end
    text(LEDtimes(1)+50, 72, 'LED',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',....
        'FontSize', font_size_large, 'FontName', 'Arial', 'Color', 'cyan',...
        'FontWeight', 'bold')
    % plot trace
    plot(ts_vip, filtfilt(SOS,G,traces_vip(3, :)), 'Color', 'black');
    xlim([2600 6600]); ylim([-10 72]);
    text(mean(xlim), -12, 'Example VIP mean IPSC',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
    % scale bars
    line([6250 6500],[50 50],'LineWidth', width_scalebar, 'Color', 'black')
    line([6250 6250],[50 60],'LineWidth', width_scalebar, 'Color', 'black')
    text(mean([6250 6500]), 49, '250ms',...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');
    text(6250, 55, '10pA ',...
         'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');
axis off

%% ------------------ plot LED calibration curves -------------------------
% ChR2
ndnf_cre_ap_data = table2array(LED_ndnf_cre(:,2:end-2));
[cre_mean_dats, ~, cre_sem_dats, ~] = statistics(ndnf_cre_ap_data, 2);
subplot(6,5,6); hold on;
    line([0 10], [1 1], 'LineStyle', '-','LineWidth', width_scalebar,'Color',[.85 .85 .85])
    plot(LED_ndnf_cre.LED_mW_per_mm2, cre_mean_dats, '.-', 'Color', 'black', 'LineWidth', width_scalebar); hold on
    line([dio_chr2 dio_chr2], [-0.05 1.6], 'LineStyle', ':', 'Color', 'red', 'LineWidth', width_scalebar)
    plot([LED_ndnf_cre.LED_mW_per_mm2 LED_ndnf_cre.LED_mW_per_mm2]',...
         [cre_mean_dats+cre_sem_dats cre_mean_dats-cre_sem_dats]',...
         'LineWidth', width_scalebar, 'Color', 'black')
	plotAesthetics(gca, 1, font_size_large); 
	xlabel('LED power [mW/mm^2]'); ylabel(['AP number'])
    text(0.5, 1.5, 'n=8',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold');
    ylim([-0.05 1.6])
   
% ChrimsonR
ndnf_flpo_ap_data = table2array(LED_ndnf_flpo(:,2:end-1));
[flpo_mean_dats, ~, flpo_sem_dats, ~] = statistics(ndnf_flpo_ap_data, 2);  
subplot(6,5,11); hold on;
    line([0 2.3], [1 1], 'LineStyle', '-','LineWidth',width_scalebar,'Color',[.85 .85 .85])
    plot(LED_ndnf_flpo.LED_mW_per_mm2, flpo_mean_dats, '.-', 'Color', 'black', 'LineWidth', width_scalebar); hold on
    line([fdio_chrimsonr fdio_chrimsonr], [-0.05 1.6], 'LineStyle', ':', 'Color', 'red', 'LineWidth', width_scalebar)
    plot([LED_ndnf_flpo.LED_mW_per_mm2 LED_ndnf_flpo.LED_mW_per_mm2]',...
         [flpo_mean_dats+flpo_sem_dats flpo_mean_dats-flpo_sem_dats]',...
         'LineWidth', width_scalebar, 'Color', 'black')
    plotAesthetics(gca, 1, font_size_large); 
    xlabel('LED power [mW/mm^2]'); ylabel(['AP number'])
    text(0.1, 1.5, 'n=5',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold');
    ylim([-0.05 1.6])
      
%% ------------------------- bar n connected ------------------------------
subplot(6,5,7)
    b = bar(n_connect_percent, "BarLayout", "stacked");
    yLims = ylim;
    for iter = 1:4
        text(iter, -13, ['n=' num2str(sum(n_connect(iter,:)))],...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
            'FontSize', font_size_large, 'FontName', 'Arial',...
            'FontWeight', 'bold');
    end
    set(b, "FaceAlpha", 0.8);
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS', 'SST'})
    plotAesthetics(gca, 1, font_size_large);
    ylabel('Percent of cells')

% statistics
isConnected = [ones(n_connect(1,1),1); zeros(n_connect(1,2),1);...
               ones(n_connect(2,1),1); zeros(n_connect(2,2),1);...
               ones(n_connect(3,1),1); zeros(n_connect(3,2),1);...
               ones(n_connect(4,1),1); zeros(n_connect(4,2),1)];
cellType = [ones(sum(n_connect(1,:)),1);...
            ones(sum(n_connect(2,:)),1)*2;...
            ones(sum(n_connect(3,:)),1)*3;...
            ones(sum(n_connect(4,:)),1)*4];
[TABLE,CHI2,P] = crosstab(isConnected, cellType); % p = 0.258
xlim([0.5 4.5]); ylim([0 130])
text(mean(xlim), 101, 'n.s.',...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
    'FontSize', font_size_large, 'FontName', 'Arial',...
	 'FontWeight', 'bold')
text(0.6, 115, 'Connected',...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
     'FontSize', font_size_large, 'FontName', 'Arial', 'Color', [12, 19, 235]./255,...
	 'FontWeight', 'bold');
text(0.6, 125, 'Not connected',...
     'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
     'FontSize', font_size_large, 'FontName', 'Arial', 'Color', 'red',...
	 'FontWeight', 'bold');
 
 
%% ----------------------------- traces -----------------------------------
subplot(5,30,19:22); % VIP
    plot(ts_vip, mean(traces_vip_filt,2), 'Color', 'black', "LineWidth", 1);
    xlim([2600 3600]); ylim([-15 40]); axis off
	text(mean(xlim), -5, 'VIP',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
	text(mean(xlim), -10, ['n=' num2str(sum(ndnf_to_vip.isConnected))],...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	'FontWeight', 'bold');
subplot(5,30,23:26); %PV
    plot(ts_pv, mean(traces_pv_filt,2), 'Color', 'black', "LineWidth", 1);
    xlim([2600 3600]); ylim([-15 40]); axis off
    text(mean(xlim), -5, 'PV',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
	text(mean(xlim), -10, ['n=' num2str(sum(ndnf_to_pv.isConnected))],...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
subplot(5,30,27:30); % non-FS
    plot(ts_nonfs, mean(traces_nonfs_filt,2), 'Color', 'black', "LineWidth", 1);
    xlim([2600 3600]); ylim([-15 40]); axis off
	text(mean(xlim), -5, 'non-FS',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
	text(mean(xlim), -10, ['n=' num2str(sum(ndnf_to_nonfs.isConnected))],...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');

 
% scale bars
line([3200 3450],[30 30],'LineWidth', width_scalebar, 'Color', 'black')
line([3200 3200],[30 40],'LineWidth', width_scalebar, 'Color', 'black')
text(3325, 29, '250ms',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
text(3200, 35, '10pA ',...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
	 'FontWeight', 'bold');
 

%% ---------------------- box plots + statistics --------------------------
box_data_ac = [ndnf_to_vip;
               ndnf_to_pv;
               ndnf_to_nonfs;
               ndnf_to_sst];
box_groups_ac = [ones(size(ndnf_to_vip,1),1)*1;
                 ones(size(ndnf_to_pv,1),1)*2;
                 ones(size(ndnf_to_nonfs,1),1)*3;
                 ones(size(ndnf_to_sst,1),1)*4];

%% amplitude
subplot(6,5,8);
    bp_ampl = boxplot(box_data_ac.amplitude_1_pA, box_groups_ac); hold on
    set(bp_ampl, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(5-jter,:),'FaceAlpha',.8);
    end
    sc_ampl = scatter(box_groups_ac, box_data_ac.amplitude_1_pA, 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel('Peak amplliitude [pA]'); 
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS', 'SST'});
	plotAesthetics(gca, 1, font_size_large); 
    % add n observations
    yLims = ylim;
    for iter = 1:4
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups_ac==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p_ampl, tbl, stats] = kruskalwallis(box_data_ac.amplitude_1_pA, box_groups_ac, 'off');
	mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off"); % p_pv_sst = 0.0305; p_vip_sst = 0.0043
    counter = 1;
    for iter = 1:size(mc,1)
        if mc(iter,end) < 0.05
            yPos = yLims(2) + diff(yLims) * 0.15 * counter;
            line([mc(iter,1) mc(iter,2)], [yPos yPos], 'LineWidth', 1, 'Color', 'black');
            if mc(iter, end) < 0.001; stars = "***";
            elseif mc(iter, end) < 0.01; stars = "**";
            else; stars = "*"; end
            text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
            counter  = counter + 1;
        end
    end
	ylim([yLims(1) yPos])
%     text(mean(xlim), yLims(2)+(counter)*.15*diff(yLims), ['p = ' num2str(p_ampl)],... % p = 0.0053
%          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
%          'FontName', 'Arial', 'FontSize', font_size_large,...
%          'FontWeight', 'bold');

%% charge
subplot(6,5,9);
    bp_charge= boxplot(box_data_ac.charge_1_1000ms_pC, box_groups_ac); hold on    
    set(bp_charge, 'Color', 'black', 'LineStyle', '-');
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(5-jter,:),'FaceAlpha',.8);
    end
    sc_charge = scatter(box_groups_ac, box_data_ac.charge_1_1000ms_pC, 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS', 'SST'});
    ylabel('Charge [pC]')
	plotAesthetics(gca, 1, font_size_large); 
    % n observations
    yLims = ylim;
    for iter = 1:4
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups_ac==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p_charge, tbl, stats] = kruskalwallis(box_data_ac.charge_1_1000ms_pC, box_groups_ac, 'off');
	mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off"); % p_pv_sst = 0.0365; p_vip_sst = 0.00429
    counter = 1;
    for iter = 1:size(mc,1)
        if mc(iter,end) < 0.05
            yPos = yLims(2) + diff(yLims) * 0.15 * counter;
            line([mc(iter,1) mc(iter,2)], [yPos yPos], 'LineWidth', 1, 'Color', 'black');
            if mc(iter, end) < 0.001; stars = "***";
            elseif mc(iter, end) < 0.01; stars = "**";
            else; stars = "*"; end
            text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
            counter  = counter + 1;
        end
    end
	ylim([yLims(1) yPos])
%     text(mean(xlim), yLims(2)+(counter)*.15*diff(yLims), ['p = ' num2str(p_charge)],... % p = 0.0232
%          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
%          'FontName', 'Arial', 'FontSize', font_size_large,...
%          'FontWeight', 'bold');

%% format data for kinetics
box_data_kin = [ndnf_to_vip;...
                ndnf_to_pv;...
                ndnf_to_nonfs];
box_groups_kin = [ones(size(ndnf_to_vip,1),1)*1;...
                  ones(size(ndnf_to_pv,1),1)*2;...
                  ones(size(ndnf_to_nonfs,1),1)*3];
useInd = box_data_kin.amplitude_1_pA > ampl_cutoff &...
         box_data_kin.isConnected == 1;
     
%% time to peak
subplot(6,5,10); % t max
    bp_tmax = boxplot(box_data_kin.tmaxAmpl_1_ms(useInd), box_groups_kin(useInd)); hold on
    set(bp_tmax, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(4-jter,:),'FaceAlpha',.8);
    end
    sc_max = scatter(box_groups_kin(useInd), box_data_kin.tmaxAmpl_1_ms(useInd), 'MarkerEdgeColor', 'black','SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS'});
    ylabel('Latency to peak [ms]');
    plotAesthetics(gca, 1, font_size_large);
    % n observations
	yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups_kin(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
	p_tmax = kruskalwallis(box_data_kin.tmaxAmpl_1_ms(useInd), box_groups_kin(useInd), 'off'); % p = 0.1583
% 	text(mean(xlim), yLims(2), ['p = ' num2str(p_tmax)],...
%          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
%          'FontName', 'Arial', 'FontSize', font_size_large,...
%          'FontWeight', 'bold');
	text(mean(xlim), yLims(2), 'n.s.',...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');

%% rise time
subplot(6,5,12);% rise time
    bp_rise = boxplot(box_data_kin.dtRise_1_ms(useInd), box_groups_kin(useInd)); hold on
    set(bp_rise, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(4-jter,:),'FaceAlpha',.8);
    end
    sc_rise = scatter(box_groups_kin(useInd), box_data_kin.dtRise_1_ms(useInd), 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS'});
    ylabel('20-80 rise time [ms]')
    plotAesthetics(gca, 1, font_size_large); 
    % n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups_kin(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    p_rise = kruskalwallis(box_data_kin.dtRise_1_ms(useInd), box_groups_kin(useInd), 'off'); % p = 0.2900
% 	text(mean(xlim), yLims(2), ['p = ' num2str(p_rise)],...
%          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
%          'FontName', 'Arial', 'FontSize', font_size_large,...
%          'FontWeight', 'bold');
	text(mean(xlim), yLims(2), 'n.s.',...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');

%% decay time
subplot(6,5,13);
    bp_decay = boxplot(box_data_kin.dtDecay_1_ms(useInd), box_groups_kin(useInd)); hold on
    set(bp_decay, 'Color', 'black', 'LineStyle', '-');
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(4-jter,:),'FaceAlpha',.8);
    end
   	sc_deacy = scatter(box_groups_kin(useInd), box_data_kin.dtDecay_1_ms(useInd), 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'VIP', 'PV', 'non-FS'});
	ylabel('80-20 decay time [ms]');
    plotAesthetics(gca, 1, font_size_large); 
    % n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups_kin(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    p_decay = kruskalwallis(box_data_kin.dtDecay_1_ms(useInd), box_groups_kin(useInd), 'off'); % p_decay = 0.6814
% 	text(mean(xlim), yLims(2), ['p = ' num2str(p_decay)],...
%          'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
%          'FontName', 'Arial', 'FontSize', font_size_large,...
%          'FontWeight', 'bold');
	text(mean(xlim), yLims(2), 'n.s.',...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');

%% ----------------------- PPR / short term depression ---------------------
% format data
box_data_ppr = [ndnf_to_vip;...
                ndnf_to_pv;...
                ndnf_to_nonfs];
box_groups_ppr = [ones(size(ndnf_to_vip,1),1)*1
                  ones(size(ndnf_to_pv,1),1)*2;...
                  ones(size(ndnf_to_nonfs,1),1)*3];
PPR_data = [box_data_ppr.amplitude_1_pA,...
            box_data_ppr.amplitude_2_pA,...
            box_data_ppr.amplitude_3_pA,...
            box_data_ppr.amplitude_4_pA];
PPR_groups = box_groups_ppr(box_data_ppr.isConnected==1);
PPR_data = PPR_data(box_data_ppr.isConnected==1,:);
PPR_data(isnan(PPR_data)) = 0;
[vip_mean, ~, vip_sem, ~] = statistics(PPR_data(PPR_groups==1,:), 1);
[pv_mean, ~, pv_sem, ~] = statistics(PPR_data(PPR_groups==2,:), 1);
[nonfs_mean, ~, nonfs_sem, ~] = statistics(PPR_data(PPR_groups==3,:), 1);


% plot
subplot(6,5,14);
    plot([1 2 3 4]-.05, pv_mean, '-', 'LineWidth', width_scalebar, 'Color',turquois); hold on
    plot([1 2 3 4], nonfs_mean, '-', 'LineWidth', width_scalebar, 'Color',green); hold on
    plot([1 2 3 4]+.05, vip_mean, '-', 'LineWidth', width_scalebar, 'Color', orange); hold on
    plot([1 2 3 4; 1 2 3 4]-.05, [pv_mean+pv_sem; pv_mean-pv_sem], 'LineWidth',width_scalebar, 'Color',turquois);
    plot([1 2 3 4; 1 2 3 4], [nonfs_mean+nonfs_sem; nonfs_mean-nonfs_sem], 'LineWidth',width_scalebar, 'Color',green);
    plot([1 2 3 4; 1 2 3 4]+.05, [vip_mean+vip_sem; vip_mean-vip_sem], 'LineWidth',width_scalebar, 'Color', orange);
% aesthetics
plotAesthetics(gca, 1, font_size_large);
ylabel('Amplitude [pA]'); xlabel('LED pulse (1Hz)')
set(gca, "XTick", 1:4)
xlim([0.75 4.25])
% statistics
counter = 1;
for iter = 1:3
    [p,tbl,stats] = friedman(PPR_data(PPR_groups==iter,:), 1, "off");
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    for jter = 1:6
        if mc(jter,end) < 0.05
            line([mc(jter,1) mc(jter,2)], [43 43]+counter*2, 'LineWidth', width_scalebar, 'Color', colCodes_con_ndnf(iter,:));
            counter  = counter + 1;
        end
    end
end
% n observations
yLims = ylim;
for iter = 1:3
    text(iter+.5, yLims(1)+0.23*-diff(yLims),....
        ['n=' num2str(sum(PPR_groups==iter))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", 'FontSize', font_size_large,...
        'FontWeight', 'bold', 'Color', colCodes_con_ndnf(iter,:));
end
% legend
text(4, 45, 'VIP',...
 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
 'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_ndnf(1,:),...
 'FontWeight', 'bold');
text(4, 40, 'PV',...
 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
 'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_ndnf(2,:),...
 'FontWeight', 'bold');
text(4, 35, 'non-FS',...
 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
 'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_ndnf(3,:),...
 'FontWeight', 'bold');


%% PPR
subplot(6,5,15); 
    bp_ppr = boxplot(PPR_data(:,2)./PPR_data(:,1), PPR_groups); hold on
    set(bp_ppr, 'Color', 'black')
    plotAesthetics(gca, 1, font_size_large);
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(4-jter,:),'FaceAlpha',.8);
    end
    sc_ppr = scatter(PPR_groups, PPR_data(:,2)./PPR_data(:,1), 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
% Aesthetics
h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
ylim([-0.05 1.3+counter/100])
set(gca, "XTickLabel", {'PV', 'non-FS', 'VIP'})
ylabel('Paired-pulse ratio')
% n observations
yLims = ylim;
for iter = 1:3
    text(iter, yLims(1)+0.15*-diff(yLims),....
        ['n=' num2str(sum(PPR_groups==iter))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", 'FontSize', font_size_large,...
        'FontWeight', 'bold');
end
% statistics
counter = 1;
[p,tbl,stats] = kruskalwallis(PPR_data(:,2)./PPR_data(:,1), PPR_groups,"off");
mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
yLims = ylim;
counter = 1;
if p < 0.05
    for iter = 1:size(mc,1)
        if mc(iter, end) < 0.05
            yPos = yLims(2) + diff(yLims) * 0.15 * counter;
            line([mc(iter, 1) mc(iter, 2)],[yPos yPos],"Color", "black", "LineWidth", width_scalebar);
            if mc(iter, end) < 0.001; stars = "***";
            elseif mc(iter, end) < 0.01; stars = "**";
            else; stars = "*"; end
            text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
            counter = counter + 1;
        end
    end
end
ylim([yLims(1) yPos])


%% save figure
set(F1_ndnf_connectivity,'renderer','Painters')
% print('F1_ndnf_connectivity.pdf','-dpdf')