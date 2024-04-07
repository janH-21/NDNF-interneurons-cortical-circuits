%% ----------------- prepare workspace and load data ----------------------
clear all; close all; clc;
cd("/main");
load('analysis/misc/aesthetics.mat')
si = 20;
% colors
turquois = [42, 162, 199]./255;
green = [103, 194, 58]./255;
cols = [[0 0 0]; turquois; green];
% LED calibration data
load('./analysis/Figure 3/LED_ndnfnpy_intersect_ChR2.mat');
fon_con_chr2 = 4.57;
% figure
F3_ndnfnpy_connectivity = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ----------------------- load connectvity data --------------------------
regexPattern = '_meanTrace.mat';
dirPath = './analysis/connectivity/';
addpath(genpath(dirPath)) 
% load mean traces by cell type
ndnfnpy_to_pc = readtable([dirPath 'NdnfNpy_PC_5_results.csv']);
[traces_pc, ts_pc] = readTraces([dirPath 'NDNFNPY_PC_meanTraces/'], regexPattern, 700000, 20);
[traces_pc_filt, ~, ~] = filtData(traces_pc', 100, 500);
ndnfnpy_to_pv = readtable([dirPath 'NdnfNpy_PV_5_results.csv']);
[traces_pv, ts_pv] = readTraces([dirPath 'NDNFNPY_PV_meanTraces/'], regexPattern, 700000, 20);
[traces_pv_filt, ~, ~] = filtData(traces_pv', 100, 500);
ndnfnpy_to_nonfs = readtable([dirPath 'NdnfNpy_nonFS_5_results.csv']);
[traces_nonfs, ts_nonfs] = readTraces([dirPath 'NDNFNPY_nonFS_meanTraces/'], regexPattern, 700000, 20);
[traces_nonfs_filt, ~, ~] = filtData(traces_nonfs', 100, 500);

% amplitude cutoff: recordings with smaller amplituide not included: SNR to
% high to reliable estimate parameters other than amplitude
ampl_cutoff = 10; % [pA]

% n connected
n_connect = [[sum(ndnfnpy_to_pc.isConnected == 1),...                    
                sum(ndnfnpy_to_pv.isConnected == 1),...
                sum(ndnfnpy_to_nonfs.isConnected == 1)];...
               [sum(ndnfnpy_to_pc.isConnected == 0),...
                sum(ndnfnpy_to_pv.isConnected == 0),...
                sum(ndnfnpy_to_nonfs.isConnected == 0)]]';
n_connect_percent = n_connect ./ sum(n_connect,2) * 100;

% load mean traces for pseudo-paired recordings
pairs = readtable([dirPath 'NDNFNPY_PC_PV_pairs_5.csv']);
[traces_pairs_pc, ts_pairs_pc] = readTraces([dirPath 'NDNFNPY_PC_PV_pairs_meanTraces/PC'], regexPattern, 700000, 20);
[traces_pairs_pc_filt, ~, ~] = filtData(traces_pairs_pc', 100, 500);
[traces_pairs_fs, ts_pairs_fs] = readTraces([dirPath 'NDNFNPY_PC_PV_pairs_meanTraces/PV'], regexPattern, 700000, 20);
[traces_pairs_fs_filt, ~, ~] = filtData(traces_pairs_fs', 100, 500);

% load mean PC CGP traces
CGP_pc = readtable([dirPath 'NdnfNpy_PC_CGP_5_results.csv']);
[traces_CGP_pre, ts_CGP_pre] = readTraces([dirPath 'NDNFNPY_PC_CGP_meanTraces/pre'], regexPattern, 700000, 20);
[traces_CGP_pre_filt, ~, ~] = filtData(traces_CGP_pre', 100, 500);
[traces_CGP_post, ts_CGP_post] = readTraces([dirPath 'NDNFNPY_PC_CGP_meanTraces/post'], regexPattern, 700000, 20);
[traces_CGP_post_filt, ~, ~] = filtData(traces_CGP_post', 100, 500);


%% --------------------- LED calibration ----------------------------------
ndnfnpy_cre_ap_data = table2array(LED_ndnf_npy(:,2:end-2));
[cf_mean_dats, ~, cf_sem_dats, ~] = statistics(ndnfnpy_cre_ap_data, 2);
subplot(5,5,2);
    line([0 10], [1 1], 'LineStyle', '-','LineWidth', 1,'Color',[.85 .85 .85]); hold on
    plot(LED_ndnf_npy.LED_mW_per_mm2, cf_mean_dats, '.-', 'Color', 'black', 'LineWidth', 1); hold on
    line([fon_con_chr2 fon_con_chr2], [-0.05 1.1], 'LineStyle', ':', 'Color', 'red', 'LineWidth', 1)
    plot([LED_ndnf_npy.LED_mW_per_mm2 LED_ndnf_npy.LED_mW_per_mm2]',...
         [cf_mean_dats+cf_sem_dats cf_mean_dats-cf_sem_dats]',...
         'LineWidth', 1, 'Color', 'black')
	plotAesthetics(gca, 1, font_size_large); 
	xlabel('LED power [mW/mm^2]'); ylabel('AP number')
    ylim([-0.05 1.1])
    text(0, 1.1, ' n=11',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold')
    
    
%% ------------------------- n connected ----------------------------------
subplot(5,5,3);
    b = bar(n_connect_percent, "BarLayout", "stacked");
    % aesthetics
    set(gca, 'YTick', 0:50:100) 
    xlim([0.5 3.5]); ylim([0 130]);
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'})
    plotAesthetics(gca, 1, font_size_large);
    ylabel('Percent of cells')
    set(b, "FaceAlpha", 0.8);
    % statistics
    isConnected = [ones(n_connect(1,1),1); zeros(n_connect(1,2),1);....
                   ones(n_connect(2,1),1); zeros(n_connect(2,2),1);....
                   ones(n_connect(3,1),1); zeros(n_connect(3,2),1)];
    cellType = [ones(sum(n_connect(1,:)),1);...
                ones(sum(n_connect(2,:)),1)*2;...
                ones(sum(n_connect(3,:)),1)*3];
    [TABLE,CHI2,P] = crosstab(isConnected, cellType); % p = 0.2003
    yLims = ylim;
    text(mean(xlim), 101, 'n.s.',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold')
    % legend
    text(0.5, 125, ' Connected',...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
         'FontSize', font_size_large, 'FontName', 'Arial', 'Color', [12, 19, 235]./255,...
         'FontWeight', 'bold');
    text(0.5, 135, ' Not connected',...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
         'FontSize', font_size_large, 'FontName', 'Arial', 'Color', 'red',...
         'FontWeight', 'bold');
    % n observations
    for iter = 1:3
        text(iter, yLims(1)+0.18*-diff(yLims), ['n=' num2str(sum(n_connect(iter,:)))],...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
            'FontSize', font_size_large, 'FontName', 'Arial',...
            'FontWeight', 'bold');
    end
 
    
%% ---------------------------- traces ------------------------------------
% PCs
subplot(5,30,19:22); 
    plot(ts_pc, mean(traces_pc_filt,2), 'Color', 'black', 'LineWidth', 1);
    xlim([2600 3600]); ylim([-5 55]); axis off
	text(mean(xlim), -1, 'PC',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
    text(mean(xlim), -8, ['n=' num2str(sum(ndnfnpy_to_pc.isConnected))],...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
% FS
subplot(5,30,23:26); 
    plot(ts_pv, mean(traces_pv_filt,2), 'Color', 'black', 'LineWidth', 1);
    xlim([2600 3600]); ylim([-5 55]); axis off
	text(mean(xlim), -1, 'PV',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
	text(mean(xlim), -8, ['n=' num2str(sum(ndnfnpy_to_pv.isConnected))],...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
% non-FS
subplot(5,30,27:30); 
    plot(ts_nonfs, mean(traces_nonfs_filt,2), 'Color', 'black', 'LineWidth', 1);
    xlim([2600 3600]); ylim([-5 55]); axis off
	text(mean(xlim), -1, 'non-FS',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
	text(mean(xlim), -8, ['n=' num2str(sum(ndnfnpy_to_nonfs.isConnected))],...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');

% scale bars
line([3200 3400],[40 40],'LineWidth', width_scalebar, 'Color', 'black')
line([3200 3200],[40 50],'LineWidth', width_scalebar, 'Color', 'black')
text(3300, 39, '200ms',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
text(3180, 45, '20pA',...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
 
 
%% -------------------- box plots + statistics ----------------------------
box_data = [ndnfnpy_to_pc;...
            ndnfnpy_to_pv;...
            ndnfnpy_to_nonfs];
box_groups = [ones(size(ndnfnpy_to_pc,1),1);...
              ones(size(ndnfnpy_to_pv,1),1)*2;...
              ones(size(ndnfnpy_to_nonfs,1),1)*3];
useInd = box_data.amplitude_1_pA > ampl_cutoff &...
         box_data.isConnected == 1; %[pA]

     
%% Amplitude
subplot(5,5,6); 
    bp_ampl = boxplot(box_data.amplitude_1_pA, box_groups); hold on
    set(bp_ampl, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
    sc_ampl = scatter(box_groups, box_data.amplitude_1_pA, 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel('Peak amplliitude [pA]'); 
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'});
    plotAesthetics(gca, 1, font_size_large);
    % add n observaations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p,tbl,stats] = kruskalwallis(box_data.amplitude_1_pA, box_groups, 'off'); % p = 5.3922e-07
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off"); % p_mc = 0 / 0.12 / 0.94
    counter = 1; yLims = ylim;
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
    
    
%% charge
subplot(5,5,7); 
    bp_charge= boxplot(box_data.charge_1_1000ms_pC, box_groups); hold on    
    set(bp_charge, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
    sc_charge = scatter(box_groups, box_data.charge_1_1000ms_pC, 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'});
	ylabel('Charge [pC]')
	plotAesthetics(gca, 1, font_size_large);
    % add n  observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p,tbl,stats] = kruskalwallis(box_data.charge_1_1000ms_pC, box_groups, 'off');
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off"); % p_mc = 0.009, 0.4, 0.94
     counter = 1; yLims = ylim;
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
    
    
%% t max
subplot(5,5,8); 
    bp_tmax = boxplot(box_data.tmaxAmpl_1_ms(useInd), box_groups(useInd)); hold on
    set(bp_tmax, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
    sc_max = scatter(box_groups(useInd), box_data.tmaxAmpl_1_ms(useInd), 'MarkerEdgeColor', 'black','SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'});
    ylabel('Latency to peak [ms]');
    plotAesthetics(gca, 1, font_size_large); 
    % add n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistcs
    [p,tbl,stats] = kruskalwallis(box_data.tmaxAmpl_1_ms(useInd), box_groups(useInd), 'off');
    text(mean(xlim), yLims(2), 'n.s.',...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');

    
%% rise time
subplot(5,5,9);
    bp_rise = boxplot(box_data.dtRise_1_ms(useInd), box_groups(useInd)); hold on
    set(bp_rise, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
    sc_rise = scatter(box_groups(useInd), box_data.dtRise_1_ms(useInd), 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'});
    ylabel('20-80 rise time [ms]');
    plotAesthetics(gca, 1, font_size_large);
    % add n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
	[p,tbl,stats] = kruskalwallis(box_data.dtRise_1_ms(useInd), box_groups(useInd), 'off'); % p = 0.31
    text(mean(xlim), yLims(2), 'n.s.',...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');

    
%% decay time
subplot(5,5,10); 
    bp_decay = boxplot(box_data.dtDecay_1_ms(useInd), box_groups(useInd)); hold on
    set(bp_decay, 'Color', 'black');
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
    sc_deacy = scatter(box_groups(useInd), box_data.dtDecay_1_ms(useInd), 'MarkerEdgeColor', 'black', 'SizeData', marker_size_box);
    % aesthetis
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'});
    ylabel('80-20 decay time [ms]');
    plotAesthetics(gca, 1, font_size_large);    
    % add n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(box_groups(useInd)==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
	[p,tbl,stats] = kruskalwallis(box_data.dtDecay_1_ms(useInd), box_groups(useInd), 'off'); 
	text(mean(xlim), yLims(2), 'n.s.',...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');

    
%% ------------------- PPR / short term depression ------------------------
% statistics
PPR_data = [box_data.amplitude_1_pA,...
            box_data.amplitude_2_pA,...
            box_data.amplitude_3_pA,...
            box_data.amplitude_4_pA];
PPR_data = PPR_data(box_data.isConnected == 1,:);
PPR_group = box_groups(box_data.isConnected == 1);
% PPR_data(isnan(PPR_data)) = 0;
[pc_mean, ~, pc_sem, ~] = statistics(PPR_data(PPR_group==1,:), 1);
[fs_mean, ~, fs_sem, ~] = statistics(PPR_data(PPR_group==2,:), 1);
[nonfs_mean, ~, nonfs_sem, ~] = statistics(PPR_data(PPR_group==3,:), 1);
% plot time course
subplot(5,5,11);
    % mean trace
    plot([1 2 3 4]-.05, pc_mean, '-', 'LineWidth', 1, 'Color',colCodes_con_inter(1,:)); hold on
    plot([1 2 3 4], fs_mean, '-', 'LineWidth', 1, 'Color',colCodes_con_inter(2,:)); hold on
    plot([1 2 3 4]+.05, nonfs_mean, '-', 'LineWidth', 1, 'Color',colCodes_con_inter(3,:)); hold on
    % SEM
    plot([1 2 3 4; 1 2 3 4]-.05, [pc_mean+pc_sem; pc_mean-pc_sem], 'LineWidth',1, 'Color',colCodes_con_inter(1,:));
    plot([1 2 3 4; 1 2 3 4], [fs_mean+fs_sem; fs_mean-fs_sem], 'LineWidth',1, 'Color',colCodes_con_inter(2,:));
    plot([1 2 3 4; 1 2 3 4]+.05, [nonfs_mean+nonfs_sem; nonfs_mean-nonfs_sem], 'LineWidth',1, 'Color',colCodes_con_inter(3,:));
    % add statistics
    counter = 1;
    for iter = 1:3
        [p,tbl,stats] = friedman(PPR_data(PPR_group==iter,:), 1, "off"); 
        mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
        for jter = 1:6
            if mc(jter,end) < 0.05
                line([mc(jter,1) mc(jter,2)], [73 73]+counter*2, 'LineWidth', 1, 'Color', cols(iter,:));
                counter  = counter + 1;
            end
        end
    end
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    ylabel('Amplitude [pA]'); xlabel('LED pulse'); set(gca, "XTick", 1:4);
    xlim([0.75 4.25])
    text(4, 65, 'PC',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_inter(1,:),...
     'FontWeight', 'bold');
    text(4, 57, 'PV',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_inter(2,:),...
     'FontWeight', 'bold');
    text(4, 49, 'non-FS',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large, "Color", colCodes_con_inter(3,:),...
     'FontWeight', 'bold');
    % add n  observations
    yLims = ylim;
    for iter = 1:3
        text(iter+.5, yLims(1)+0.18*-diff(yLims),....
            ['n=' num2str(sum(PPR_group==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold', "Color", colCodes_con_inter(iter,:));
    end

    
%% plot PPR
subplot(5,5,12);
    bp_ppr = boxplot(PPR_data(:,2)./PPR_data(:,1), PPR_group); hold on
    set(bp_ppr, 'Color', 'black')
    plotAesthetics(gca, 1, font_size_large);
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_con_inter(3,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_con_inter(2,:),'FaceAlpha',.8);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_con_inter(1,:),'FaceAlpha',.6);
	sc_ppr = scatter(PPR_group, PPR_data(:,2)./PPR_data(:,1), 'MarkerEdgeColor', 'black','SizeData', marker_size_box);
    % Aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", {'PC', 'PV', 'non-FS'})
    ylabel('Paired-pulse ratio');
    plotAesthetics(gca, 1, font_size_large);
    % add n observations
    yLims = ylim;
    for iter = 1:3
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(PPR_group==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", "FontSize", font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p,tbl,stats] = kruskalwallis(PPR_data(:,2)./PPR_data(:,1), PPR_group,"off"); % p = 0.1363
    text(mean(xlim), yLims(2), 'n.s.',...
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center',...
         'FontName', 'Arial', 'FontSize', font_size_large,...
         'FontWeight', 'bold');
 

%% ---------------------------- pairs -------------------------------------
% plot 
subplot(5,10,27);
    bp_decay = boxplot([pairs.pc_ampl_1_pA pairs.pv_ampl_1_pA]); hold on
	set(bp_decay, "Color", "black");
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),cols(3-jter,:),'FaceAlpha',.5);
    end
    plot([1 2], [pairs.pc_ampl_1_pA pairs.pv_ampl_1_pA],...
        '-', "Color", "black");
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    plotAesthetics(gca, 1, font_size_large);
    set(gca, "XTickLabel", ["PC", "PV"]);
    ylabel("Amplitude [pA]");
    % statistics: wilcoxon sign-rank test
    p = signrank(pairs.pc_ampl_1_pA, pairs.pv_ampl_1_pA); % p = 0.0039
    line([1 2],[130 130],'Color', 'black', 'LineWidth', 1);
    text(1.5, 130, '**',...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
    % add n observations
    yLims = ylim;
    text(1.5, yLims(1)+0.1*-diff(yLims),....
        ['n=' num2str(size(pairs,1))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", "FontSize", font_size_large,...
        'FontWeight', 'bold');

    
%% plot  traces
traces_pairs_pc_filt_2 = traces_pairs_pc_filt(1.3e5:1.8e5,:);
pairs_pc_mean = mean(traces_pairs_pc_filt_2,2);
pairs_pc_mean = pairs_pc_mean - median(pairs_pc_mean);
traces_pairs_fs_filt_2 = traces_pairs_fs_filt(1.3e5:1.8e5,:);
pairs_fs_mean = mean(traces_pairs_fs_filt_2,2);
pairs_fs_mean = pairs_fs_mean - median(pairs_fs_mean);
ts_pc_plot = (1:size(traces_pairs_pc_filt_2,1)).*si/1000; % ts_pairs_pc(1.3e5:1.8e5);
subplot(5,10,25:26);
    plot(ts_pc_plot, pairs_pc_mean, 'LineWidth', 1, 'Color', cols(1,:)); hold on
    plot(ts_pc_plot, pairs_fs_mean, 'LineWidth', 1, 'Color', cols(2,:));
    line([800 1000], [25 25], 'LineWidth', width_scalebar, 'Color', 'black');
    line([800 800], [25 33], 'LineWidth', width_scalebar, 'Color', 'black');
    text(900, 24.5, '200ms',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
     'FontWeight', 'bold');
	text(795, 30, '5pA',...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
     'FontWeight', 'bold');
    text(1000, 15, 'PC',...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
     'FontWeight', 'bold');
    text(1000, 16, 'PV',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large, 'Color', cols(2,:),...
     'FontWeight', 'bold');
 	text(mean(xlim), -2, 'mean traces',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
    axis off
    pair = gca;
    pos_pair = pair.Position;
  	set(pair, "Position", pos_pair + [0 0 -.03 0])
    
    
%% ---------------------------- CGP data ----------------------------------
% get indices: every second is post-CGP
pre_ind = zeros(size(CGP_pc,1),1);
post_ind = ones(size(CGP_pc,1),1);
groups = ones(size(CGP_pc,1),1)*2;
for iter = 1:2:length(pre_ind)
    pre_ind(iter) = pre_ind(iter) + 1;
    post_ind(iter) = post_ind(iter) - 1;
    groups(iter) = groups(iter) - 1;
end
% plot results
subplot(5,10,30);
    bp_decay = boxplot(CGP_pc.dtDecay_1_ms, groups); hold on
    set(bp_decay, "Color", "black");
    h = findobj(gca,'Tag','Box');
    patch(get(h(2),'XData'),get(h(2),'YData'),cols(1,:),'FaceAlpha',.5);
    patch(get(h(1),'XData'),get(h(1),'YData'),[181, 24, 47]./255,'FaceAlpha',.5);
    plot([1 2], [CGP_pc.dtDecay_1_ms(boolean(pre_ind))...
                 CGP_pc.dtDecay_1_ms(boolean(post_ind))],...
        '-', "Color", "black");
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", ["pre", "post"]);
    ylabel("80-20 decay time [ms]");
    plotAesthetics(gca, 1, font_size_large);
    % statistics: wilcoxon sign-rank test
    p_CGP = signrank(CGP_pc.dtDecay_1_ms(boolean(pre_ind)),...
                     CGP_pc.dtDecay_1_ms(boolean(post_ind))); % p = 0.0391
    line([1 2],[485 485], "Color", "black", "LineWidth", 1)
    text(1.5, 485, '*',...
         "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
         "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
    % add n observations
    yLims = ylim;
    text(1.5, yLims(1)+0.1*-diff(yLims),....
        ['n=' num2str(sum(groups==1))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", "FontSize", font_size_large,...
        'FontWeight', 'bold');
 
% plot normalized traces
traces_CGP_pre_filt_2 = traces_CGP_pre_filt(1.3e5:1.8e5,:);
traces_CGP_pre_filt_3 = traces_CGP_pre_filt_2 ./ max(traces_CGP_pre_filt_2);
pre_CPG_mean = mean(traces_CGP_pre_filt_3,2);
traces_CGP_post_filt_2 = traces_CGP_post_filt(1.3e5:1.8e5,:);
traces_CGP_post_filt_3 = traces_CGP_post_filt_2 ./ max(traces_CGP_post_filt_2);
post_CPG_mean = mean(traces_CGP_post_filt_3,2);
ts_CGP_plot = ts_CGP_pre(1.3e5:1.8e5);
subplot(5,10,28:29);
    plot(ts_CGP_plot, pre_CPG_mean, 'LineWidth', 1, 'Color', cols(1,:)); hold on
    plot(ts_CGP_plot, post_CPG_mean, 'LineWidth', 1, 'Color', [181, 24, 47]./255);
    line([3300 3500], [0.85 0.85], 'LineWidth', width_scalebar, 'Color', 'black');
    text(3400, 0.84, '200ms',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
     'FontWeight', 'bold');
    text(3500, 0.5, 'pre CGP',...
     'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large,...
     'FontWeight', 'bold');
    text(3500, 0.48, 'post CGP',...
     'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
     'FontName', 'Arial', 'FontSize', font_size_large, 'Color', [181, 24, 47]./255,...
     'FontWeight', 'bold');
 	text(mean(xlim), -0.05, 'norm. mean traces',...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center',...
        'FontName', 'Arial', 'FontSize', font_size_large,...
        'FontWeight', 'bold');
    axis off;
    cgp = gca;
    pos_cgp = cgp.Position;
    set(cgp, "Position", pos_cgp + [0 0 -.03 0])


%% save figure
set(F3_ndnfnpy_connectivity,'renderer','Painters')
% print('F5_main_intersect_connectivity_v006.pdf','-dpdf')