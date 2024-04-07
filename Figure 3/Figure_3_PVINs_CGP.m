%% prepare workspace and load data
%clear all; close all; clc;
cd("/main/")
load('./analysis/misc/aesthetics.mat')
cols = [[0 0 0]; turquois; green];

%% load data
% load result tables
table_ndnf_fs = readtable('./analysis/connectivity/Ndnf_PV_CGP_5_results.csv');
table_ndnfnpy_fs = readtable('./analysis/connectivity/NdnfNpy_PV_CGP_5_results.csv');

% format IPSC chars
exp_idx = [ones(size(table_ndnf_fs.dtDecay_1_ms)); ones(size(table_ndnfnpy_fs.dtDecay_1_ms))*2];
all_fs = [table_ndnf_fs.dtDecay_1_ms; table_ndnfnpy_fs.dtDecay_1_ms];
cgp_idx = ones(size(all_fs));
cgp_idx(2:2:end) = 2;

% load traces
regexPattern = '_meanTrace.mat';
dirPath = './analysis/connectivity/';
addpath(genpath(dirPath));
[traces_ndnf_pre, ts_ndnf_pre] = readTraces([dirPath 'NDNF_PV_CGP_meanTraces/pre'], regexPattern, 700000, 20);
[traces_ndnf_pre_filt, ~, ~] = filtData(traces_ndnf_pre', 100, 500);
[traces_ndnfnpy_pre, ts_ndnfnpy_pre] = readTraces([dirPath 'NDNFNPY_PV_CGP_meanTraces/pre'], regexPattern, 700000, 20);
[traces_ndnfnpy_pre_filt, ~, ~] = filtData(traces_ndnfnpy_pre', 100, 500);    
[traces_ndnf_post, ts_ndnf_post] = readTraces([dirPath 'NDNF_PV_CGP_meanTraces/post'], regexPattern, 700000, 20);
[traces_ndnf_post_filt, ~, ~] = filtData(traces_ndnf_post', 100, 500);
[traces_ndnfnpy_post, ts_ndnfnpy_ppost] = readTraces([dirPath 'NDNFNPY_PV_CGP_meanTraces/post'], regexPattern, 700000, 20);
[traces_ndnfnpy_post_filt, ~, ~] = filtData(traces_ndnfnpy_post', 100, 500);  

% format traces
traces_pre = [traces_ndnf_pre_filt(1.3e5:1.8e5,:), traces_ndnfnpy_pre_filt(1.3e5:1.8e5,:)];
traces_post = [traces_ndnf_post_filt(1.3e5:1.8e5,:), traces_ndnfnpy_post_filt(1.3e5:1.8e5,:)];
traces_pre_norm = traces_pre ./ max(traces_pre, [], 1);
traces_pre_mean = mean(traces_pre_norm,2);
traces_post_norm = traces_post ./ max(traces_post, [], 1);
traces_post_mean = mean(traces_post_norm,2);
ts_CGP_plot = ts_ndnf_pre(1.3e5:1.8e5);


%% boxplot
subplot(1,2,1)
    bp_cgp_fs = boxplot(all_fs, cgp_idx); hold on
        set(bp_cgp_fs, "Color", "black");
        h = findobj('LineStyle','--'); set(h, 'LineStyle','-');
        h = findobj(gca,'Tag','Box');
        patch(get(h(2),'XData'),get(h(2),'YData'),"black",'FaceAlpha',.5);
        patch(get(h(1),'XData'),get(h(1),'YData'),"red",'FaceAlpha',.5);
    plot([1 1 1 1 1 1; 2 2 2 2 2 2],...
         [all_fs(exp_idx==1 & cgp_idx==1)'; all_fs(exp_idx==1 & cgp_idx==2)'],...
         ".-", "Color", colCodes_Clu(3,:))
    plot([1 1 1; 2 2 2],...
         [all_fs(exp_idx==2 & cgp_idx==1)'; all_fs(exp_idx==2 & cgp_idx==2)'],...
         ".-", "Color" , colCodes_Clu(4,:))
    % aesthetics
    plotAesthetics(gca,1,font_size_large);
    set(gca, "XTickLabel", {'pre','post'});
    ylabel("80-20 decay time [ms]")
    ylim([0 320]); yLims = ylim;
    % statistics
    [H,P] = ttest(all_fs(cgp_idx==1), all_fs(cgp_idx==2)); % p = 0.0048
    text(mean(xlim), 280, '**',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
        'FontSize', font_stat_star, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', [.3 .3 .3])
    % n observations
    text(2.5, 270, "Ndnf > FS IN, n = 6",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "FontWeight", "bold", "Color", colCodes_Clu(3,:));
    text(2.5, 250, "NdnfNpy > FS IN, n = 3",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "FontWeight", "bold", "Color", colCodes_Clu(4,:));


%% plot traces
subplot(1,2,2)
    plot(ts_CGP_plot, traces_pre_mean, 'LineWidth', 1, 'Color', cols(1,:)); hold on
    plot(ts_CGP_plot, traces_post_mean, 'LineWidth', 1, 'Color', [181, 24, 47]./255);
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

