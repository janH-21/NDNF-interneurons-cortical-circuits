%% prepare workspace
clear all; close all; clc
cd("/main/")
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat');

%% CURRENT STEPS: proportion of BF+ cells 
% read data
ctrl_neg = sum(resultTable_currentSteps.cellType == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 0);
ctrl_pos = sum(resultTable_currentSteps.cellType == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 1);

% currently hard coded, replace with direct extraction from data
syn_step_neg = 8;
syn_step_pos = 16;
ctrl_28Hz_neg = 10;
ctrl_2Hz_pos = 18;
syn_28Hz_neg = 5;
syn_28Hz_pos = 10;

% format data
occur_step = [ctrl_pos ctrl_neg;
              syn_step_pos syn_step_neg];
occur_28Hz = [ctrl_2Hz_pos ctrl_28Hz_neg;
              syn_28Hz_pos syn_28Hz_neg];
        
occur_step_perc = (occur_step(:,1)./(occur_step(:,1) + occur_step(:,2)))*100;
occur_step_perc(:,2) = 100 - occur_step_perc;

occur_28Hz_perc = (occur_28Hz(:,1)./(occur_28Hz(:,1) + occur_28Hz(:,2)))*100;
occur_28Hz_perc(:,2) = 100 - occur_28Hz_perc;

% plot results steps
subplot(1,2,1)
b = bar(occur_step_perc,"BarLayout", "stacked");
    set(b(1), "FaceColor", colCodes_bar(1,:), "FaceAlpha", .6);
	set(b(2), "FaceColor", colCodes_bar(2,:), "FaceAlpha", .6);
    % aesthetics
    ylabel("Percent of cells");
     set(gca, "XTickLabel", ["Ctrl", "Syn. blockers"],...
              "YTick", [0 50 100]);
    xtickangle(35);
    plotAesthetics(gca,1,14);
    % legend
    text(0.8, 140, 'pF+',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', 14, 'FontWeight', 'bold',...
        'FontName', 'Arial', 'Color', colCodes_bar(1,:));
    text(0.8, 150, 'pF-',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', 14, "FontWeight", "bold",...
        'FontName', 'Arial', 'Color', colCodes_bar(2,:));
    xlim([0.7 2.5]); ylim([0 107])  
    title("Persistent Firing Steps")
    % fisher exact post hoc + bonferroni
    [h, p, stats] = fishertest(array2table(occur_step));
    subtitle(['p=' num2str(p)])
    % add n observations
    for iter = 1:size(occur_step,1)
        n = sum(occur_step(iter,1:2));
        text(iter, 102, ['n=' num2str(n)],...
            "FontSize", 14, "FontName", "Arial",...
            "VerticalAlignment", "bottom", "HorizontalAlignment", "center")
    end

% plot results 28Hz
subplot(1,2,2)
b = bar(occur_28Hz_perc,"BarLayout", "stacked");
    set(b(1), "FaceColor", colCodes_bar(1,:), "FaceAlpha", .6);
	set(b(2), "FaceColor", colCodes_bar(2,:), "FaceAlpha", .6);
    % aesthetics
    ylabel("Percent of cells");
    set(gca, "XTickLabel", ["Ctrl", "Syn. blockers"],...
              "YTick", [0 50 100]);
    xtickangle(35);
    plotAesthetics(gca,1,14);
    % legend
    text(0.8, 140, 'pF+',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', 14, 'FontWeight', 'bold',...
        'FontName', 'Arial', 'Color', colCodes_bar(1,:));
    text(0.8, 150, 'pF-',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', 14, "FontWeight", "bold",...
        'FontName', 'Arial', 'Color', colCodes_bar(2,:));
    xlim([0.7 2.5]); ylim([0 107])  
    title("Persistent Firing 28Hz")
    % fisher exact post hoc + bonferroni
    [h, p, stats] = fishertest(array2table(occur_28Hz));
    subtitle(['p=' num2str(p)])
    % add n observations
    for iter = 1:size(occur_step,1)
        n = sum(occur_28Hz(iter,1:2));
        text(iter, 102, ['n=' num2str(n)],...
            "FontSize", 14, "FontName", "Arial",...
            "VerticalAlignment", "bottom", "HorizontalAlignment", "center")
    end