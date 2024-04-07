%% prepare workspace and load data
clear all; close all; clc;
cd('/main/');
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat');


%% --------------------- cross-cell FR comparison -------------------------
% format data
bp_FR_data = [resultTable_currentSteps.maxFR_Hz_evoked(~isnan(resultTable_currentSteps.maxFR_Hz_evoked) &...
                                                      resultTable_currentSteps.useBFanalysis == 1 &...
                                                      resultTable_currentSteps.isBF == 1 &...
                                                      resultTable_currentSteps.cellType == 1);... 
             resultTable_currentSteps.maxFR_Hz_evoked(~isnan(resultTable_currentSteps.maxFR_Hz_evoked) &...
                                                      resultTable_currentSteps.useBFanalysis == 1 &...
                                                      resultTable_currentSteps.isBF == 0 &...
                                                      resultTable_currentSteps.cellType == 1)];

bp_FR_group  = [ones(size(resultTable_currentSteps.maxFR_Hz_evoked(~isnan(resultTable_currentSteps.maxFR_Hz_evoked) &...
                                                      resultTable_currentSteps.useBFanalysis == 1 &...
                                                      resultTable_currentSteps.isBF == 1 &...
                                                      resultTable_currentSteps.cellType == 1)))*1;...
                ones(size(resultTable_currentSteps.maxFR_Hz_evoked(~isnan(resultTable_currentSteps.maxFR_Hz_evoked) &...
                                                      resultTable_currentSteps.useBFanalysis == 1 &...
                                                      resultTable_currentSteps.isBF == 0 &...
                                                      resultTable_currentSteps.cellType == 1)))*2];

% plot results
    bp = boxplot(bp_FR_data, bp_FR_group);
    % Aesthetics
    set(bp, "Color", "black", "LineWidth", 1.5);
	h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),blue_std,'FaceAlpha',.6);
    end
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
	ylabel("Max. evoked firing rate [Hz]"); 
    set(gca, "XTickLabel", ["Ndnf pF+","Ndnf pF-"]);
    plotAesthetics(gca, 1, 14);
    % statistics
    p_FR_ndnf = ranksum(bp_FR_data(bp_FR_group==1), bp_FR_data(bp_FR_group==2));
    line([1 2], [145 145], "Color", "black", "LineWidth", 2);
    text(1.5, 146, "***", "FontSize", 20, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "FontWeight", "bold");
	yLims = ylim;
    ylim([yLims(1) 150])
    % n observationns
    for iter = 1:2
        n = sum(bp_FR_group == iter);
        text(iter, 17, ['n=' num2str(n)], "FontSize", 12, "FontName", "Arial",...
             "VerticalAlignment", "top", "HorizontalAlignment", "center")
    end
    
    
