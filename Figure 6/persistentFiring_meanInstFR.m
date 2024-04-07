%% prepare workspace and load data
% run cross-cell current step master script first
close all; clc;
cd("/main/")
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat')

% avg. max inter-step ifreq box
counter = 1;
ifreqs_maxBF = {};
mean_ifreq_max_BF = [];
for iter = 1:length(allData)
    % choose only NDNF cell with BF and eligible for analysis
    if resultTable_currentSteps.cellType(iter) == 1 &&...
       resultTable_currentSteps.isBF(iter) == 1 &&...
       resultTable_currentSteps.useBFanalysis(iter) == 1
        % find episode with max. # of spikes
        [~, maxBF_ind] = max(allData(iter).analysis.nBarrageSpikes);
        % get ifreq for max BF episode
        maxBF_1_Ts = allData(iter).analysis.APstats.threshTs{maxBF_ind};
        maxBF_2_Ts = allData(iter).analysis.APstats.threshTs{maxBF_ind+1};
        maxBF_1_ind = maxBF_1_Ts > 2062;
        maxBF_2_ind = maxBF_2_Ts < 1062;
        maxBF_all_Ts = [maxBF_1_Ts(maxBF_1_ind)...
                        maxBF_2_Ts(maxBF_2_ind)+4000];
        ifreqs_maxBF{counter} = 1000./diff(maxBF_all_Ts); % ifreq in ms
        mean_ifreq_max_BF(counter) = mean(ifreqs_maxBF{counter});
        counter = counter + 1;
    end   
end

% plot mean ifreq of inter-step with max. AP#; mean = 51.2245 Hz; median = 50.1405 Hz
bp = boxplot(mean_ifreq_max_BF);
    set(bp, "Color", "black");
    plotAesthetics(gca, 1, font_size_large);
    % aesthetics
    ylabel("Mean inst. FR [Hz]"); xlabel("Inter-step with max. AP#")
    h = findobj(gca,'Tag','Box');
    patch(get(h,'XData'),get(h,'YData'),blue_std,'FaceAlpha',.8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", []); 
    % n observations
    text(1,-15,['n=' num2str(length(mean_ifreq_max_BF))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",....
        "FontSize", 14, "FontName", "Arial");
