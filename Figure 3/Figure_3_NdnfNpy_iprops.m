%% prepare workspace and load data
% NOTE: Run Figure 2 first
close all; clc;
cd("/main/");
load("./analysis/misc/aesthetics.mat")

% intersect PCA score (NOTE: run Figure 2 first, clear everything else except allData & resultTable)
score_NdnfNPY = (normData_NdnfNPY-mu)  * coeff;

% figure
F3_ndnfnpy_ephys = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% -------------------- example recording cell 1 ----------------------------
cellNo = 121;
ts = 1:length(allData(cellNo).recordingData.rawData_mV(:,1)); ts = ts./50;
subplot(7,7,1)
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,1), "Color", "blue");hold on
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,allData(cellNo).analysis.supraThresholdSweep), "Color", "Black");   
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,allData(cellNo).analysis.supraThresholdSweep-1), "Color", red);
    text(1100,-110,"-100pA", "Color", "blue",...
        "HorizontalAlignment", "center", "VerticalAlignment","bottom",...
        "FontSize", font_size_small, "FontName", "Arial")
    text(1550,-110,[num2str(allData(cellNo).recordingData.Isteps_pA(allData(cellNo).analysis.supraThresholdSweep-1)) 'pA'],...
        "Color", red, "HorizontalAlignment", "center", "VerticalAlignment","bottom",....
        "FontSize", font_size_small, "FontName", "Arial")
    text(2000,-110,[num2str(allData(cellNo).recordingData.Isteps_pA(allData(cellNo).analysis.supraThresholdSweep)) 'pA'],...
        "Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
	text(1355,-110,"/","Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
    text(1780,-110,"/","Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
    ylim([-100 35]); axis off
    xlim([900 2200])
    % scale
    line([1300 1300],[0 20],"Color", "Black", "LineWidth", width_scalebar)
    line([1300 1600],[0 0],"Color", "Black", "LineWidth", width_scalebar)
    text(1450, -1, "300ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    text(1290, 10, "20mV ",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
% position
ephys_pos = get(gca, "Position");
set(gca, "Position", [ephys_pos(1) ephys_pos(2)-0.01 ephys_pos(3) ephys_pos(4)+0.03])


%% -------------------- example recording cell 2 ----------------------------
cellNo = 179;
ts = 1:length(allData(cellNo).recordingData.rawData_mV(:,1)); ts = ts./50;
subplot(7,7,2)
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,1), "Color", "blue");hold on
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,allData(cellNo).analysis.supraThresholdSweep), "Color", "Black");   
    plot(ts,allData(cellNo).recordingData.rawData_mV(:,allData(cellNo).analysis.supraThresholdSweep-1), "Color", red);
	text(1100,-110,"-100pA", "Color", "blue", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
    text(1550,-110,[num2str(allData(cellNo).recordingData.Isteps_pA(allData(cellNo).analysis.supraThresholdSweep-1)) 'pA'],...
        "Color", red, "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
    text(2000,-110,[num2str(allData(cellNo).recordingData.Isteps_pA(allData(cellNo).analysis.supraThresholdSweep)) 'pA'],...
        "Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
	text(1355,-110,"/","Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
	text(1780,-110,"/","Color", "black", "HorizontalAlignment", "center", "VerticalAlignment","bottom", "FontSize", font_size_small, "FontName", "Arial")
    ylim([-100 35]); axis off
    xlim([900 2200])
% position
ephys_pos = get(gca, "Position");
set(gca, "Position", [ephys_pos(1) ephys_pos(2)-0.01 ephys_pos(3) ephys_pos(4)+0.03])

    
%%  -------------------------- distance -----------------------------------
% calculate centroids of Ward's method clustering
centroid_Clu1 = mean(normData_Ndnf(ClustIDs==1,:),1);
centroid_Clu2 = mean(normData_Ndnf(ClustIDs==2,:),1);
centroid_ndnf = mean(normData_Ndnf,1);
% find Ndnf Npy cells suitable for analysis & centroid
centroid_ndnf_NPY = mean(normData_NdnfNPY,1);
% find distances of intersectional cells to centroids 
dists = NaN(size(normData_NdnfNPY,1),3);
for iter = 1:size(normData_NdnfNPY,1)
    dists(iter, 1) = norm(centroid_Clu1 - normData_NdnfNPY(iter,:));
    dists(iter, 2) = norm(centroid_Clu2 - normData_NdnfNPY(iter,:));
    dists(iter, 3) = dists(iter, 2)/dists(iter, 1);
end
subplot(7,7,[8 9])
    % histogram
    hDist = histogram(log2(dists(:, 3)), 'BinWidth', 0.2); hold on
    % single data
    scatter(log2(dists(:,3)), zeros(size(dists(:,3)))+6, "MarkerEdgeColor","black", "SizeData", .5);
    line([mean(log2(dists(:, 3))) mean(log2(dists(:, 3)))], [0 6.9], 'Color', red, 'LineStyle', ':', "LineWidth", 1);
    % aesthetics
    set(hDist, "FaceColor", colCodes_Clu(4,:), "FaceAlpha", .8)
	xSize = max(abs(xlim)); xlim([-xSize xSize]); ylim([0 8]);
    ylabel('Number of cells'); xlabel('log_2(dist. Clu2 / dist. Clu1)')
    plotAesthetics(gca, 1, font_size_large);
    ylim([0 7]);
    % statistics
    [h,p] = kstest(log2(dists(:, 3))); % p = 0.1132       
    text(0, 7, "n.s.",...
         "FontSize", font_size_large, "FontName", "Arial",...
         "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
         "FontWeight", "bold");
    % n observations
    text(0, -2.5, ['n=' num2str(size(score_NdnfNPY,1))] ,...
         "FontName", "Arial", "FontSize", font_size_large, "Color", 'black',...
         "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
         "FontWeight", "bold"); 


%% -------------------------- spike onset CDF -----------------------------
subplot(7,7,[11 12])
    h_OnsetClu1 = histcounts(rawData_Ndnf(ClustIDs == 1,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetClu1, "Color", colCodes_Clu(1,:), "LineWidth", 1); hold on
    h_OnsetClu2 = histcounts(rawData_Ndnf(ClustIDs == 2,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetClu2, "Color", colCodes_Clu(2,:), "LineWidth", 1);
    h_OnsetAll = histcounts(rawData_Ndnf(:,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetAll, "Color", colCodes_Clu(3,:), "LineWidth", 1);
    h_OnsetIntersect = histcounts(rawData_NdnfNPY(:,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetIntersect, "Color", colCodes_Clu(4,:), "LineWidth", 1);
    plotAesthetics(gca, 1, font_size_large); ylim([-0.05 1.05])
    xlabel("First spike onset [ms]"); ylabel("Cum. pprobability")
    text(950, 0.4, ['all Ndnf, n=', num2str(length(ClustIDs))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(3,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(950, 0.3, ['Clu1, n=', num2str(sum(ClustIDs==1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(950, 0.2, ['Clu2, n=', num2str(sum(ClustIDs==2))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(950, 0.1, ['NdnfNpy, n=', num2str(size(rawData_NdnfNPY,1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(4,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")

    
%% ------------------------- boxplot distances ----------------------------
subplot(7,7,14)
    bp = boxplot(dists(:,1:2));
    set(bp, "Color", "black"); hold on
    plot(repmat([1 2]', 1, size(dists(:,1:2),1)), dists(:,1:2)', "-", "Color", colCodes_Clu(4,:))
    h = findobj(gca,'Tag','Box');
	patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_Clu(2,:),'FaceAlpha',.6);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_Clu(1,:),'FaceAlpha',.6);
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    ylabel("Dist. to centroid [AU]");
  	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    set(gca, "XTickLabel", ["Clu1", "Clu2"]); 
    % n observations
    text(1.5,.13,['n=' num2str(size(score_NdnfNPY,1))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",....
        "FontSize", font_size_large, "FontName", "Arial",...
        "FontWeight", "bold");
    % statistics
    p = signrank(dists(:,1), dists(:,2));
    text(1.5, 1.4, "n.s.",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "FontWeight", "bold")

    
%% --------------------------- PCA plots ------------------------------- %%
% ---------------------------- PCA 1 vs 2 ------------------------------- %
subplot(7,6,[13 20])
% plot single data
    scatter(score(ClustIDs == 1,1), score(ClustIDs == 1,2), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,1), score(ClustIDs == 2,2), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    scatter(score_NdnfNPY(:,1), score_NdnfNPY(:,2), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(4,:), 'MarkerFaceColor', colCodes_Clu(4,:))
    % population centroids
    plot(mean(score(:,1)), mean(score(:,2)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(3,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 1,1)), mean(score(ClustIDs == 1,2)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,1)), mean(score(ClustIDs == 2,2)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    plot(mean(score_NdnfNPY(:,1)), mean(score_NdnfNPY(:,1)), 'h',...
         'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(4,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 1'); ylabel('PC 2'); 
    xlim([-1.1 1.1]); ylim([-1.1 1.1])
    plotAesthetics(gca, 1,font_size_large);
    axis square

% ----------------------------- PCA 1 vs 3 ------------------------------ %
subplot(7,6,[13 20]+2)
% plot single data
    scatter(score(ClustIDs == 1,1), score(ClustIDs == 1,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,1), score(ClustIDs == 2,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    scatter(score_NdnfNPY(:,1), score_NdnfNPY(:,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(4,:), 'MarkerFaceColor', colCodes_Clu(4,:))
    % population centroids
    plot(mean(score(:,1)), mean(score(:,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(3,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 1,1)), mean(score(ClustIDs == 1,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,1)), mean(score(ClustIDs == 2,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    plot(mean(score_NdnfNPY(:,3)), mean(score_NdnfNPY(:,3)), 'h',...
         'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(4,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 1'); ylabel('PC 3'); 
    xlim([-1.1 1.1]); ylim([-1.1 1.1])
    plotAesthetics(gca, 1,font_size_large);
    axis square
    
% ---------------------------- PCA 2 vs 3 ------------------------------- %
subplot(7,6,[13 20]+4)
% plot single data
    scatter(score(ClustIDs == 1,2), score(ClustIDs == 1,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,2), score(ClustIDs == 2,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    scatter(score_NdnfNPY(:,2), score_NdnfNPY(:,3), marker_size_pca,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(4,:), 'MarkerFaceColor', colCodes_Clu(4,:))
    % population centroids
    plot(mean(score(:,2)), mean(score(:,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(3,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 1,2)), mean(score(ClustIDs == 1,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,2)), mean(score(ClustIDs == 2,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    plot(mean(score_NdnfNPY(:,2)), mean(score_NdnfNPY(:,3)), 'h',...
         'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(4,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 2'); ylabel('PC 3'); 
    xlim([-1.1 1.1]); ylim([-1.1 1.1])
    plotAesthetics(gca, 1,font_size_large);
    axis square
	% legend
	text(-1, 1, "all Ndnf",...
        "FontName", "Arial", "FontSize", font_size_large, "Color", colCodes_Clu(3,:),...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(-1, .9, ['Cluster 1, n='  num2str(sum(ClustIDs == 1))],...
        "FontName", "Arial", "FontSize", font_size_large, "Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(-1, .8, ['Cluster 2, n='  num2str(sum(ClustIDs == 2))],...
        "FontName", "Arial", "FontSize", font_size_large, "Color", colCodes_Clu(2,:),...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(-1, .7, ['Ndnf/Npy, n=' num2str(size(score_NdnfNPY,1))] ,...
        "FontName", "Arial", "FontSize", font_size_large, "Color", colCodes_Clu(4,:),...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    plot(-1, -1, 'h', 'Color', 'black', 'MarkerSize', 6, 'MarkerFaceColor', "white", 'LineWidth', 1)
    text(-.93, -1, "Centroids",...
        "FontName", "Arial", "FontSize", font_size_large, "Color", "black",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    

%% save figure
set(gcf,'renderer','Painters')
%print('F3_ndnfnpy_ephys.pdf','-dpdf')