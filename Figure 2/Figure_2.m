%% prepare workspace and load data
close all; clc;
cd('/main/');
load('./analysis/misc/aesthetics.mat')

% indices for current step data
useTbl = [resultTable_currentSteps.cellNo(cleanDataInd)...
          resultTable_currentSteps.cellType(cleanDataInd)...
          resultTable_currentSteps.useBFanalysis(cleanDataInd)...
          resultTable_currentSteps.isBF(cleanDataInd)]; 

%% choose data for analysis
rawData = collectData; % include only data with valid sag analysis according to set criterion
clusterData = normalize(rawData, 1, 'range'); % normalize Ndnf data together with intersectional data for this analysis --> important for population comparisons

% separate into Ndnf and intersectonal data
% Ndnf
rawData_Ndnf = rawData(useTbl(:,2)==1, :);
normData_Ndnf = clusterData(useTbl(:,2)==1, :);
useTbl_Ndnf = useTbl(useTbl(:,2)==1, :);
% Ndnf/NPY
rawData_NdnfNPY = rawData(useTbl(:,2)==2, :);
normData_NdnfNPY = clusterData(useTbl(:,2)==2, :);
useTbl_NdnfNPY = useTbl(useTbl(:,2)==2, :);


% data dist pia
load("./analysis/Figure 2/distPia.mat");
% figure
F2_ndnf_Isteps = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ---------------------- recording paradigm ------------------------------
subplot(10,5,6);
    % baseline
    line([0 1062], [0 0], "Color", "Black"); hold on
    line([2062 4000], [0 0], "Color", "Black")
    % step
    for iter = -100:10:50
        line([1062 2062], [iter iter], "Color", "Black")
    end
    line([1062 1062], [-100 50], "Color", "Black")
    line([2062 2062], [-100 50], "Color", "Black")
    % step dotted
    for iter = 60:10:80
        line([1062 2062], [iter iter], "Color", "Black", "LineStyle", ":")
    end
    line([1062 1062], [50 80], "Color", "Black", "LineStyle", ":")
    line([2062 2062], [50 80], "Color", "Black", "LineStyle", ":")
    % aesthetics
    ylim([-110 150]); xlim([0 4000]);
    axis off
    % scale
    line([0 500], [-10 -10], "Color", "Black", "LineWidth", 1)
    text(250, -11, "500ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
  	line([2300 2300], [70 80], "Color", "Black", "LineWidth", 1)
    text(2360, 75, "\DeltaI = 10pA",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(2900, 0, "I = 0pA",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "left", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    text(2360, -100, "I_{min} = -100pA",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontWeight", "bold")

    pos = get(gca, "Position");
    set(gca, "Position", pos + [0 0 0 .045])

    
%% ---------------------- example recordings ------------------------------
cellNos = [143 74 67];
for iter = 1:3
    subplot(5,30, [7:10]-4+iter*4);
    cellNo = cellNos(iter);
    ts = 1:length(allData(cellNo).recordingData.rawData_mV(:,1)); ts = ts./50;
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
end
% scale
subplot(5,30, 7:10);
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
% legend
subplot(5,30, 15:18);
text(2300, 25, "suprathreshold",...
    "FontSize", font_size_large, "FontName", "Arial",...
    "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
    "FontWeight", "bold")
text(2300, 15, "subthreshold",...
    "FontSize", font_size_large, "FontName", "Arial", "Color", red,...
    "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
    "FontWeight", "bold")
text(2300, 5, "first",...
    "FontSize", font_size_large, "FontName", "Arial", "Color", "blue",...
    "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
    "FontWeight", "bold")


%% ----------------- ward based on normalized data ------------------------
Z = linkage(normData_Ndnf,'ward','euclidean');
ClustIDs = cluster(Z,'maxclust',2);
% calculate optimal order
D = pdist(normData_Ndnf);
leafOrder = optimalleaforder(Z,D);
% plot
subplot(5,30,[31 49]);
    d = dendrogram(Z, 0, 'Reorder',leafOrder,'ColorThreshold',4);  
    xlabel('cell ID');
%    set(d, "Color", "black")
    set(gca, 'FontSize', font_size_large,...
             'XTickLabel', '',...
             "TickLength", [0 0]);
    set(d, 'LineWidth', 1);
    plotAesthetics(gca,1,font_size_large)
    % labels &  n observations
    xLims = xlim; yLims = ylim;
    text(xLims(1), yLims(2)*0.9, ' Cluster 2',...
        'FontSize', font_size_large, "FontName", "Arial",...
        'HorizontalAlignment', 'left', "VerticalAlignment", "top",...
        "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "FontWeight", "bold")
    text(xLims(2), yLims(2)*0.9, 'Cluster 1 ',...
        'FontSize', font_size_large, "FontName", "Arial",... 
        'HorizontalAlignment', 'right', "VerticalAlignment", "top",...
        "Color", colCodes_Clu(1,:),...
        "FontWeight", "bold")
    text(xLims(1), yLims(2)*0.83, ['    n=', num2str(sum(ClustIDs==2))],...
        'FontSize', font_size_large, "FontName", "Arial",...
        'HorizontalAlignment', 'left', "VerticalAlignment", "top",...
        "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "FontWeight", "bold")
    text(xLims(2), yLims(2)*0.83, ['n=', num2str(sum(ClustIDs==1)) '    '],...
        'FontSize', font_size_large, "FontName", "Arial",...
        'HorizontalAlignment', 'right', "VerticalAlignment", "top",...
        "Color", colCodes_Clu(1,:),...
        "FontWeight", "bold")
	text(sum(ClustIDs==2)+.5, yLims(2)*0.9, "All Ndnf INs",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", "FontSize", font_size_large,...
        "FontWeight", "bold");
    text(sum(ClustIDs==2)+.5, yLims(2)*0.83, ['n=' num2str(length(ClustIDs))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontName", "Arial", "FontSize", font_size_large,...
        "FontWeight", "bold");
    h=findobj('Color','red'); set(h, 'Color', colCodes_Clu(2,:));
    h=findobj('Color','cyan'); set(h, 'Color', colCodes_Clu(1,:));
    ylabel("Euclidean distance [AU]"); xlabel('Cell ID');
    
% add ward results to resultTable
ClustIDs_ResultTable = NaN(size(resultTable_currentSteps,1),1);
for iter = 1:size(useTbl_Ndnf,1)
    cellNos = useTbl_Ndnf(iter,1);
    ClustIDs_ResultTable(cellNos) = ClustIDs(iter);
end
resultTable_currentSteps.ClustIDs = ClustIDs_ResultTable;

% position
ward_pos = get(gca, "Position");
set(gca, "Position", [ward_pos(1) ward_pos(2) ward_pos(3)-0.01 ward_pos(4)])


%% ------------------------------- PCA ------------------------------------
% CALCULATE BASED ON NDNF CELLS
[coeff,score,latent,tsquared,explained,mu] = pca(normData_Ndnf);

sz1 = ones(sum(ClustIDs == 1),1)*marker_size_pca;
sz2 = ones(sum(ClustIDs == 2),1)*marker_size_pca;

% PC1 vs PC2
subplot(5,30,61:65);
    % single cells
    scatter(score(ClustIDs == 1,1), score(ClustIDs == 1,2), sz1,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,1), score(ClustIDs == 2,2), sz2,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    % population centroids
    plot(mean(score(ClustIDs == 1,1)), mean(score(ClustIDs == 1,2)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,1)), mean(score(ClustIDs == 2,2)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 1'); ylabel('PC 2'); 
    xlim([-1 1]); ylim([-1 1]);
    plotAesthetics(gca, 1, font_size_large); 
	% position
    pca_pos = get(gca, "Position");
    set(gca, "Position", [pca_pos(1) pca_pos(2)+0.03 pca_pos(3) pca_pos(4)-0.01])
    
% PC1 vs PC3
subplot(5,30,68:72);
    % single cells
    scatter(score(ClustIDs == 1,1), score(ClustIDs == 1,3), sz1,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,1), score(ClustIDs == 2,3), sz2,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    % population centroids
    plot(mean(score(ClustIDs == 1,1)), mean(score(ClustIDs == 1,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,1)), mean(score(ClustIDs == 2,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 1'); ylabel('PC 3');
    xlim([-1 1]); ylim([-1 1]);
    plotAesthetics(gca, 1, font_size_large);
	% position
    pca_pos = get(gca, "Position");
    set(gca, "Position", [pca_pos(1) pca_pos(2)+0.03 pca_pos(3) pca_pos(4)-0.01])
    
% PC2 vs PC3    
subplot(5,30,75:79);
    % single cells
    scatter(score(ClustIDs == 1,2), score(ClustIDs == 1,3), sz1,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(1,:), 'MarkerFaceColor', colCodes_Clu(1,:)); hold on
    scatter(score(ClustIDs == 2,2), score(ClustIDs == 2,3), sz2,...
            'Marker', 'o', 'MarkerEdgeColor', colCodes_Clu(2,:), 'MarkerFaceColor', colCodes_Clu(2,:))
    % population centroids
    plot(mean(score(ClustIDs == 1,2)), mean(score(ClustIDs == 1,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(1,:), 'LineWidth', 1)
    plot(mean(score(ClustIDs == 2,2)), mean(score(ClustIDs == 2,3)),...
         'h', 'Color', 'black', 'MarkerSize', 7, 'MarkerFaceColor', colCodes_Clu(2,:), 'LineWidth', 1)
    % aesthetics
    xlabel('PC 2'); ylabel('PC 3'); 
    xlim([-1 1]); ylim([-1 1]);
    plotAesthetics(gca, 1, font_size_large); 
    % legend
    text(1, -0.6, ['Clu1, n=', num2str(sum(ClustIDs==1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(1, -0.75, ['Clu2, n=', num2str(sum(ClustIDs==2))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(1, -0.9, "Centroids",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", 'black',...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	plot(0, -.9, 'h', 'Color', 'black',...
        'MarkerSize', 5, 'MarkerFaceColor', 'white', 'LineWidth', 1)
  
    % position
    pca_pos = get(gca, "Position");
    set(gca, "Position", [pca_pos(1) pca_pos(2)+0.03 pca_pos(3) pca_pos(4)-0.01])
    
    
%% -------------------------- spike onset CDF -----------------------------
subplot(15,30,[111 145]);
    h_OnsetCLu1 = histcounts(rawData_Ndnf(ClustIDs == 1,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetCLu1, "Color", colCodes_Clu(1,:), "LineWidth", 1); hold on
    h_OnsetClu2 = histcounts(rawData_Ndnf(ClustIDs == 2,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetClu2, "Color", colCodes_Clu(2,:), "LineWidth", 1);
    h_OnsetAll = histcounts(rawData_Ndnf(:,1),"BinWidth", 1, "Normalization", "cdf", "BinLimits", [0 957]);
    stairs(h_OnsetAll, "Color", colCodes_Clu(3,:), "LineWidth", 1);
    plotAesthetics(gca, 1, font_size_large);
    xlabel("First Spike Onset [ms]"); ylabel("Cum. probability")
    text(950, 0.3, ['all, n=', num2str(length(ClustIDs))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(3,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(950, 0.2, ['Clu1, n=', num2str(sum(ClustIDs==1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
	text(950, 0.1, ['Clu2, n=', num2str(sum(ClustIDs==2))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")

% position
onset_pos = get(gca, "Position");
set(gca, "Position", [onset_pos(1) onset_pos(2) onset_pos(3)-0.01 onset_pos(4)])

   
%% ---------------------------- distance to pia ---------------------------
clear fNames
clustID_pia = NaN(size(distPia,1),1);
% extract recording file names
for iter = 1:length(allData)
    fNames(iter) = string(allData(iter).recordingFileName(1:end-4));
end
% remove those not used for clustering
fNames = fNames(cleanDataInd)'; % clean data ind from cross cells collection script
fNames = fNames(resultTable_currentSteps.cellType(cleanDataInd) ==1); % NPY+ cells
% find matches
for iter = 1:size(distPia,1)
    % extract and repair file names
    pia_name = char(table2array(distPia(iter,1)));
    pia_name = string([pia_name(1:end-8) pia_name(end-3:end)]);
    for jter = 1:length(fNames)
        if pia_name == fNames(jter)
            % assign cluster ID
            clustID_pia(iter) = ClustIDs(jter);
        end
    end
end
% data for plotting
dPia = table2array(distPia(~isnan(clustID_pia),3));
clustID_pia_2 = clustID_pia(~isnan(clustID_pia));
% plot data
subplot(15,30,[111 145]+5);
    bp_pia = boxplot(dPia, clustID_pia_2);
    set(bp_pia, "Color", "black"); 
    h = findobj(gca,'Tag','Box');
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_Clu(1,:),'FaceAlpha',.8);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_Clu(2,:),'FaceAlpha',.8);
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel("Distance to pia [µm]");
    set(gca, "XTickLabel", ["Clu1", "Clu2"]);
    plotAesthetics(gca, 1, font_size_large);
    % statistics
    p_pia = ranksum(dPia(clustID_pia_2==1), dPia(clustID_pia_2==2)); % p = 0.4802
    yLims = ylim;
    text(1.5, yLims(2), "n.s.",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontWeight", "bold")
    % legend
    yLims = ylim;
    text(1, yLims(1)+0.15*-diff(yLims), ['n=', num2str(sum(clustID_pia_2==1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", 'black',...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(2, yLims(1)+0.15*-diff(yLims), ['n=', num2str(sum(clustID_pia_2==2))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", 'black',...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    
% position
pia_pos = get(gca, "Position");
set(gca, "Position", [pia_pos(1)+0.04 pia_pos(2) pia_pos(3)-0.04 pia_pos(4)])


%% ------------------------ avg clu1 clu2 sag/hump ------------------------
% variables
cluInds = resultTable_currentSteps.ClustIDs(~isnan(resultTable_currentSteps.ClustIDs));
nClust = sum(~isnan(resultTable_currentSteps.ClustIDs));
sagSweep = NaN(nClust, 200000);
subThreshSweep = NaN(nClust, 200000);
% extact data
counter = 1;
for iter = 1:size(resultTable_currentSteps,1)
    if ~isnan(resultTable_currentSteps.ClustIDs(iter))
        sagSweep_ind = allData(iter).analysis.sag.base75SweepInd;
        sagSweep(counter, :) = allData(iter).recordingData.rawData_mV(:,sagSweep_ind);
        subThreshSweep_ind = allData(iter).analysis.supraThresholdSweep - 1;
        subThreshSweep(counter, :) = allData(iter).recordingData.rawData_mV(:,subThreshSweep_ind);
        counter = counter + 1;
    end
end
% average +- SEM
[mean_sag_clu1, std_sag_clu1, ~, ~] = statistics(sagSweep(cluInds == 1,:), 1);
[mean_sub_clu1, std_sub_clu1, ~, ~] = statistics(subThreshSweep(cluInds == 1,:), 1);
[mean_sag_clu2, std_sag_clu2, ~, ~] = statistics(sagSweep(cluInds == 2,:), 1);
[mean_sub_clu2, std_sub_clu2, ~, ~] = statistics(subThreshSweep(cluInds == 2,:), 1);

% plot results
ts = (1:length(std_sag_clu1))./1000*20;
subplot(15,30,[111 145]+60); hold on 
    % clu1 hump
    plot(ts, mean_sub_clu1, "Color", colCodes_Clu(1,:), "LineWidth", 1)
	ylim([-72 -30]); axis off
	% scale
    line([3000 4000],[-60 -60],"Color","black","LineWidth",width_scalebar)
    line([3000 3000],[-60 -50],"Color","black","LineWidth",width_scalebar)
    text(3750, -60.5, "1000ms ",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    text(3000, -55, " 10mV",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "LEFT", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    % legend
    text(1500, -77, "Clu1",...
        "FontSize", font_size_large, "FontName", "Arial","Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontWeight", "bold")
    text(1500, -77, ['n=', num2str(sum(ClustIDs==1))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(1,:),...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    
subplot(15,30,[111 145]+120); hold on
    % clu1 sag
	plot(ts, mean_sag_clu1, "Color", colCodes_Clu(1,:), "LineWidth", 1)
	ylim([-85 -65]); axis off
	% scale
    line([3000 4000],[-76 -76],"Color","black","LineWidth",width_scalebar)
    line([3000 3000],[-76 -71],"Color","black","LineWidth",width_scalebar)
    text(3750, -76.5, "1000ms ",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    text(3000, -73.5, " 5mV",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "LEFT", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    % arrow
	x = 1085;
	pgon = polyshape([x x-150 x+150],[-78.5 -80 -80]); 
    plot(pgon, "EdgeColor", "black", "FaceColor", "black",...
               "FaceAlpha", 1, "EdgeAlpha", 1); 


subplot(15,30,[111 145]+65); hold on 
    % ramp clu2
    plot(ts, mean_sub_clu2, "Color", colCodes_Clu(2,:), "LineWidth", 1)
	ylim([-72 -30]); axis off
    % legend
    text(1500, -77, "Clu2",...
        "FontSize", font_size_large, "FontName", "Arial","Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontWeight", "bold")
    text(1500, -77, ['n=', num2str(sum(ClustIDs==2))],...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_Clu(2,:)-[.2 0 .2],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
	% arrow
	x = 1535;
	pgon = polyshape([x x-150 x+150],[-45.5 -42 -42]); 
    plot(pgon, "EdgeColor", "black", "FaceColor", "black",...
               "FaceAlpha", 1, "EdgeAlpha", 1);

subplot(15,30,[111 145]+125); hold on
    % sag clu2
	plot(ts, mean_sag_clu2, "Color", colCodes_Clu(2,:), "LineWidth", 1)
	ylim([-85 -65]); axis off

    
    
%% save figure
set(F2_ndnf_Isteps,'renderer','Painters')
%print('F2_intrinsicProperties.pdf','-dpdf')
