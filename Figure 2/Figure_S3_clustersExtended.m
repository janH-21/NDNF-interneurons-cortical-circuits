%% prepare workspace and load data
close all; clc;
cd("/main/")
load('./analysis/misc/aesthetics.mat');
load('./analysis/misc/varNames.mat')
% working directory
% figure
FS3_clusters = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% --------------------- WARD CLUSTER EVALUATION --------------------------
evalDB = evalclusters(normData_Ndnf,'linkage','DaviesBouldin','KList', 1:size(normData_Ndnf,1));
evalSilhouette = evalclusters(normData_Ndnf,'linkage', 'silhouette','KList', 1:size(normData_Ndnf,1));
evalGap  = evalclusters(normData_Ndnf,'linkage','gap','KList', 1:size(normData_Ndnf,1));
evalCH = evalclusters(normData_Ndnf,'linkage','CalinskiHarabasz','KList', 1:size(normData_Ndnf,1));


%% plot criteria 
%% DB
subplot(7,6,1)
    plt_DB_1 = plot(evalDB);
	set(plt_DB_1, "Color", blue_std, "LineWidth", 1, "Marker", ".")
    % plot(1:5, evalDB.CriterionValues(1:5), ".-", "LineWidth", 1);
    ylabel("Values"); xlabel("Cluster count");
    xlim([0.5 5.5]);
    text(1.9,2.1,"Local optimum", "HorizontalAlignment", "right", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "Rotation", 90,...
        "FontWeight", "bold")
    plotAesthetics(gca, 1, font_size_large);
    title("Davis-Bouldin criterion")
    line([2 2], ylim, "Color", "red", "LineStyle", ":")
subplot(7,6,[2 3])
    plt_DB_2 = plot(evalDB); 
    set(plt_DB_2, "Color", blue_std, "LineWidth", 1, "Marker", "none")
    % plot(1:5, evalDB.CriterionValues(1:5), ".-", "LineWidth", 1);
    ylabel("Values"); xlabel("Cluster count");
    xlim([0 size(normData_Ndnf,1)]+.5);
    text(145,2.1,"Global optimum", "HorizontalAlignment", "right", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "Rotation", 90,...
        "FontWeight", "bold")
    plotAesthetics(gca, 1, font_size_large);
    line([evalDB.OptimalK evalDB.OptimalK], ylim, "Color", "red", "LineStyle", "-")
    
    
%% Silhouette
subplot(7,6,4)
    plt_sh_1 = plot(evalSilhouette);
    set(plt_sh_1, "Color", blue_std, "LineWidth", 1, "Marker", ".")
    ylabel("Values"); xlabel("Cluster count");
    xlim([0.5 5.5]); 
    plotAesthetics(gca, 1, font_size_large);
	title("Silhouette criterion")
    line([2 2], [0.2 0.35], "Color", "red", "LineStyle", ":")
subplot(7,6,[5 6])
    plt_sh_2 = plot(evalSilhouette);
    set(plt_sh_2, "Color", blue_std, "LineWidth", 1, "Marker", "none")
    ylabel("Values"); xlabel("Cluster count");
    xlim([0 size(normData_Ndnf,1)]+0.5); 
    plotAesthetics(gca, 1, font_size_large);
    line([evalSilhouette.OptimalK evalSilhouette.OptimalK], ylim, "Color", "red", "LineStyle", "-")
    
%% gap
subplot(7,6,[7])
    plt_gap_1 = plot(evalGap);  %plot(1:5, evalGap.CriterionValues(1:5), ".-", "LineWidth", 1);
    set(plt_gap_1, "Color", blue_std, "LineWidth", 1, "Marker", ".")
    ylabel("Values"); xlabel("Cluster count");
    xlim([0.5 5.5]); 
    plotAesthetics(gca, 1, font_size_large);
    title("Gap value criterion"); 
	line([2 2], ylim, "Color", "red", "LineStyle", ":")
subplot(7,6,[8 9])
    plot(1:size(normData_Ndnf,1), evalGap.CriterionValues, "-", "LineWidth", 1);
    ylabel("Values"); xlabel("Cluster count");
    xlim([0 size(normData_Ndnf,1)]+0.5); 
    plotAesthetics(gca, 1, font_size_large);
	line([evalGap.OptimalK evalGap.OptimalK], ylim, "Color", "red", "LineStyle", "-")
        
%% Calinski-Harabasz criterion
subplot(7,6,10)
    plt_sh_1 = plot(evalCH);
    set(plt_sh_1, "Color", blue_std, "LineWidth", 1, "Marker", ".")
    ylabel("Values"); xlabel("Cluster count");
    xlim([0.5 5.5]); 
    plotAesthetics(gca, 1, font_size_large);
	title("Calinski-Harabasz criterion")
    line([2 2], [0 40], "Color", "red", "LineStyle", ":")
subplot(7,6,[11 12])
    plt_sh_2 = plot(evalCH);
    set(plt_sh_2, "Color", blue_std, "LineWidth", 1, "Marker", "none")
    ylabel("Values"); xlabel("Cluster count");
    xlim([0 size(normData_Ndnf,1)]+0.5); 
    plotAesthetics(gca, 1, font_size_large);
    line([evalCH.OptimalK evalCH.OptimalK], ylim, "Color", "red", "LineStyle", "-")
    
    
%% -------------------------PCA var explained -----------------------------
% Panel H+I: PCA % Variance Explained
% cumulative variance explained
cumExplained = NaN(size(explained));
for iter = 1:length(explained)
    cumExplained(iter) = sum(explained(1:iter));
end
subplot(7, 4, [10 11])
    bpCumExplained = bar(cumExplained, "FaceAlpha", .8); 
    set(gca,'YTick', [0:20:100])
    ylabel(["Cum. variance" "explained [%]"]); xlabel('Principle components')
    plotAesthetics(gca, 1, font_size_large);
	cumVar = gca;
    pos_cumVar = cumVar.Position;
    set(cumVar, "Position", pos_cumVar + [.04 0 -.04 0])
    
    
%% ---------------------- PCA feature correlation -------------------------
n_cols = size(score, 2);
n_rows = size(normData_Ndnf,2);
corrmatrix = NaN(n_rows, n_cols);
subplot(21,6,[19 38]+6*6)
    for iter = 1:n_rows
        for jter = 1:n_cols
            tmp_coef = corrcoef(normData_Ndnf(:,iter), score(:,jter));
            tmp_coef = tmp_coef(1,2);
            corrmatrix(iter,jter) = tmp_coef; 
            coef_magn = abs(tmp_coef);
            if tmp_coef < 0; colCode = [[255, 104, 23]./255 coef_magn];
            else; colCode = [0 0.4470 0.7410 coef_magn]; end
            vi = viscircles([jter iter],coef_magn/2, "Color", colCode);
            vi.Children(2).LineWidth = 1;
            vi.Children(1).LineWidth = 1;
        end
    end
    axis equal
    xlim([0 n_cols+1])
    ylim([0 n_rows+1])
    set(gca, "XTick", 1:n_cols,...
             "XTickLabel", 1:14,...
             "YTick", 1:n_rows,...
             "YTickLabel", 1:14,...
             "YDir", "reverse")
    xtickangle(90)
    plotAesthetics(gca,1,font_size_large)
    xlabel("Principle component")
    ylabel("Intrinsic property")
    text(.5, -1, "Positive correlation", "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontName", "Arial", "FontSize", font_size_large, "Color", [0 0.4470 0.7410],...
        "FontWeight", "bold")
    text(.5, -.2, "Negative correlation", "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontName", "Arial", "FontSize", font_size_large, "Color", [255, 104, 23]./255,...
        "FontWeight", "bold")
    
    
%% ------------------------ Correlation circle ----------------------------
subplot(21, 6, [21 40]+6*6)
vi = viscircles([0 0], 1, "Color", grey_light); hold on
vi.Children(2).LineWidth = 1;
vi.Children(1).LineWidth = 1;
    for iter = 1:14
        xpos = corrmatrix(iter,1);
        ypos = corrmatrix(iter,2);
        if ypos > 0; vert = "bottom";
        else; vert = "top"; end
        if xpos > 0; horz = "left";
        else horz = "right"; end
        text(xpos, ypos, num2str(iter),...
            "VerticalAlignment", vert,...
            "HorizontalAlignment", horz,...
            "FontSize", font_size_large, "FontName", "Arial",...
            "FontWeight", "bold"); 
        line([0 xpos], [0, ypos], "Color", "black")
    end
    line([-1.1 1.1], [0 0], "Color", grey_dark, "LineStyle", ':', "LineWidth", 1)
	line([0 0], [-1.1 1.1], "Color", grey_dark, "LineStyle", ':', "LineWidth", 1)
    axis square;
	xlim([-1.1 1.1]); ylim([-1.1 1.1])
    xlabel("Correlation with PC1");
    ylabel("Correlation with PC2")
    plotAesthetics(gca,1,font_size_large)
	circ = gca;
    pos_circ = circ.Position;
    set(circ, "Position", pos_circ + [-.02 0 0 0])
    
%% ------------------------ ward on features ------------------------------
% linkage evaluation
Z = linkage(corrmatrix,'ward','euclidean');
ClustIDs = cluster(Z,'maxclust',2);
% calculate optimal order
D = pdist(corrmatrix);
leafOrder = optimalleaforder(Z,D);
% plot
subplot(21, 6, [23 30]+6*6);
    d = dendrogram(Z, 0, 'Reorder',leafOrder);
    ylabel("Euclid. dist. [AU]")
    xlabel("Intrinsic property")
    plotAesthetics(gca,1,font_size_large)
    set(d, "Color", "black")
    
text(-2, -.5, "Intrinsic properties:",...
	"VerticalAlignment", "middle", "HorizontalAlignment", "left",...
	"FontSize", font_size_large, "FontName", "Arial",...
	"FontWeight", "bold")
for iter = 1:7
    text(-2, -.4*iter-.5, [num2str(iter) ' - ', varNames{iter}],...
        "VerticalAlignment", "middle", "HorizontalAlignment", "left",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(6, -.4*iter-.5, [num2str(iter+7) ' - ', varNames{iter+7}],...
        "VerticalAlignment", "middle", "HorizontalAlignment", "left",...
        "FontSize", font_size_small, "FontName", "Arial")
end


%% -------------------------- save figure ---------------------------------
set(gcf,'renderer','Painters')
%print('FS3_cluster_eval_extended.pdf','-dpdf')