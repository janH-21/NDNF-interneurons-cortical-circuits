%% ---------------- prepare workspace and load data -----------------------
clear all; close all; clc;
cd("/main/");
load('analysis/misc/aesthetics.mat')
magenta = [217, 17, 197]./255;
yellow = [219, 182, 18]./255;
F3_RNAscopeValidation = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");
binWidth = 5;


%% ------------------------ prepare data ----------------------------------
intersectData = readtable("analysis/Figure 3/allData_intersect.csv");
ndnf_npy_mRNA = table2array(intersectData(...
        intersectData.counter_2==1 & ...
        intersectData.counter_3==1,9));
ndnf_only_mRNA = table2array(intersectData(...
        intersectData.counter_2==0 &...
        intersectData.counter_3==1,9));
all = table2array(intersectData(...
        intersectData.counter_1==1 & ...
        intersectData.counter_2==1 & ...
        intersectData.counter_3==1,9));   
EYFP = table2array(intersectData(...
        intersectData.counter_1==1,9));
ndnf_all = table2array(intersectData(...
        intersectData.counter_3==1,9));
npy_all = table2array(intersectData(...
        intersectData.counter_2==1,9));


%% ------------------------ intersect concept -----------------------------
subplot(6, 5, 1)
    rectangle('Position',[0 0 1.1 1.1],'Curvature',[1 1], "FaceColor", [magenta 0.2]); hold on
    rectangle('Position',[0.7 0 1.1 1.1],'Curvature',[1 1], "FaceColor", [yellow 0.2]); hold on
    text(0.25, 0.5, ["NPY" "-cre"], "FontName", "Arial", "FontSize", font_size_large_minus,...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle")
    text(1.65, 0.58, ["Ndnf" "-ires" "-FlpO"], "FontName", "Arial", "FontSize", font_size_large_minus,...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle")
    text(0.9, 1.45, "cOn/fOn-ChR2-EYFP", "FontName", "Arial", "FontSize", font_size_large_minus,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    text(0.9, 1.3, "cOn/fOn-EYFP", "FontName", "Arial", "FontSize", font_size_large_minus,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    %line([0.9 0.9],[1.2 0.5], "LineWidth", 1, "Color", "Black")
    quiver(.9, 1.2, 0, -.65, 0, "LineWidth", 1, "Color", "black", "MaxHeadSize", 1)
    axis equal; axis off
    
    
%% ------------------------ intersect example -----------------------------
intersectExample = imread('analysis/Figure 3/RNAscope_example_01/composition.png');
subplot(6, 5, [6 17])
    imagesc(imrotate(intersectExample, 0));
    axis equal; axis off; ylim([0 1040])
    % legend
    text(10, 1010, "EYFP",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "green",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle", "FontWeight", "bold")
    text(130, 1010, "DAPI",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", [66, 135, 245]./255,...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle", "FontWeight", "bold")
    text(270, 1010, "Ndnf",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", yellow,...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle", "FontAngle", "italic", "FontWeight", "bold")
    text(410, 1010, "Npy",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", magenta,...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle", "FontAngle", "italic", "FontWeight", "bold")
    text(65, 300, "50µm", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold") 
    % panel labels
    text(483, 25, "all channels",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", [.95 .95 .95],...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle")
    text(483, 360, "Ndnf ",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", yellow,...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontAngle", "italic")
    text(483, 690, "Npy ",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", magenta,...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontAngle", "italic")
    % cell 1
    text(60, 220, "1", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(60, 220+1*335, "1", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(60, 220+2*330, "1", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    % cell 2
    text(140, 260, "2", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(140, 260+1*335, "2", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(140, 260+2*330, "2", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    % cell 3
    text(420, 300, "3", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(420, 300+1*335, "3", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")
    text(420, 300+2*330, "3", "FontSize", font_size_large, "Color", [.95 .95 .95], "FontWeight", "bold")


%% --------------------------- marker: count ------------------------------
subplot(6,5,[3 4])
    % Npy
    h_npy = histcounts(npy_all, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_npy)', smooth(smooth(h_npy))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_npy), smooth(smooth(h_npy)), "LineWidth", 1, "Color", magenta); 
    
    % EYFP
	h_eyfp = histcounts(EYFP, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_eyfp)', smooth(smooth(h_eyfp))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_eyfp), smooth(smooth(h_eyfp)), "LineWidth", 1, "Color", "green");    
    
	% Ndnf
    h_ndnf = histcounts(ndnf_all, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area(1:binWidth:binWidth*length(h_ndnf)', smooth(smooth(h_ndnf))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf), smooth(smooth(h_ndnf)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    line([106.2 106.2], [0 42], "Color", "black", "LineStyle", ":", "LineWidth", 1)
    ylabel("Count"); xlabel("Dist. Pia [µm]")
    hold off
    ylim([0 24])
    
    
%% ------------------------- marker: density ------------------------------
subplot(6,5,[8 9])
    % Npy
    h_npy = histcounts(npy_all, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_npy)', smooth(smooth(h_npy))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_npy), smooth(smooth(h_npy)), "LineWidth", 1, "Color", magenta); 
    
    % EYFP
	h_eyfp = histcounts(EYFP, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_eyfp)', smooth(smooth(h_eyfp))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_eyfp), smooth(smooth(h_eyfp)), "LineWidth", 1, "Color", "green");    
    
	% Ndnf
    h_ndnf = histcounts(ndnf_all, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area(1:binWidth:binWidth*length(h_ndnf)', smooth(smooth(h_ndnf))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf), smooth(smooth(h_ndnf)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    line([106.2 106.2], [0 0.01], "Color", "black", "LineStyle", ":", "LineWidth", 1)
    ylabel("Density"); xlabel("Dist. Pia [µm]")
    hold off
    
    
%% --------------------- marker combination: counts -----------------------
% 1. full distribution
subplot(6, 5,[13 14])
    % EYFP
    h_EYFP = histcounts(EYFP, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_EYFP)', smooth(smooth(h_EYFP))',...
        'EdgeColor', 'none', 'FaceColor', "blue", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_EYFP), smooth(smooth(h_EYFP)), "LineWidth", 1, "Color", "blue"); hold on
    % Ndnf & Npy
    h_ndnf_npy_mRNA = histcounts(ndnf_npy_mRNA, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area(1:binWidth:binWidth*length(h_ndnf_npy_mRNA)', smooth(smooth(h_ndnf_npy_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_npy_mRNA), smooth(smooth(h_ndnf_npy_mRNA)), "LineWidth", 1, "Color", magenta); 
    % EYFP & Ndnf & Npy
	h_all = histcounts(all, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_all)', smooth(smooth(h_all))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_all), smooth(smooth(h_all)), "LineWidth", 1, "Color", "green");    
    % Ndnf & ~Npy
    h_ndnf_only_mRNA = histcounts(ndnf_only_mRNA, "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_ndnf_only_mRNA)', smooth(smooth(h_ndnf_only_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_only_mRNA), smooth(smooth(h_ndnf_only_mRNA)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    line([106.2 106.2], [0 30], "Color", "black", "LineStyle", ":", "LineWidth", 1)
    ylabel("Count"); xlabel("Dist. Pia [µm]")
    ylim([0 17])
    hold off

% 2. zoom on L1
subplot(6,5,15)
    % EYFP
    h_EYFP = histcounts(EYFP(EYFP<106.2), "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_EYFP)', smooth(smooth(h_EYFP))',...
        'EdgeColor', 'none', 'FaceColor', "blue", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_EYFP), smooth(smooth(h_EYFP)), "LineWidth", 1, "Color", "blue"); hold on
    % Ndnf & Npy
    h_ndnf_npy_mRNA = histcounts(ndnf_npy_mRNA(ndnf_npy_mRNA<106.2), "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 110]);
    area(1:binWidth:binWidth*length(h_ndnf_npy_mRNA)', smooth(smooth(h_ndnf_npy_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_npy_mRNA), smooth(smooth(h_ndnf_npy_mRNA)), "LineWidth", 1, "Color", magenta); 
    % EYFP & Ndnf & Npy
	h_all = histcounts(all(all<106.2), "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_all)', smooth(smooth(h_all))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_all), smooth(smooth(h_all)), "LineWidth", 1, "Color", "green");    
    % Ndnf & ~Npy
    h_ndnf_only_mRNA = histcounts(ndnf_only_mRNA(ndnf_only_mRNA<106.2), "BinWidth",binWidth, "Normalization", "count", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_ndnf_only_mRNA)', smooth(smooth(h_ndnf_only_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_only_mRNA), smooth(smooth(h_ndnf_only_mRNA)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    ylim([0 17])
    plotAesthetics(gca, 1, font_size_large);
    ylabel("Count"); xlabel("Dist. Pia [µm]")
    hold off
   

%% ------------------- marker combinations: density -----------------------
% 1. full distribution
subplot(6,5,[18 19])
    h_EYFP = histcounts(EYFP, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_EYFP)', smooth(smooth(h_EYFP))',...
        'EdgeColor', 'none', 'FaceColor', "blue", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_EYFP), smooth(smooth(h_EYFP)), "LineWidth", 1, "Color", "blue"); hold on
    % Ndnf & Npy
    h_ndnf_npy_mRNA = histcounts(ndnf_npy_mRNA, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area(1:binWidth:binWidth*length(h_ndnf_npy_mRNA)', smooth(smooth(h_ndnf_npy_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_npy_mRNA), smooth(smooth(h_ndnf_npy_mRNA)), "LineWidth", 1, "Color", magenta); 
    % EYFP & Ndnf & Npy
	h_all = histcounts(all, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_all)', smooth(smooth(h_all))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_all), smooth(smooth(h_all)), "LineWidth", 1, "Color", "green");    
    % Ndnf & ~Npy
    h_ndnf_only_mRNA = histcounts(ndnf_only_mRNA, "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 1200]);
    area( 1:binWidth:binWidth*length(h_ndnf_only_mRNA)', smooth(smooth(h_ndnf_only_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_only_mRNA), smooth(smooth(h_ndnf_only_mRNA)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    line([106.2 106.2], [0 0.01], "Color", "black", "LineStyle", ":", "LineWidth", 1)
    ylabel("Density"); xlabel("Dist. Pia [µm]")
    hold off

% 2. zoom on L1
subplot(6,5,20)
    % EYFP
    h_EYFP = histcounts(EYFP(EYFP<106.2), "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_EYFP)', smooth(smooth(h_EYFP))',...
        'EdgeColor', 'none', 'FaceColor', "blue", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_EYFP), smooth(smooth(h_EYFP)), "LineWidth", 1, "Color", "blue"); hold on
    % Ndnf & Npy
    h_ndnf_npy_mRNA = histcounts(ndnf_npy_mRNA(ndnf_npy_mRNA<106.2), "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 110]);
    area(1:binWidth:binWidth*length(h_ndnf_npy_mRNA)', smooth(smooth(h_ndnf_npy_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', magenta, 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_ndnf_npy_mRNA), smooth(smooth(h_ndnf_npy_mRNA)), "LineWidth", 1, "Color", magenta); 
    % EYFP & Ndnf & Npy
	h_all = histcounts(all(all<106.2), "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_all)', smooth(smooth(h_all))',...
        'EdgeColor', 'none', 'FaceColor', "green", 'FaceAlpha', 0.04); hold on
    plot(1:binWidth:binWidth*length(h_all), smooth(smooth(h_all)), "LineWidth", 1, "Color", "green");    
    % Ndnf & ~Npy
    h_ndnf_only_mRNA = histcounts(ndnf_only_mRNA(ndnf_only_mRNA<106.2), "BinWidth",binWidth, "Normalization", "pdf", "BinLimits", [0 110]);
    area( 1:binWidth:binWidth*length(h_ndnf_only_mRNA)', smooth(smooth(h_ndnf_only_mRNA))',...
        'EdgeColor', 'none', 'FaceColor', yellow, 'FaceAlpha', 0.04); hold on
               plot(1:binWidth:binWidth*length(h_ndnf_only_mRNA), smooth(smooth(h_ndnf_only_mRNA)), "LineWidth", 1, "Color", yellow); 
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    ylabel("Density"); xlabel("Dist. Pia [µm]")
    hold off
    
    
%% ------------------ distance pia mRNA scope data ------------------------
% box grouping
box_groups = [ones(size(ndnf_only_mRNA));...
              ones(size(ndnf_npy_mRNA))*2;...
              ones(size(all))*3;...
              ones(size(EYFP))*4];
box_data = [ndnf_only_mRNA;...
            ndnf_npy_mRNA;...
            all;...
            EYFP];
% box: distance to pia 
subplot(6,5,[21 23])
    box_depth_mRNA = boxplot(box_data, box_groups, "Orientation", "horizontal");
	set(box_depth_mRNA, "Color", "black");
	h = findobj(gca,'Tag','Box');
	patch(get(h(4),'XData'),get(h(4),'YData'), yellow ,'FaceAlpha',.8);
	patch(get(h(3),'XData'),get(h(3),'YData'), magenta,'FaceAlpha',.8);
	patch(get(h(2),'XData'),get(h(2),'YData'), 'green','FaceAlpha',.8);
	patch(get(h(1),'XData'),get(h(1),'YData'), 'blue','FaceAlpha',.8);
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    xlabel("Dist. to pia [µm]");
  	h=findobj('LineStyle','--'); set(h, 'LineStyle','-'); 
    set(gca, "YTickLabel", {"Ndnf only", "Npy only", "EYFP + Ndnf + Npy", "all EYFP"})
	% statistics
    yLims = ylim;
    [p, ~, stats] = kruskalwallis(box_data, box_groups, "off"); % p = 0.0149
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    line([1050 1050], [1 2], "Color", "black", "LineWidth", width_scalebar)
    text(1070, 1.35, "*",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontName", "Arial", "FontSize", font_stat_star)
    % n observations
    for iter = 1:4
        text(1100, iter,....
            ['n=' num2str(sum(box_groups==iter))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
            "FontName", "Arial", "FontSize", font_size_small);
    end


%% ------------------- intersect Venn diagram -----------------------------
% counter_1 = EYFP, counter_2 = npy, counter_3 = ndnf
venn_ndnf_only = sum(intersectData.counter_1==0 &...
                     intersectData.counter_2==0 &...
                     intersectData.counter_3==1);
venn_npy_only = sum(intersectData.counter_1==0 &...
                    intersectData.counter_2==1 &...
                    intersectData.counter_3==0);                
venn_EYFP_only = sum(intersectData.counter_1==1 &...
                     intersectData.counter_2==0 &...
                     intersectData.counter_3==0);
venn_ndnf_npy = sum(intersectData.counter_1==0 &...
                     intersectData.counter_2==1 &...
                     intersectData.counter_3==1);
venn_ndnf_EYFP = sum(intersectData.counter_1==1 &...
                     intersectData.counter_2==0 &...
                     intersectData.counter_3==1);
venn_npy_EYFP = sum(intersectData.counter_1==1 &...
                    intersectData.counter_2==1 &...
                    intersectData.counter_3==0);
venn_all = sum(intersectData.counter_1==1 &...
               intersectData.counter_2==1 &...
               intersectData.counter_3==1);

% plot venn diagramm
subplot(6, 5, 24)
    rectangle('Position',[0 1 1.2 1.2],'Curvature',[1 1], "FaceColor", [magenta 0.2]); hold on
    rectangle('Position',[0.66 1 1.2 1.2],'Curvature',[1 1], "FaceColor", [yellow 0.2]);
    rectangle('Position',[0.33 0.5 1.2 1.2],'Curvature',[1 1], "FaceColor", [0 1 0 0.2]);
    % numbers
    text(0.93, 1.5, num2str(venn_all), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    text(0.6, 1.66, num2str(venn_npy_only), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle")
    text(1.33, 1.66, num2str(venn_ndnf_only), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle")
    text(0.93, 0.76, num2str(venn_EYFP_only), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    text(0.5, 1.23, num2str(venn_npy_EYFP), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    text(0.93, 1.9, num2str(venn_ndnf_npy), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    text(1.33, 1.23, num2str(venn_ndnf_EYFP), "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle")
    % labels
    axis equal; axis off; hold off
    
    
%% ------------------------ statistics bar plot ---------------------------
TP = venn_all;
FP = venn_EYFP_only + venn_ndnf_EYFP + venn_npy_EYFP;
TN = venn_ndnf_only + venn_npy_only;
FN = venn_ndnf_npy;
sens = TP / (TP + FN);
spec = TN / (TN + FP);

subplot(6, 5, [25])
    bar([sens spec], "FaceAlpha", .6);
    plotAesthetics(gca,1,font_size_large);
    set(gca, "XTickLabel", {"Sensitivity    ", "    Specificity"})
	pos5 = get(gca, "Position");
    set(gca, "Position", pos5 + [0 0 0.05 -.015])
	text(1, .75, "0.74", "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "FontWeight", "bold")
    text(2, .98, "0.97", "FontName", "Arial", "FontSize", font_size_large,...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom", "FontWeight", "bold")
    

%% save figure
set(gcf,'renderer','Painters')
%print('F3_S1_intersect_rnaScope_v2_binWidth5.pdf','-dpdf')
