%% prepare workspace and load data
close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat');
load('./analysis/misc/IPSC_lowpassFilter.mat');
% figure
FS2_qualityControl = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ------- Panel A: Vm rest vs. AP threshold Quality Controls -------------
allNDNF_dats = collectData(resultTable_currentSteps.cellType(cleanDataInd) == 1,:); % same as in intrinsicProperties/MANOVA
subplot(15,15,[4 37])
    h_thresh = histogram(allNDNF_dats(:,3), "BinWidth", 5);
    set(gca, "XAxisLocation", "top", "YAxisLocation", "right");
    xlabel("AP threshold [mV]"); ylabel("Cell Count");
	xlim([h_thresh.BinEdges(1) h_thresh.BinEdges(end)]);
    plotAesthetics(gca, 1, font_size_large);
	a = gca;
    pos1 = a.Position;
subplot(15,15,[53 86])
    h_VmRest = histogram(allNDNF_dats(:,10), "BinWidth", 5, "Orientation", "horizontal");
    set(gca, "YDir", "reverse", "XDir", "reverse");
    xlabel("Cell Count"); ylabel("Vm_{rest} [mV]"); 
    ylim([h_VmRest.BinEdges(1) h_VmRest.BinEdges(end)]);
    xlim([0 55]);
    plotAesthetics(gca, 1, font_size_large);
	b = gca;
    pos2 = b.Position;
    dx = pos1(1)-(pos2(1)+pos2(3));
    dy = pos1(2)-(pos2(2)+pos2(4));
	set(gca, "Position", [pos2(1)+dx pos2(2)+dy pos2(3) pos2(4)]);
text(-63, -54, ["# Cells"],...
	"FontSize", font_size_large, "FontName", "Arial",...
	"HorizontalAlignment", "center", "VerticalAlignment", "top",...
	"FontWeight", "bold")
subplot(15, 15, [57 90]) 
    h2D_thresh_VmRest = histcounts2(allNDNF_dats(:,10), allNDNF_dats(:,3), "BinWidth", 5);
    hm = heatmap(h2D_thresh_VmRest,...
                 "FontName", "Arial", "FontSize", font_size_large);
    set(hm, "Position", [pos1(1)+0.0001 pos2(2)+.0154 pos1(3) pos2(4)]);     
    set(hm, "XDisplayLabels", repmat(' ', size(get(hm, "XDisplayLabels"))))
    set(hm, "YDisplayLabels", repmat(' ', size(get(hm, "YDisplayLabels"))))
    for iter = 1:3 % control heatmap font weight
        h_tmp=hm.NodeChildren(iter); 
        h_tmp.FontWeight='bold'; 
    end
    
    
%% ------------------ Panel B-D: Sag Quality Control ----------------------
[d,si,h]=abfload('./data/sag_VC/2022_11_24_NdnfCre_line2882_mouse246_0014.abf');

% Panel B: Protocol Command
subplot(21,6,[55 68])
    t = (1:size(d,1))*si/1000;
    for iter = 1:size(d,3)
        plot(t, d(:,2,iter), "Color", "black"); hold on
    end; hold off
    xlim([0 9500]); ylim([-110 -45]);
    % scale
    line([6000 6000],[-80 -70],"Color","Black","LineWidth",width_scalebar) % 10mV
    line([6000 7000],[-80 -80],"Color","Black","LineWidth",width_scalebar) % 1000ms
	text(6500, -81, "1000ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
	text(5900, -75, "10mV",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(-1500,-47, "-50mV", "FontName", "Arial", "FontSize", font_size_large)
	axis off;
    
% Panel C: Protocol Example Recording
subplot(21,6,[57 70])
    for iter = 1:size(d,3)
        plot(t, filtfilt(SOS,G,d(:,1,iter)), "Color", "black"); hold on
    end; hold off
    ylim([-450 100]); xlim([900 5500]);
    % scale
    line([5000 5000],[-200 -100],"Color","Black","LineWidth",width_scalebar) % 100pA
    line([5000 5500],[-200 -200],"Color","Black","LineWidth",width_scalebar) % 500ms
	text(5250, -203, "500ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
	text(4900, -150, "100pA",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    axis off;

load("./analysis/Figure 2/sagData_VC_CC.mat");
% Panel D: Sag CC vs. VC Correlation
subplot(21,6,[59 72])
    sagFrac_CC = table2array(sag_VC_CC_table(1:10,3));
    sagFrac_VC = table2array(sag_VC_CC_table(1:10,5));
    sc_sag = scatter(sagFrac_CC, sagFrac_VC, "LineWidth", 1, "SizeData", 8);
    ls_sag = lsline; axis square; set(ls_sag, "LineWidth", 1);
    [R_sag,P_sag] = corrcoef(sagFrac_CC(1:10)', sagFrac_VC(1:10)'); % P_sag = 4.1269e-5
    text(41,1, strjoin({'R =', num2str(R_sag(2,1))}), "FontName", "Arial", "FontSize", font_size_large,...
        "FontWeight", "bold", "VerticalAlignment", "bottom", "HorizontalAlignment", "right")
    text(20, 30, "n=10", "FontName", "Arial", "FontSize", font_size_large,...
        "FontWeight", "bold", "VerticalAlignment", "top", "HorizontalAlignment", "center")
    text(20,28,"***", "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
         "FontSize", font_stat_star, "FontName", "Arial");
    xlabel("CC sag fraction [%]"); ylabel("VC sag fraction [%]")
    plotAesthetics(gca, 1, font_size_large);
  
    
%% --------------------- Bursting LS Example ------------------------------
subplot(5,3,[10])
    [d,si,h]=abfload('./analysis/Figure 2/2021_11_04_NdnfFlpO_PAVcre_line2883_mouse239_0000.abf');
    t = (1:size(d,1))*si/1000;
    plot(t, d(:,1,23), "Color", "black")
    % scale
    line([500 500],[0 20],"Color","Black","LineWidth",width_scalebar) % 20mV
    line([500 1000],[0 0],"Color","Black","LineWidth",width_scalebar) % 500ms
        text(460, 10, "20mV ",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(750, -3, "500ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
	text(mean(xlim), -80, ["2 of 166 Ndnf INs (1.2%)" "1 of 39 Ndnf/Npy INs (2.56%)"],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    axis square; axis off;
subplot(5,3,11)
    plot(t, d(:,1,23), "Color", "black")
    xlim([1850 2050]); ylim([-70 40]);
    % scale
    line([1875 1875],[0 20],"Color","Black","LineWidth",width_scalebar) % 20mV
    line([1875 1900],[0 0],"Color","Black","LineWidth",width_scalebar) % 25ms
	text(1870, 10, "20mV",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(1887.5, -3, "25ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    axis square; axis off;
    
    
%% --------------------- Rebound spike example ----------------------------
subplot(5,3,[12])
    [d,si,h]=abfload('./analysis/Figure 2/2021_11_18_NdnfFlpO_PAVcre_line2883_mouse241_0005.abf');
    t = (1:size(d,1))*si/1000;
    plot(t, d(:,1,1), "Color", "black")
    % sale
    line([500 500],[0 20],"Color","Black","LineWidth",width_scalebar) % 20mV
    line([500 1000],[0 0],"Color","Black","LineWidth",width_scalebar) % 500ms
    text(460, 10, "20mV ",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontWeight", "bold")
    text(750, -3, "500ms",...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    axis square; axis off;
	text(mean(xlim), -105, ["3 of 166 Ndnf INs (1.81%)"],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
    

%% save figure
set(gcf,'renderer','Painters')
% print('FS2_qualityControl.pdf','-dpdf')

