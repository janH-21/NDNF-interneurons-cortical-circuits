%% prepare workspace and load data
clear all; close all; clc;
cd("/main/");
load("./analysis/misc/aesthetics.mat");
% 250ms 2Hz 28Hz example
[d_250ms, si, ~] = abfload('./data/persistentFiring_NDNF_250ms28Hz/2022_10_27_NdnfCre_line3016_mouse165_0043.abf');
% figure
FS11_example = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ---------------------- pF example 250ms --------------------------------
d_plot = squeeze(reshape(d_250ms(:,1,1:10),[],10));
subplot(7,6,[1 12]);
    for iter = 1:10
        plot(d_plot(:,iter)-iter*100, "Color", "black"); hold on
    end
    ylim([-1000, 0])
    text(-10000, mean(ylim), "Sweep #",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "Rotation", 90)
    xLims = xlim; xlim([-5000 xLims(2)]);
    % scale
    line([482000 482000],[-100 0],"Color", "black", "LineWidth", width_scalebar);
	line([482000 492000],[-100 -100],"Color", "black", "LineWidth", width_scalebar);
	text(482000, -50, "100mV ", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontWeight", "bold");
    text(487000, -110, "200ms", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold");
    axis off


%% save figure
set(gcf,'renderer','Painters')
% print('FS11_example.pdf','-dpdf')