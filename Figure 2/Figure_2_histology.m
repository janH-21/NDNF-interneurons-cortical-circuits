%% prepare workspace and load data
close all; clc;
cd("/main/");

figure
ndnf = imread('./analysis/Figure 2/2022-09-27_slide01_slice03_20x_L1-01-Orthogonal Projection-03_2_c1-2.tif');
ndnf = imrotate(ndnf, 60);
imagesc(ndnf);
    axis equal;
    ylim([1100 1850])
    xlim([1000 2700])
    axis off
    line([2500-321 2500],[1050 1050], "Color", "white", "LineWidth", 3)
    text(2500-321/2, 1040, "100µm",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "Color", "white", "FontName", "Arial", "FontWeight", "bold", "FontSize", font_size_large)
    