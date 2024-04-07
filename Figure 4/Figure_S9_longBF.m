%% prepare workspace and load data
close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat')

% load recording
fName = './analysis/Figure 4/2022_04_14_NdnfCre_line2882_mouse160_0032.abf';
[d,si,h]=abfload(fName);       
t = (1:size(d,1))*si/1000; % time stamps in ms

% figure
FS93_longBF = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ---------------------------- Panel A -----------------------------------
subplot(7,1,1)
    plot(t, d, "Color", "black"); hold on
    % scale bar
    line([max(t)-30000 max(t)],[-75 -75],"Color", "black", "LineWidth", 1)
    text(max(t)-15000, -85, "30s",...
        "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")

% arrows
x = 1000;
pgon = polyshape([x x-3000 x+3000],[5 20 20]); 
plot(pgon, "EdgeColor", "red", "FaceColor", "red",...
           "FaceAlpha", 1, "EdgeAlpha", 1);
text(4000, 20, " Panel B",...
    "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
    "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "FontWeight", "bold")
x = 481000;
pgon = polyshape([x x-3000 x+3000],[5 20 20]); 
plot(pgon, "EdgeColor", "red", "FaceColor", "red",...
           "FaceAlpha", 1, "EdgeAlpha", 1);
text(484000, 20, " Panel D",...
    "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
    "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "FontWeight", "bold")
       
% aesthetics
ylim([-78 20]); xlim([-2000 max(t)]); axis off

%% --------------------------- Panel B ------------------------------------
subplot(7,2,3)
    plot(t, d, "Color", "black"); hold on
    % aesthetics
    xlim([0 2000]); ylim([-78 20]); axis off
    % scale
    line([0 100],[-75 -75],"Color", "black", "LineWidth", 1)
    text(50, -85, "100ms",...
        "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    % arrows
    x = 835;
    pgon = polyshape([x x-25 x+25],[5 20 20]); 
    plot(pgon, "EdgeColor", "red", "FaceColor", "red",...
               "FaceAlpha", 1, "EdgeAlpha", 1);
    text(860, 20, " Panel C",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "FontWeight", "bold");

%% --------------------------- Panel D ------------------------------------
subplot(7,2,4)
    plot(t, d, "Color", "black"); hold on
    % aesthetics
    xlim([480000 482000]); ylim([-78 20]); axis off
    line([480000 480100],[-75 -75],"Color", "black", "LineWidth", 1)
    text(480050, -85, "100ms", "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    % arrows
    x = 480685;
    pgon = polyshape([x x-25 x+25],[5 20 20]); 
    plot(pgon, "EdgeColor", "red", "FaceColor", "red",...
               "FaceAlpha", 1, "EdgeAlpha", 1);
    text(480710, 20, " Panel E",...
        "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "FontWeight", "bold");

%% --------------------------- Panel C ------------------------------------

subplot(7,2,5)
    plot(t, d, "Color", "black")
    xlim([800 870]); ylim([-78 20])
    line([800 805],[-75 -75],"Color", "black", "LineWidth", 1)
    text(802.5, -85, "5ms", "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    axis off
	% scale bar
    line([869 869],[0 20],"Color", "black", "LineWidth", 1)
    text(869, 10, "20mV ",...
        "HorizontalAlignment", "right",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    
%% --------------------------- Panel E ------------------------------------
subplot(7,2,6)
    plot(t, d, "Color", "black")
    xlim([480650 480720]); ylim([-78 20])
    line([480650 480655],[-75 -75],"Color", "black", "LineWidth", 1)
    text(480652.5, -85, "5ms", "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    axis off


%% save figure
set(gcf,'renderer','Painters')
% print('FS9_longBF.pdf','-dpdf')