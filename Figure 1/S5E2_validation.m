%% prepare workspace
clear all; close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat');
load('./analysis/Figure 1/resultTable_INs_FR_genetics.mat')
col_S5E2 = [0.1647 0.6353 0.7804];
col_FS = [0.1 0.1 1];
col_nonFS = [0.4039 0.7608 0.2275];
F_S5E2 = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% histology
subplot(5,2, [1 3])
    S5E2 = imread('./analysis/Figure 1/2022-09-27_03_03_20x_AAVs_48umThickness-01_z03c1-2.tif');
    S5E2_rot = imrotate(S5E2, -65);
    S5E2_plot = S5E2_rot(1000:4100, 2500:4500, :);
    imagesc(S5E2_plot);
    axis equal;
    line([1900 1900-320.5128]+500,[250 250],"LineWidth",2,"Color", "black")
	text(mean([1900 1900-320.5128])+500, 240, "100µm",...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
         'FontSize', font_size_large, 'FontName', 'Arial',...
         'FontWeight', 'bold', 'Color', 'black')
   	axis off

    
%% IO-plot data and groups
path_NDNF_PV = './analysis/Figure 1/currentSteps_NDNF>PV/';
path_NDNF_nonFS = './analysis/Figure 1/currentSteps_NDNF>nonFS/';
path_NDNFNPY_PV = './analysis/Figure 1/currentSteps_NDNFNPY>PV/';
path_NDNFNPY_nonFS = './analysis/Figure 1/currentSteps_NDNFNPY>nonFS/';
%load("result")
% load data
[IO_NDNF_PV, fNames_NDNF_PV] = get_IO_traces(path_NDNF_PV);
    NDNF_PV_groups = [2 2 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1];
[IO_NDNF_nonFS, fNames_NDNF_nonFS] = get_IO_traces(path_NDNF_nonFS);
    NDNF_nF_groups = [3 3 3 3 3 3 3 3 3 3 3 3 3 3];
[IO_NDNFNPY_PV, fNames_NDNFNPY_PV] = get_IO_traces(path_NDNFNPY_PV);
    NDNFNPY_PV_groups = [1 2 2 1 1 1 1 1 1 1 1 1 1 1 1 2 2 1 1 2 1 2 1];
[IO_NDNFNPY_nonFS, fNames_NDNFNPY_nonFS] = get_IO_traces(path_NDNFNPY_nonFS);
    NDNFNPY_nF_groups = [3 3 3 3 3 3 3 3 3 3 3 3];
% format data
dats_NDNF_nonFS = NaN(length(IO_NDNF_nonFS),180);
dats_NDNFNPY_nonFS = NaN(length(IO_NDNFNPY_nonFS),180);
dats_NDNF_PV = NaN(length(IO_NDNF_PV),180);
dats_NDNF_FS = NaN(length(IO_NDNF_PV),180);
dats_NDNFNPY_PV = NaN(length(IO_NDNFNPY_PV),180);
dats_NDNFNPY_FS = NaN(length(IO_NDNFNPY_PV),180);
for iter = 1:length(IO_NDNF_nonFS)
    dats_NDNF_nonFS(iter, 1:length(IO_NDNF_nonFS{iter})) = IO_NDNF_nonFS{iter}; end
for iter = 1:length(IO_NDNFNPY_nonFS)
    dats_NDNFNPY_nonFS(iter, 1:length(IO_NDNFNPY_nonFS{iter})) = IO_NDNFNPY_nonFS{iter}; end
for iter = 1:length(IO_NDNF_PV) %#ok<ALIGN>
    if NDNF_PV_groups(iter) == 1; dats_NDNF_PV(iter, 1:length(IO_NDNF_PV{iter})) = IO_NDNF_PV{iter}; end
    if NDNF_PV_groups(iter) == 2; dats_NDNF_FS(iter, 1:length(IO_NDNF_PV{iter})) = IO_NDNF_PV{iter}; end; end
for iter = 1:length(IO_NDNFNPY_PV) %#ok<ALIGN>
    if NDNFNPY_PV_groups(iter) == 1; dats_NDNFNPY_PV(iter, 1:length(IO_NDNFNPY_PV{iter})) = IO_NDNFNPY_PV{iter}; end
    if NDNFNPY_PV_groups(iter) == 2; dats_NDNFNPY_FS(iter, 1:length(IO_NDNFNPY_PV{iter})) = IO_NDNFNPY_PV{iter}; end; end
    
    
%% plot average IO cure, PV and FS separate
subplot(5,2,4)
    [~, n_FS, ~, ~,] = plot_IO_meanSEM([dats_NDNF_FS; dats_NDNFNPY_FS], col_FS, .2);
    [~, n_PV, ~, ~,] = plot_IO_meanSEM([dats_NDNF_PV; dats_NDNFNPY_PV], col_S5E2, .2);
    [~, n_nonFS, ~, ~,] = plot_IO_meanSEM([dats_NDNF_nonFS; dats_NDNFNPY_nonFS], col_nonFS, .2);
    % legend
    text(20, 300, ['S5E2+ INs, n=' num2str(sum(resultTable.finalGroups == 1))],...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', col_S5E2)
	text(20, 270, ['S5E2- INs: FR > 150Hz, n=' num2str(sum(resultTable.finalGroups == 2))],...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', col_FS)
	text(20, 240, ['S5E2- INs: FR < 150Hz, n=' num2str(sum(resultTable.finalGroups == 3))],...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', col_nonFS)
    % aesthetics
    plotAesthetics(gca, 1, font_size_large); xlim([0 810]);
    xlabel("Current injected [pA]"); ylabel("Firing rate [Hz]")
    title("Combined data from INs from F1 and F5")
    line(xlim, [150 150], "LineStyle", ":", "Color", "red", "LineWidth", 1)
    hold off
       
    
%% example recording 
% NDNF > PV, identified by S5E2), FR max = 299Hz, blue = 1st, red= max, black = first AP, scale: 20mV/100ms
exampleRecording = load("./analysis/Figure 1/currentSteps_NDNF>PV/2022_09_03_NdnfCre_line2882_mouse211_0004_barrageDataSingleCell.mat");
ts = (1:size(exampleRecording.saveData.recordingData.rawData_mV, 1))*exampleRecording.saveData.recordingData.si_us/1000;
subplot(5,2,2)
    plot(ts, exampleRecording.saveData.recordingData.rawData_mV(:,1), "Color", "blue"); hold on
    plot(ts, exampleRecording.saveData.recordingData.rawData_mV(:,75), "Color", "red");
    plot(ts, exampleRecording.saveData.recordingData.rawData_mV(:,31), "Color", "black");
    % scale
    line([2300 2300],[30 50],"Color", "black", "LineWidth", 1)
    line([2300 2400],[30 30],"Color", "black", "LineWidth", 1)
	text(2350, 27, "100ms",...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'black')
	text(2290, 40, "20mV ",...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'black')
    % legend
    text(900, -90, "First sweep -100pA",...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'blue')
	text(900, -105, "Supra threshold sweep +210pA",...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'black')
	text(900, -120, "FR_{max} sweep (299Hz) +650pA",...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'red')
    % aesthetics
	text(1600, 45, "S5E2+ IN example recording",...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'FontSize', font_size_large, 'FontName', 'Arial',...
        'FontWeight', 'bold', 'Color', 'black')
    xlim([800 2400]); axis off; hold off

    
    
%% scatter plot
subplot(5,6,13:14)
    scatter(resultTable.AP_halfWidth_ms(resultTable.finalGroups == 1),...
            resultTable.maxFR_Hz_evoked(resultTable.finalGroups == 1),...
                '.', "CData", col_S5E2); hold on
	scatter(resultTable.AP_halfWidth_ms(resultTable.finalGroups == 2),...
            resultTable.maxFR_Hz_evoked(resultTable.finalGroups == 2),...
            '.', "CData", [0 0 1]);
	scatter(resultTable.AP_halfWidth_ms(resultTable.finalGroups == 3),...
            resultTable.maxFR_Hz_evoked(resultTable.finalGroups == 3),...
            '.', "CData", col_nonFS);
	xlabel("AP width at half-maximum amplitude [ms]"); ylabel("Maximum evoked FR [Hz]")
    title("Electrophysiology")
    line([0 1.2], [150 150], "LineStyle", ":", "Color", [1 0 0], "LineWidth", 1)
    plotAesthetics(gca, 1, font_size_small); axis square;
    
    
%% box plots
grSep = resultTable.finalGroups;
grCom = resultTable.finalGroups;
grCom(grCom == 2) = 1;
grCom(grCom == 3) = 2;

subplot(5,6,15)
bp1 = boxplot(resultTable.maxFR_Hz_evoked, grSep);
    set(bp1, "Color", "black");
    h = findobj(gca,'Tag','Box');
    patch(get(h(1), 'XData'), get(h(1),'YData'), col_nonFS, 'FaceAlpha', .8);
    patch(get(h(2), 'XData'), get(h(2),'YData'), col_FS, 'FaceAlpha', .8);
    patch(get(h(3), 'XData'), get(h(3),'YData'), col_S5E2, 'FaceAlpha', .8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel("FR_{max} [Hz]"); set(gca, "XTickLabel", {"S5E2", "FS", "non-FS"})
    plotAesthetics(gca, 1, font_size_small);
    
subplot(5,6,16)
    bp2 = boxplot(resultTable.maxFR_Hz_evoked, grCom);
    set(bp2, "Color", "black");
    h = findobj(gca,'Tag','Box');
    patch(get(h(1), 'XData'), get(h(1),'YData'), col_nonFS, 'FaceAlpha', .8);
    patch(get(h(2), 'XData'), get(h(2),'YData'), col_S5E2, 'FaceAlpha', .8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
	ylabel("FR_{max} [Hz]"); set(gca, "XTickLabel", {"S5E2+ & FS", "non-FS"})
    plotAesthetics(gca, 1, font_size_small);

mean_Dats = []; semDats = [];
subplot(5,6,17)
    for iter = 1:3 %#ok<ALIGN>
        mean_Dats(iter) = mean(resultTable.maxFR_Hz_evoked(grSep==iter));
        sem_Dats(iter) = std(resultTable.maxFR_Hz_evoked(grSep==iter)) / sqrt(sum(grSep==iter)); end
    [mean_Dats' sem_Dats']
    plot(1, mean_Dats(1), "o", "MarkerEdgeColor", col_S5E2, "MarkerFaceColor", col_S5E2, "MarkerSize", 2); hold on
    plot([1 1], [mean_Dats(1)+sem_Dats(1) mean_Dats(1)-sem_Dats(1)], "LineWidth", .5, "Color", col_S5E2)
    plot(2, mean_Dats(2), "o", "MarkerEdgeColor", col_FS, "MarkerFaceColor", col_FS, "MarkerSize", 2); hold on
    plot([2 2], [mean_Dats(2)+sem_Dats(2) mean_Dats(2)-sem_Dats(2)], "LineWidth", .5, "Color", col_FS)
    plot(3, mean_Dats(3), "o", "MarkerEdgeColor", col_nonFS, "MarkerFaceColor", col_nonFS, "MarkerSize", 2); hold on
    plot([3 3], [mean_Dats(3)+sem_Dats(3) mean_Dats(3)-sem_Dats(3)], "LineWidth", .5, "Color", col_nonFS)
    ylabel("FR_{max} [Hz]"); set(gca, "XTick", [1 2 3], "XTickLabel", {"S5E2", "FS", "non-FS"});
    plotAesthetics(gca, 1, font_size_small); ylim([40 300]); xlim([0 4])
    
mean_Dats = []; sem_Dats = [];
subplot(5,6,18)
	for iter = 1:2 %#ok<ALIGN>
        mean_Dats(iter) = mean(resultTable.maxFR_Hz_evoked(grCom==iter));
        sem_Dats(iter) = std(resultTable.maxFR_Hz_evoked(grCom==iter)) / sqrt(sum(grCom==iter)); end
    [mean_Dats' sem_Dats']
    plot(1, mean_Dats(1), "o", "MarkerEdgeColor", col_S5E2, "MarkerFaceColor", col_S5E2, "MarkerSize", 2); hold on
    plot([1 1], [mean_Dats(1)+sem_Dats(1) mean_Dats(1)-sem_Dats(1)], "LineWidth", .5, "Color", col_S5E2)
    plot(2, mean_Dats(2), "o", "MarkerEdgeColor", col_nonFS, "MarkerFaceColor", col_nonFS, "MarkerSize", 2); hold on
    plot([2 2], [mean_Dats(2)+sem_Dats(2) mean_Dats(2)-sem_Dats(2)], "LineWidth", .5, "Color", col_nonFS)
    ylabel("FR_{max} [Hz]"); set(gca, "XTick", [1 2], "XTickLabel", {"S5E2+ & FS", "non-FS"});
    plotAesthetics(gca, 1, font_size_small); ylim([40 300]); xlim([0 3])
    
    
%% compare PV+, S5E2+ and, S5E2->150Hz and S5E2-<150Hz
colCodes_con_ndnf = [[0         0         1     ];
                     [0.1647    0.6353    0.7804];
                     [0.4039    0.7608    0.2275];
                     [0.9882    0.6157    0.0118]];
% gather data
table_PV = load('./analysis/currentSteps/resultTable_PV.mat');
table_S5E2 = load('./analysis/currentSteps/resultTable_S5E2.mat');
table_INs = load('./analysis/currentSteps/resultTable_INs.mat');
table_VIP = load('./analysis/currentSteps/resultTable_VIP.mat');
resultTable_all = [table_PV.resultTable_currentSteps;...
                   table_S5E2.resultTable_currentSteps;...
                   table_INs.resultTable_currentSteps;...
                   table_VIP.resultTable_currentSteps];
resultTable = resultTable_all(resultTable_all.depolyBlock == 1, :);

% plot results
subplot(5,4,14)
    bp = boxplot(resultTable.maxFR_Hz_evoked, resultTable.genetics); hold on
    set(bp, 'Color', 'black', "LineWidth", 1);
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),colCodes_con_ndnf(5-jter,:),'FaceAlpha',.8);
    end
    % aesthetics
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel('Max. FR [Hz]'); 
    plotAesthetics(gca, 1, font_size_large); 
    % add n observations
    yLims = ylim;
    genetics = {"PV-cre", "S5E2", "none/mDlx", "VIP-cre"};
    for iter = 1:4
        text(iter, yLims(1)+0.1*-diff(yLims),....
            ['n=' num2str(sum(strcmp(resultTable.genetics, genetics{iter})))],...
            "HorizontalAlignment", "center", "VerticalAlignment", "top",...
            "FontName", "Arial", 'FontSize', font_size_large,...
            'FontWeight', 'bold');
    end
    % statistics
    [p, tbl, stats] = kruskalwallis(resultTable.maxFR_Hz_evoked, resultTable.genetics, 'off');
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off"); % p_pv_sst = 0.0305; p_vip_sst = 0.0043
    counter = 1;
    for iter = 1:size(mc,1)
        if mc(iter,end) < 0.05
            yPos = yLims(2) + diff(yLims) * 0.075 * counter;
            line([mc(iter,1) mc(iter,2)], [yPos yPos], 'LineWidth', 1, 'Color', 'black');
            if mc(iter, end) < 0.001; stars = "***";
            elseif mc(iter, end) < 0.01; stars = "**";
            else; stars = "*"; end
            text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
            counter  = counter + 1;
        end
    end
    yPos = yLims(2) + diff(yLims) * 0.075 * counter;
    ylim([yLims(1) yPos])
    title("Firing rate by IN genetics")
    subtitle("Only cells which reach depolarization block are included.")    
    
    
%% save figure
set(F_S5E2,'renderer','Painters')
%print('FS1_S5E2.svg','-dsvg')


%% subroutines
function [IO_traces,fNames] = get_IO_traces(fPath)
    fNames = getFilesFromDir(fPath, '_barrageDataSingleCell.mat', 'names');
    IO_traces = {};
    for iter = 1:length(fNames)
        load(strcat(fPath, fNames(iter)));
        IO_traces{iter} = saveData.analysis.FR_Hz_step;
        clear saveData
    end
end

function [dats_mean, notIsNaN, dats_std, dats_sem] = plot_IO_meanSEM(dats, col, a)
	dats_mean = mean(dats,1,"omitnan");
    notIsNaN = sum(~isnan(dats),1);
    dats_std = std(dats,[],1,"omitnan");
    dats_sem = dats_std ./ sqrt(notIsNaN);
    I = 10:10:800;
    I2 = [I, fliplr(I)];
    sem_area = [[dats_mean(1:80) + dats_sem(1:80)] fliplr([dats_mean(1:80) - dats_sem(1:80)])];
    fill(I2, sem_area, col, "FaceAlpha", a, "EdgeAlpha", .1); hold on
    plot(I, dats_mean(1:80), "-", "LineWidth", 1, "Color", col);
end

