%% prepare workspace and load data
close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat');
load('./analysis/currentSteps/resultTable_pF_groups.mat');

% load pF example cell
[d_1, si_1, ~]=abfload('./analysis/Figure 4//2020_07_10_NdnfVIP_JLM28301_0004.abf'); 
d_1 = squeeze(d_1(:,1,:));
% load pF example: pF before evoked firing
d_preSpkPF = abfload('./analysis/Figure 4/2023_01_04_NdnfCre_line2882_mouse293_0009.abf');
ts = (1:size(d_preSpkPF,1))*20/1000;

% cell types resultTable_currentSteps_all
% L23_PC          1
% L23_FS_IN       2
% L23_non-FS_IN   3
% L23_VIP         4

% figure
FS8_persistentFiringExtended = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% ---------------------- spike/spikelet overlay --------------------------
pf_example = reshape(d_1(:,14:35),[],2);
pf_example_pfSpike = pf_example(5.1065*10^5-50:5.1065*10^5+600 ,2); % 110
pf_example_spikelet = pf_example(5.1128*10^5.-50:5.1128*10^5+600 ,2); % 60
subplot(7,5,1)    
    plot(pf_example_spikelet(40:600),...
         "Color", colCodes_spklt(2,:), "LineWidth", 1); hold on
    plot(pf_example_pfSpike(90:650),...
         "Color", colCodes_spklt(3,:), "LineWidth", 1);
    xlim([0 561])
	% scale bar
	line([400 400],[-10 10],"Color","black","LineWidth",width_scalebar);
 	line([400 500],[-10 -10],"Color","black","LineWidth",width_scalebar);
    text(450, -11, "2ms", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    text(395, 0, "20mV ", "VerticalAlignment", "middle", "HorizontalAlignment", "right",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    axis off

    
%% -------------------- ectopic spikes sanity check -----------------------
maxBFspikes = [];
firstCount = [];
lastCount = [];
for iter = 1:length(allData)
    if resultTable_currentSteps.useBFanalysis(iter) == 1
        % get step indices
        firstPosStep = find(allData(iter).recordingData.Isteps_pA > 0, 1);
        lastStep = size(allData(iter).recordingData.rawData_mV,2);
        % extract spike count
        firstCount(iter) = max([allData(iter).analysis.APstats.nSpikesTotal{:}] -...
                                   [allData(iter).analysis.APstats.nSpikesDuringStep{:}]);
        firstCount(iter) = mean([allData(iter).analysis.APstats.nSpikesTotal{firstPosStep:firstPosStep+2}] -...
                                   [allData(iter).analysis.APstats.nSpikesDuringStep{firstPosStep:firstPosStep+2}]);
        lastCount(iter) = mean([allData(iter).analysis.APstats.nSpikesTotal{lastStep-2:lastStep}] -...
                                  [allData(iter).analysis.APstats.nSpikesDuringStep{lastStep-2:lastStep}]); 
        maxBFspikes(iter) = max(allData(iter).analysis.nBarrageSpikes);
    end
end
use_ind = (resultTable_currentSteps.cellType == 1 &...
           resultTable_currentSteps.useBFanalysis == 1 &....
           resultTable_currentSteps.isBF ==1);

% plot results
subplot(7,5,6)
    boxData = [firstCount; maxBFspikes; lastCount]';
    plot(1:3, boxData(use_ind,:), 'LineWidth', .5, 'Color', [0.1 0.1 0.1]); hold on
    bp = boxplot(boxData(use_ind,:));
    set(bp, 'LineWidth', 1, "Color", "black")
    set(gca, 'XTickLabels', {"first 3",'max.','last 3'})
    ylabel('Inter-step AP count')
    plotAesthetics(gca, 1, font_size_large)
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),blue_std,'FaceAlpha',.6);
    end
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    % statistics
    [p,tbl,stats] = friedman(boxData(use_ind,:), 1,'off');
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    yLims = ylim;
    counter = 1;
    for iter = 1:size(mc,1)
        if mc(iter,end) < 0.05
            line([mc(iter,1) mc(iter,2)],[yLims(2) yLims(2)]+0.05*(counter-1)*diff(yLims),...
                 "Color", "Black", "LineWidth", width_scalebar);
            if mc(iter,end) < 0.001; stars = '***';
            elseif mc(iter,end) < 0.01; stars = '**';
            elseif mc(iter,end) < 0.05; stars = '*';
            end
            text(mean(mc(iter,1:2)), yLims(2)+0.05*(counter-1)*diff(yLims), stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                "FontSize", 12, "FontName", "Arial", "FontWeight", "bold")
            counter = counter + 1;
        end
    end
    ylim([yLims(1) yLims(2)+0.05*(counter)*diff(yLims)]);
    text(mean(xlim), yLims(2)+0.05*(counter)*diff(yLims), ['n=', num2str(sum(use_ind))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")


%%  ------------------------ total ectopic spks ---------------------------
% count histogram
hist_data_ind = resultTable_currentSteps.useBFanalysis == 1 &...
                resultTable_currentSteps.cellType == 1;
subplot(7,5,[7 8]); % BinWidth = 10
    h1 = histogram(resultTable_currentSteps.nEctopicSupraTotal(hist_data_ind),...
        "BinWidth", 10, "BinLimits", [1 6500],...
        "Normalization", "count", "Orientation", "vertical");
    set(h1, "FaceColor", blue_std, "EdgeColor", blue_std,'FaceAlpha',.6)
	plotAesthetics(gca,1,font_size_large);
    ylabel("Cell count"); xlabel("Total ectopic AP count");
    ylim([0 max(h1.Values)+1]);
    h = gca; h.XAxis.TickLength = [0 0];
    % legend
    text(350, 51, ['n=' num2str(sum(hist_data_ind))],...
         "HorizontalAlignment", "left", "VerticalAlignment", "top",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
	text(6500, 7.8, "bin width = 10", "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
subplot(7,5,3)
    h2 = histogram(resultTable_currentSteps.nEctopicSupraTotal(hist_data_ind),...
        "BinWidth", 1,"Normalization", "count", "Orientation", "vertical");
    set(h2, "FaceColor", blue_std, "EdgeColor", blue_std,'FaceAlpha',.6)
   	plotAesthetics(gca,1,font_size_large);
    ylabel("Cell count"); xlabel("Total ectopic AP count");
    xlim([-1 21]); ylim([0 max(h2.Values)+1]);
    h = gca; h.XAxis.TickLength = [0 0];
    % legend
    line([5 5], ylim, "Color", "red", "LineStyle", ":", "LineWidth", 1)
    text(3.5, 16, "Cutoff", "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red", "FontWeight", "bold",...
        "Rotation", 90)
	text(20, 45, "bin width = 1", "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
subplot(7,5,2)
    no_ectopic_spikes = h2.Values(1);
    ectopic_spikes = sum(h2.Values) - no_ectopic_spikes;
    p_chart = pie([no_ectopic_spikes ectopic_spikes]);
    set(findobj(p_chart,'type','text'),'FontName','Arial', 'FontWeight', 'bold')
    plotAesthetics(gca,1,font_size_large)
    p_chart(1).FaceColor = blue_std; p_chart(3).FaceColor = blue_std;
    p_chart(1).FaceAlpha = .7; p_chart(3).FaceAlpha = 1;
    gca_pie = gca;
    pie_pos = gca_pie.Position;
    set(gca_pie, "Position", [pie_pos(1)+.02 pie_pos(2) pie_pos(3)-.05 pie_pos(4)-0.02]);
    text(-.25, .25, "0 ", "HorizontalAlignment", "right", "VerticalAlignment", "middle",...
        "FontSize", font_size_large+3, "FontName", "Arial", "Color", "white", "FontWeight", "bold")
	text(.1, -.25, ">0", "HorizontalAlignment", "left", "VerticalAlignment", "middle",...
        "FontSize", font_size_large+3, "FontName", "Arial", "Color", "white", "FontWeight", "bold")
	text(0, -1.3, ['n=' num2str(sum(hist_data_ind))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    title(["% of cells with" "# of ectopic APs"])
    
    
%% ------------ total ectopic spikes cross cells boxes --------------------
boxData = [resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.cellType == 1 & resultTable_currentSteps.useBFanalysis == 1);...
           resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.ClustIDs == 1 & resultTable_currentSteps.useBFanalysis == 1);...
           resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.ClustIDs == 2 & resultTable_currentSteps.useBFanalysis == 1);...
           resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.cellType == 2 & resultTable_currentSteps.useBFanalysis == 1);...
           resultTable.nEctTotal(resultTable.group == "PC");...
           resultTable.nEctTotal(resultTable.group == "PV");...
           resultTable.nEctTotal(resultTable.group == "non-FS")];
boxGroup = [ones(size(resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.cellType == 1 & resultTable_currentSteps.useBFanalysis == 1)));...
            ones(size(resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.ClustIDs == 1 & resultTable_currentSteps.useBFanalysis == 1)))*2;...
            ones(size(resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.ClustIDs == 2 & resultTable_currentSteps.useBFanalysis == 1)))*3;...
            ones(size(resultTable_currentSteps.nEctopicSupraTotal(resultTable_currentSteps.cellType == 2 & resultTable_currentSteps.useBFanalysis == 1)))*4;...
            ones(size(resultTable.nEctTotal(resultTable.group == "PC")))*5;...
            ones(size(resultTable.nEctTotal(resultTable.group == "PV")))*6;...
            ones(size(resultTable.nEctTotal(resultTable.group == "non-FS")))*7];

subplot(7,5,[4 10])
    bp = boxplot(boxData, boxGroup);
    % aesthetics
    set(bp, "Color", "black", "LineWidth", 1);
    h = findobj(gca,'Tag','Box');
    for jter=1:length(h)
        patch(get(h(jter),'XData'),get(h(jter),'YData'),blue_std,'FaceAlpha',.6);
    end
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    plotAesthetics(gca, 1, font_size_large);
    set(gca, "XTickLabel", ["Ndnf", "Clu1", "Clu2", "  NdnfNpy", "PC ", "PV ", "non-FS "]);
    ylabel("Total etopic AP count");
    xtickangle(35)
    % statistics
    [p, ~, stats] = kruskalwallis(boxData, boxGroup, "off") % p = 9.0836e-11
    c = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    % plot significances
    counter = 1;
    for iter = 1:size(c,1)
        if c(iter, end) < 0.05
            line([c(iter,1) c(iter,2)], [2900 2900]+counter*40,...
                 "Color", "black", "LineWidth", 1);
            counter = counter + 1;
        end
    end
    ylim([0 3500]);
    xlim([.6 7.4]);
    % add n observations
    for iter = 1:max(boxGroup)
        n = sum(boxGroup == iter);
        text(iter, 3600, ['n=' num2str(n)], "FontSize", font_size_large, "FontName", "Arial",...
             "VerticalAlignment", "top", "HorizontalAlignment", "center", "FontWeight", "bold")
    end
    gca_box = gca;
    gca_pos = gca_box.Position;
    set(gca_box, "Position", [gca_pos(1)+.02 gca_pos(2) gca_pos(3)-.02 gca_pos(4)]);
    
    
%% --------------------- marginal distribution ----------------------------
pf_ind = resultTable_currentSteps.useBFanalysis == 1 & resultTable_currentSteps.isBF ==1 & resultTable_currentSteps.cellType == 1;
interStepAPs_max = resultTable_currentSteps.nEctopicSupraMax(pf_ind);
interStep_count = resultTable_currentSteps.nBarrageSteps(pf_ind);
nBins = 20;
% plot
subplot(15,10,[53 64])
    h1 = histogram(interStepAPs_max, 'NumBins', nBins);
    xlim([0 max(interStepAPs_max)])
    xlabel("Max. # of inter-step APs")
    ylabel("# of Cells")
    set(gca,"FontSize", font_size_large, "FontName", "Arial", "box", "off");
    set(gca, "XAxisLocation", "top", "YAxisLocation", "right");
    plotAesthetics(gca,1,font_size_large)
    a = gca;
    pos1 = a.Position;
subplot(15, 10, [71 82])
    h2 = histogram(interStep_count, 'NumBins', nBins, 'Orientation', 'horizontal');
    set(gca, 'YDir','reverse', 'XDir','reverse')
    ylim([0 max(interStep_count)])
    xlabel("# of cells")
    ylabel("# of inter-steps with APs")
    set(gca,"FontSize", font_size_large, "FontName", "Arial", "box", "off");
    text(-50, 64, '# of cells',...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top", "Rotation", 90,...
        "FontWeight", "bold")
    text(30, 110, ['n=' num2str(sum(pf_ind))],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontWeight", "bold")
	plotAesthetics(gca,1,font_size_large)
    b = gca;
    pos2 = b.Position;
    dx = pos1(1)-(pos2(1)+pos2(3));
    dy = pos1(2)-(pos2(2)+pos2(4));
	set(gca, "Position", [pos2(1)+dx pos2(2)+dy pos2(3) pos2(4)]);
subplot(15,10,[73 84])
    h2D = histcounts2(interStep_count, interStepAPs_max, [nBins nBins]);
    hm = heatmap(h2D, "GridVisible", "off",...
                 "FontName", "Arial", "FontSize", font_size_large);
    set(hm, "Position", [pos1(1)+0.0001 pos2(2)+.0154 pos1(3) pos2(4)]);     
    set(hm, "XDisplayLabels", repmat(' ', size(get(hm, "XDisplayLabels"))))
    set(hm, "YDisplayLabels", repmat(' ', size(get(hm, "YDisplayLabels"))))
    for iter = 1:3 % control heatmap font weight
        h_tmp=hm.NodeChildren(iter); 
        h_tmp.FontWeight='bold'; 
    end

    
%%  --------------------- pF before spike onset ---------------------------
% # sporadic in analysis: 157, 130, 103, 57, 52, 41, 39, 23, 13, 9, 8, 6
% # 2 more cells many sporadic spikes

subplot(15,10,[56 67])
    line([-10 4010],[-70 -70],'LineStyle', ':'); hold on
    plot(ts, d_preSpkPF(:,1,5), "Color", "black")
    % legend
    text(1500, -81, "-60A","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold")
    text(-12, -70, "-70mV","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontWeight", "bold")
    % scale bar
    line([1000  1500],[-10 -10],"Color","black", "LineWidth", width_scalebar);
    line([1000 1000],[-10 0],"Color","black", "LineWidth", width_scalebar);
    text(990, -5, "10mV ","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontWeight", "bold")
    text(1250, -11, "500ms","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold")
    ylim([-80 10]); axis off
subplot(15,10, [58 69])
    line([-10 4010],[-70 -70],'LineStyle', ':'); hold on
    plot(ts, d_preSpkPF(:,1,14), "Color", "black")
    ylim([-80 10]); axis off
    text(1500, -81, "+30pA","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold")
     x = mean(ts([61250 61750]));
     pgon = polyshape([x x-150 x+150],[0 10 10]); 
     plot(pgon, "EdgeColor", "black", "FaceColor", "black", "FaceAlpha", .8, "EdgeAlpha", .8); 
subplot(15,10,[78 89])
    plot(ts(61250:61750), d_preSpkPF(61250:61750,1,14), "Color", "black", "LineWidth", 1)
    ylim([-80 10]); axis off
    line([1234  1235],[-10 -10],"Color","black", "LineWidth", width_scalebar);
    line([1234 1234],[-10 0],"Color","black", "LineWidth", width_scalebar);
    text(1233.8, -5, "10mV ","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontWeight", "bold")
    text(1234.5, -11, "1ms","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold")
	text(mean(xlim), -70, "14 of 168 Ndnf INs (8.3%)","FontName","Arial","FontSize",font_size_large,...
         "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold")
     

%% ---------------------- spikes vs. spikelets ----------------------------
% choose cutoff for distinguishing between spikes and spikelets
AP_maxVal_cutoff = -20; % based on scatter plots
% format data
AP_ampl_collect = {};
AP_ts_collect = {};
AP_max_collect = {};
AP_ampl_all = [];
AP_max_all = [];
AP_thresh_all = [];
% normalized within cell. CAVE: produces artifacts in spikelet-only cells
AP_ampl_all_norm = [];
AP_max_all_norm = [];
% number/proportion of spikes/spikelets
nSpikes = [];
nSpikelets = [];
nBoth = [];
percentSpikelets = [];

% read data
counter = 1;
for iter = 1:length(allData)
    % only use BF+ cells
    if resultTable_currentSteps.useBFanalysis(iter) == 1 &&...
       resultTable_currentSteps.isBF(iter) == 1 &&...
       resultTable_currentSteps.cellType(iter) == 1
        % find # of sweep with first ectopic spikes
        sweep_1 = find(allData(iter).analysis.nBarrageSpikes >0, 1);
        % local vars
        AP_ts = []; AP_ampl = []; AP_max = []; AP_thresh = [];
        % iterate first pF+ sweep to last sweep in recording
        for kter = sweep_1:length(allData(iter).analysis.nBarrageSpikes)+1
            % index for APs outside current steps
            ind = allData(iter).analysis.APstats.threshTs{kter} < 1061 |...
                  allData(iter).analysis.APstats.threshTs{kter} > 2063;
            % collect data
            AP_ts = [AP_ts allData(iter).analysis.APstats.threshTs{kter}(ind) + 4000*(kter-sweep_1)];
            AP_ampl = [AP_ampl allData(iter).analysis.APstats.APamplitude{kter}(ind)];
            AP_max = [AP_max allData(iter).analysis.APstats.APmaxVal{kter}(ind)];
            AP_thresh = [AP_thresh allData(iter).analysis.APstats.threshVal{kter}(ind)];
        end
        % assign to collector vars
        AP_ampl_collect{counter} = AP_ampl;
        AP_ts_collect{counter} = AP_ts;     
        AP_max_collect{counter} = AP_max;
        AP_ampl_all = [AP_ampl_all AP_ampl];
        AP_max_all = [AP_max_all AP_max];
        AP_thresh_all = [AP_thresh_all AP_thresh];
        AP_ampl_all_norm = [AP_ampl_all_norm (AP_ampl./max(AP_ampl))];
        AP_max_all_norm = [AP_max_all_norm (AP_max./max(AP_max))];
        percentSpikelets(counter) = sum(AP_max < AP_maxVal_cutoff)/length(AP_max);       
        nBoth = length(AP_max);
        nSpikelets(counter) = sum(AP_max < AP_maxVal_cutoff);
        nSpikes(counter) = sum(AP_max >= AP_maxVal_cutoff);
        % count
        counter = counter+1;
    end
end
% sort by %spikelets
[percentSpikelets,sortInd] = sort(percentSpikelets);
percentSpikelets = [percentSpikelets; 1-percentSpikelets]*100;

% spikes vs. spikelets bars
subplot(10,5,[26 33]+5) % AP max value
    histogram(AP_max_all, "BinWidth", 1)
    xlabel("AP max. value [mV]"); ylabel("Count")
    plotAesthetics(gca,1,font_size_large);
	line([-20 -20], ylim, "Color", "red", "LineStyle", ':', "LineWidth", 1)
    % legend
    text(-75, 1390, [num2str(size(AP_max_all,2)) ' persistent APs in ' num2str(size(AP_max_collect, 2)) ' cells'],...
        "VerticalAlignment", "top", "HorizontalAlignment", "left",...
        "FontName", "Arial", "FontSize", font_size_large,"FontWeight","bold")

    
%% save figure
set(gcf,'renderer','Painters')
%print('FS8_persistentFiringExtended.pdf','-dpdf')