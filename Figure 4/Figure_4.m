%% prepare workspace and load data
close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat');
load('./analysis/currentSteps/resultTable_pF_groups.mat');
% load pF example cell
[d_1, si_1, ~]=abfload('./analysis/Figure 4/2020_07_10_NdnfVIP_JLM28301_0004.abf'); 
d_1 = squeeze(d_1(:,1,:));
% load xcorr example cell
[d_2, si_2, ~]=abfload('./analysis/Figure 4/2020_07_10_NdnfVIP_JLM28301_0008.abf');
d_2 = squeeze(d_2(:,1,:));
% figure
F4_persistentFiring = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% --------------------------- example cell -------------------------------
pf_example = reshape(d_1(:,14:35),[],2);
pf_example_zoom = pf_example(5.0034*10^5:5.4470*10^5,2);
pf_example_stepSpike = pf_example(5.0044*10^5:5.0044*10^5+400 ,2);
pf_example_pfSpike = pf_example(5.1065*10^5:5.1065*10^5+400 ,2);
pf_example_spikelet = pf_example(5.1128*10^5:5.1128*10^5+400 ,2);
yLims = [min(min(pf_example))-4 max(max(pf_example))+5];

% plot overview line 1
subplot(11,3,1:3);
    plot(pf_example(400001:end,1), "Color", "black"); hold on;
    ylim(yLims); xlim([0 size(pf_example,1)]-400000);
    axis off;
    % scale bar
    line([20000 20000]+300000,[-20 0],"Color","black","LineWidth",width_scalebar);
 	line([20000 70000]+300000,[-20 -20],"Color","black","LineWidth",width_scalebar);
    text(45000+300000, -21, "1s", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    text(19000+300000, -10, "20mV ", "VerticalAlignment", "middle", "HorizontalAlignment", "right",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    % step / inter-step legend
    % x coords 254075 303305 / 453239; y coords -78 / 53
    line([45238 103572],[-43 -43],"Color","black","LineWidth",width_scalebar);
	text(mean([53408 103875]), -42, "Step", "VerticalAlignment", "bottom", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "Color", "black");
    line([103875 253443],[-81 -81],"Color","black","LineWidth",width_scalebar);
    text(mean([103875 253443]), -82, "Inter-step intervall", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "Color", "black");
     
% plot overview line 2
subplot(11,3,4:6);
    plot(pf_example(:,2), "Color", "black"); hold on;
    ylim(yLims); xlim([0 size(pf_example,1)]);
    axis off;
    % Panel B line
    line([5.0034*10^5 5.4470*10^5],[-81 -81],"Color","black","LineWidth",width_scalebar);
	text(mean([5.0034*10^5 5.4470*10^5]), -82, "Panel B", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "Color", "black");
    % draw arrows
    x_arrows = [345002, 547041, 805577, 1042170, 1424230];
    for iter = 1:length(x_arrows)
        x = x_arrows(iter); 
        pgon = polyshape([x x-10000 x+10000],[12 21 21]); 
        plot(pgon, "EdgeColor", colCodes_bar(1,:), "FaceColor", colCodes_bar(1,:), "FaceAlpha", .8, "EdgeAlpha", .8); 
    end
    
% plot persistent firing
subplot(11,15,31:36);
    plot(pf_example_zoom(:), "Color", "black"); hold on;
    % draw arrows
    x_arrows = [250, 11090, 1.0460e+04];
    for iter = 1:length(x_arrows)
        x = x_arrows(iter); 
        pgon = polyshape([x x-450 x+450],[12 21 21]); 
        plot(pgon, "EdgeColor", colCodes_spklt(iter,:), "FaceColor", colCodes_spklt(iter,:), "FaceAlpha", 1, "EdgeAlpha", 1); 
    end
    ylim(yLims); xlim([-5000 length(pf_example_zoom)+5000]);
    axis off;
    % scale bar
    line([1500 1500],[-80 -60],"Color","black","LineWidth",width_scalebar);
 	line([1500 4000],[-80 -80],"Color","black","LineWidth",width_scalebar);
    text(2750, -81, "50ms  ", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    text(1450, -70, "20mV ", "VerticalAlignment", "middle", "HorizontalAlignment", "right",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
     
% plot evoked spike
subplot(11,15,46:47)
    plot(pf_example_stepSpike(:), "Color", colCodes_spklt(1,:), "LineWidth", width_traces_zoom)
    ylim(yLims); xlim([0 length(pf_example_stepSpike)]);
	axis off
	text(mean(xlim), yLims(1), "Evoked spike", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", colCodes_spklt(1,:));

% plot persistent spikelet
subplot(11,15,48:49)
    plot(pf_example_pfSpike(:), "Color", colCodes_spklt(3,:), "LineWidth", width_traces_zoom)
    ylim(yLims); xlim([0 length(pf_example_pfSpike)]);
    text(mean(xlim), yLims(1), "Persistent spikelet", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", colCodes_spklt(3,:));
	axis off
	% scale bar
	line([200 200],[-10 10],"Color","black","LineWidth",width_scalebar);
 	line([200 300],[-10 -10],"Color","black","LineWidth",width_scalebar);
    text(250, -11, "2ms", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
    text(195, 0, "20mV ", "VerticalAlignment", "middle", "HorizontalAlignment", "right",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");

% plot persistent spike
subplot(11,15,50:51)
    plot(pf_example_spikelet(:), "Color", colCodes_spklt(2,:), "LineWidth", width_traces_zoom)
    ylim(yLims); xlim([0 length(pf_example_spikelet)]);
	text(mean(xlim), yLims(1), "Persistent spike", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", colCodes_spklt(2,:));
	axis off

% plot phase plot  
subplot(11,15,[38 40+15])
    plot(pf_example_stepSpike(1:end-1), diff(pf_example_stepSpike)*50,...
         "Color", colCodes_spklt(1,:), "LineWidth", width_traces_zoom); hold on
    plot(pf_example_pfSpike(1:end-1), diff(pf_example_pfSpike)*50,...
         "Color", colCodes_spklt(3,:), "LineWidth", width_traces_zoom);
    plot(pf_example_spikelet(1:end-1), diff(pf_example_spikelet)*50,...
         "Color", colCodes_spklt(2,:), "LineWidth", width_traces_zoom);
    xlim([-80 20]); ylim([-100 200]);
    xlabel("V_m [mV]"); ylabel("\DeltaVm/\Deltat [mV/ms]");
    plotAesthetics(gca,1,font_size_large);


%% -------------------- proportion pF by cell type ------------------------
% format data:
pF_ocurr = NaN(7,2);
% dim 2: 1 = n_pF; 2 = n_non-pF
% dim 1: 1 = Ndnf; 2 = Ndnf Clu1; 3 = Ndnf Clu2; 4 = Ndnf/Npy; 5 = PC; 6 = FS; 7 = non-FS + VIP

% Ndnf counts
pf_occur(1,1) = sum(resultTable_currentSteps.cellType == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 1);
pf_occur(1,2) = sum(resultTable_currentSteps.cellType == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 0);
% Ndnf Clu 1 counts
pf_occur(2,1) = sum(resultTable_currentSteps.ClustIDs == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 1);
pf_occur(2,2) = sum(resultTable_currentSteps.ClustIDs == 1 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 0);
% Ndnf Clu 2 counts
pf_occur(3,1) = sum(resultTable_currentSteps.ClustIDs == 2 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 1);
pf_occur(3,2) = sum(resultTable_currentSteps.ClustIDs == 2 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 0);
% NdnfNpy counts
pf_occur(4,1) = sum(resultTable_currentSteps.cellType == 2 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 1);
pf_occur(4,2) = sum(resultTable_currentSteps.cellType == 2 &...
                    resultTable_currentSteps.useBFanalysis == 1 &...
                    resultTable_currentSteps.isBF == 0);
% PC counts
pf_occur(5,1) = sum(resultTable.group == "PC" & resultTable.nEctTotal > 5);
pf_occur(5,2) = sum(resultTable.group == "PC" & resultTable.nEctTotal <= 5);
% FS counts
pf_occur(6,1) = sum(resultTable.group == "PV" & resultTable.nEctTotal > 5);
pf_occur(6,2) = sum(resultTable.group == "PV" & resultTable.nEctTotal <= 5);
% non-FS + VIP counts
pf_occur(7,1) = sum(resultTable.group == "non-FS" & resultTable.nEctTotal > 5);
pf_occur(7,2) = sum(resultTable.group == "non-FS" & resultTable.nEctTotal <= 5);
                
% calculate proptorions
pf_occur_perc = (pf_occur(:,1)./(pf_occur(:,1) + pf_occur(:,2)))*100;
pf_occur_perc(:,2) = 100 - pf_occur_perc;       

% bar plot: cross-cell comparison occurrence
subplot(11,15,[42 45+15]);
    b = bar(pf_occur_perc,"BarLayout", "stacked");
    set(b(1), "FaceColor", colCodes_bar(1,:), "FaceAlpha", .6);
	set(b(2), "FaceColor", colCodes_bar(2,:), "FaceAlpha", .6);
    % aesthetics
	%set(b, "FaceAlpha", 0.65);
    ylabel("Percent of cells");
    set(gca, "XTickLabel", ["Ndnf", "Clu1", "Clu2", " Ndnf/Npy", "PC", "FS", "non-FS"],...
             "YTick", [0 50 100]);
    xtickangle(35);
    plotAesthetics(gca,1,font_size_large);
    % legend
    text(0.8, 140, 'pF+',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontWeight', 'bold',...
        'FontName', 'Arial', 'Color', colCodes_bar(1,:));
    text(0.8, 150, 'pF-',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, "FontWeight", "bold",...
        'FontName', 'Arial', 'Color', colCodes_bar(2,:));
    xlim([0.7 7.5])  

% format data for statistics + chi^2
chi_groups = [];
chi_data = [];
for iter = 1:size(pf_occur,1)
    tmp_group = ones(sum(pf_occur(iter,1:2)),1)*iter;
    chi_groups = [chi_groups; tmp_group];
    for jter = 1:2
        tmp_data = ones(sum(pf_occur(iter,jter)),1)*(jter-1);
        chi_data = [chi_data; tmp_data];
    end
end
[TABLE,CHI2,P] = crosstab(chi_groups, chi_data); % P = 1.9251e-12

% fisher exact post hoc + bonferroni
fisher_results = [];
counter = 1;
for iter = 1:size(pf_occur,1)
    for jter = 1:size(pf_occur,1)
        if iter < jter
            fisher_results(counter, 1) = iter;
            fisher_results(counter, 2) = jter;
            [h, p, stats] = fishertest(array2table(TABLE([iter jter],:)));
            fisher_results(counter, 3) = p;
            counter = counter + 1;
        end
    end
end
fisher_results(:,4) = fisher_results(:,3) * length(fisher_results(:,3));

% add staistics bars to plot
counter = 1;
for iter = 1:size(fisher_results,1)
    if fisher_results(iter,4) < 0.05
        line([fisher_results(iter,1) fisher_results(iter,2)],[100 100]+counter*4,...
            "Color","black","LineWidth",width_scalebar);
        counter = counter + 1;
    end
end

% add n observations
for iter = 1:size(pf_occur,1)
    n = sum(pf_occur(iter,1:2));
    text(iter, 152, ['n=' num2str(n)],...
        "FontSize", font_size_small, "FontName", "Arial",...
        "VerticalAlignment", "bottom", "HorizontalAlignment", "center")
end

    
%% ----------------------------- spikelets --------------------------------
% based on scatter plots
AP_maxVal_cutoff = -20; 
% format data
AP_ampl_collect = {};
AP_ts_collect = {};
AP_ts_relative = [];
AP_max_collect = {};
AP_ampl_all = [];
AP_max_all = [];
AP_thresh_all = [];
% per sweep
AP_ts_perSweep = {};
AP_max_perSweep = {};
AP_ampl_perSweep = {};
AP_thr_perSweep = {};
% per inter-step
AP_ts_inter = {};
AP_max_inter = {};
AP_ampl_inter = {};
AP_thr_inter = {};
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
       resultTable_currentSteps.isBF(iter) == 1
        % find # of sweep with first ectopic spikes
        sweep_1 = find(allData(iter).analysis.nBarrageSpikes >0, 1);
        % local vars
        AP_ts = []; AP_ampl = []; AP_max = []; AP_thresh = [];
        % iterate first pF+ sweep to last sweep in recording
        kter_max = length(allData(iter).analysis.nBarrageSpikes)+1;
        for kter = sweep_1:kter_max
            % index for APs outside current steps
            ind = allData(iter).analysis.APstats.threshTs{kter} < 1061 |...
                  allData(iter).analysis.APstats.threshTs{kter} > 2063;
            % collect data
            AP_ts = [AP_ts allData(iter).analysis.APstats.threshTs{kter}(ind) + 4000*(kter-sweep_1)];
            AP_ampl = [AP_ampl allData(iter).analysis.APstats.APamplitude{kter}(ind)];
            AP_max = [AP_max allData(iter).analysis.APstats.APmaxVal{kter}(ind)];
            AP_thresh = [AP_thresh allData(iter).analysis.APstats.threshVal{kter}(ind)];
            % collect sweep-wise data: {cell #}{sweep #}
            AP_ts_perSweep{counter}{kter} = allData(iter).analysis.APstats.threshTs{kter}(ind);
            AP_max_perSweep{counter}{kter} = allData(iter).analysis.APstats.APmaxVal{kter}(ind);
            AP_ampl_perSweep{counter}{kter} = allData(iter).analysis.APstats.APamplitude{kter}(ind);
            AP_thr_perSweep{counter}{kter} = allData(iter).analysis.APstats.threshVal{kter}(ind);
            % collect inter-step-wise data:
            if kter < kter_max
                % AP params
                ind_end = allData(iter).analysis.APstats.threshTs{kter} > 2063;
                ind_start = allData(iter).analysis.APstats.threshTs{kter+1} < 1061;
                AP_ts_inter{counter}{kter} = [allData(iter).analysis.APstats.threshTs{kter}(ind_end)...
                                              allData(iter).analysis.APstats.threshTs{kter+1}(ind_start)+4000];                  
                AP_max_inter{counter}{kter} = [allData(iter).analysis.APstats.APmaxVal{kter}(ind_end)...
                                               allData(iter).analysis.APstats.APmaxVal{kter+1}(ind_start)];                            
                AP_ampl_inter{counter}{kter} = [allData(iter).analysis.APstats.APamplitude{kter}(ind_end)...
                                                allData(iter).analysis.APstats.APamplitude{kter+1}(ind_start)];                         
                AP_thr_inter{counter}{kter} = [allData(iter).analysis.APstats.threshVal{kter}(ind_end)...
                                               allData(iter).analysis.APstats.threshVal{kter+1}(ind_start)];
            end
            % all
            AP_ts_relative = [AP_ts_relative allData(iter).analysis.APstats.threshTs{kter}(ind)];
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
        % count cells
        counter = counter+1;
    end
end
% sort by %spikelets
[percentSpikelets,sortInd] = sort(percentSpikelets);
percentSpikelets = [percentSpikelets; 1-percentSpikelets]*100;

% proportion: barplot spikes vs. spikelets
subplot(33,10,[131 132+20]); % percent
    b = bar(percentSpikelets',...
            "BarLayout", "stacked", "BarWidth", 1);
	set(b(1), "FaceColor", colCodes_spklt(3,:), "FaceAlpha", .9, "EdgeColor", "none");
	set(b(2), "FaceColor", colCodes_spklt(2,:), "FaceAlpha", .9, "EdgeColor", "none");
    plotAesthetics(gca, 1, font_size_large); 
    xlim([0 length(sortInd)]+.5); ylim([0 100])
    xlabel("Cell #"); ylabel("Persistent firing [%]"); 
    text(1, 101, "Spikes", "HorizontalAlignment", "left", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_spklt(2,:), "FontWeight", "bold");
    text(106, 101, "Spikelets", "HorizontalAlignment", "right", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", colCodes_spklt(3,:), "FontWeight", "bold");
    % legend lines: spikes vs. spikelets vs. both
    nSpkOnly = sum(percentSpikelets(1,:)==0);
    nBoth = sum(percentSpikelets(1,:)>0 & percentSpikelets(1,:)<100);
    line([nSpkOnly nSpkOnly]+0.5,ylim,"Color","black", "LineStyle", ":")
	line([nSpkOnly nSpkOnly]+0.5+nBoth,ylim,"Color","black", "LineStyle", ":")
 
% xcorr example
subplot(33,10, [133 135+20])
xcorr_example = reshape(d_2(3e4:end,29:30),[],1);
    plot(xcorr_example, "Color", "black")
    xlim([68 200]*10^3); ylim([min(xcorr_example)-1 max(xcorr_example)+1])
    axis off
    % scale bar
    line([81000 91000],[-5 -5],"Color","black","LineWidth",width_scalebar);
    line([81000 81000],[-5 15],"Color","black","LineWidth",width_scalebar);
	text(86000, -6, "200ms", "VerticalAlignment", "top", "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");
	text(81000, 16, "20mV ", "VerticalAlignment", "bottom", "HorizontalAlignment", "center",...
        "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold", "Color", "black");

    
%% xcorr principle illustration
subplot(33,10, [136 137+20])
	xcorr_example = reshape(d_2(3e4:end,29:30),[],1);
    xcorr_principle = xcorr_example(96700:99008);   
    % AP thresholds are at samples 172, 893, 1602
	plot([172:893]-1000, xcorr_principle(172:893), "Color", colCodes_spklt(3,:), "LineWidth", width_traces_zoom); hold on
    plot(1:893, xcorr_principle(1:893), "Color", colCodes_spklt(3,:), "LineWidth", width_traces_zoom);
    plot(893:1602, xcorr_principle(893:1602), "Color", colCodes_spklt(2,:), "LineWidth", width_traces_zoom); 
    plot(1602:2309, xcorr_principle(1602:2309), "Color", colCodes_spklt(3,:), "LineWidth", width_traces_zoom); 
    xlim([172-1000 2800]); ylim([-150 max(xcorr_principle)+1]);
    % discontinuity trace
    line([-148 -62],[-66.2 -60.2], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([-40 42],[-66.5 -60.5], "Color", "black");%, "LineWidth", width_traces_zoom)
    % -lag continuous
    line([172 172], [-70 -73], "Color", "black");%, "LineWidth", width_traces_zoom)
	line([893 893], [-70 -73], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([172 893],[-73 -73], "Color", "black")
	text(mean([172 893]), -74, "-Lag", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    % - lag discontinuous
    line([172 172]-1000, [-102 -105], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([893 893], [-102 -105], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([172-1000 -19], [-105 -105], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([83 893], [-105 -105], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([-62 24], [-108 -102], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([40 126], [-108 -102], "Color", "black");%, "LineWidth", width_traces_zoom)
    text(0, -108, "-Lag", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
    % +lag continuous
    line([893 893], [-112 -115], "Color", "black");%, "LineWidth", width_traces_zoom)
	line([1602 1602], [-112 -115], "Color", "black");%, "LineWidth", width_traces_zoom)
    line([893 1602],[-115 -115], "Color", "black")
	text(mean([893 1602]), -116, "+Lag", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
         "FontSize", font_size_large, "FontName", "Arial", "FontWeight", "bold")
	axis off 
    
    
%% cross-correlation
lag_all = [];  
cell_count = 0;
ap_count = 0;
lag_count = 0;
for iter = 1:length(AP_ts_inter) % iterate pF+ cells
    if length(unique(AP_max_collect{iter} > -20)) == 2 % i.e. choose cell only if spikes and spikelets in recording
        cell_count = cell_count + 1;
        lag_per_cell = [];
        for jter = 1:length(AP_ts_inter{iter}) % iterate sweeps in cell
        	if length(unique(AP_max_inter{iter}{jter} > -20)) == 2 % choose sweep only if it contains spikelets and spikes
                bool_ind_spks = AP_max_inter{iter}{jter} > -20;
                [~, lag] = crosscorrelation(AP_ts_inter{iter}{jter}(~bool_ind_spks),...
                                            AP_ts_inter{iter}{jter}(bool_ind_spks),10);
                lag_per_cell = [lag_per_cell lag];
                ap_count = ap_count + length(AP_max_inter{iter}{jter});
                lag_count = lag_count + length(lag);
            end
        end
        lag_all(cell_count,:) = histcounts(lag_per_cell, "BinWidth", 100, "Normalization", "count", "BinLimits", [-3000 3000]);
    end
end
lag_all_norm = normalize(lag_all, 2, "range", [0 1]);
[mean_xcorr, std_xcorr, sem_xcorr, median_xcorr] = statistics(lag_all_norm, 1);

subplot(33,10,[138 140+20]);
    b = bar(0.5:1:59.5, mean_xcorr, "BarWidth", 1, "FaceAlpha", 1); hold on
    set(b, "FaceAlpha", 0.65);
    plot(0.5:1:59.5, mean_xcorr, "o", "MarkerFaceColor", "black", "MarkerEdgeColor", "black", "MarkerSize", 1)
    plot([0.5:1:59.5; 0.5:1:59.5], [mean_xcorr-sem_xcorr; mean_xcorr+sem_xcorr],...
         "LineStyle", "-", "Color", "black", "LineWidth", 0.75)
    ylabel("Density"); xlabel("Lag [ms]")
    set(gca, "XTick", 0:10:60, "XTickLabel", -3000:1000:3000);
    plotAesthetics(gca,1,font_size_large);
    ylim([0 max(mean_xcorr)+0.1]); xlim([-0.5 60.5]);
    % legend
    text(60, 0.9, [num2str(lag_count) ' pairs of' newline num2str(ap_count) ' spike(let)s' newline 'in ' num2str(cell_count) ' cells'],...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontName', 'Arial', "FontWeight", "bold")
        
    
%% save figure
set(gcf,'renderer','Painters')
% print('F4_persistentFiring.pdf','-dpdf')