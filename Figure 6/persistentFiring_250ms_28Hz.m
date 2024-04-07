%% prepare workspace and load data
clear all; close all; clc;
cd("/main/");
load('./analysis/misc/aesthetics.mat')
% load 250ms 2Hz 28Hz example
[d_250ms, si, ~] = abfload('./data/persistentFiring_NDNF_250ms28Hz/2022_10_27_NdnfCre_line3016_mouse165_0043.abf');
% load 250ms 2Hz 28Hz data
load('./analysis/persistentFiring_NDNF_250ms28Hz/allData_250ms_2Hz_7APs.mat');
% figure
Fig6_28Hz = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");


%% separate by persistent firing
pF_pos = [2 3 5 6 8 12 13 14 16 17 18 19 22 24 25 26 27 28];
pF_neg = [1 4 7 9 10 11 15 20 21 23];


%% barplot: %pF+/- (barplot; needs workaround to make stacked with single columns)
n_total = length(pF_pos) + length(pF_neg);
percent_pos = length(pF_pos) / n_total * 100;
percent_neg = length(pF_neg) / n_total * 100;

subplot(6,5,1)
    rectangle("Position", [-.5 0 1 percent_pos], "FaceColor", [colCodes_bar(1,:) .6])
    rectangle("Position", [-.5 percent_pos 1 percent_neg], "FaceColor", [colCodes_bar(2,:) .6])
    xlim([-1.4 1.4])
    % aesthetics
    plotAesthetics(gca,1,font_size_large);
    set(gca, "YTick",0:50:100,"XTickLabel","");
    ylabel("Cells [%]");
   	% legend
    text(-1.31, 98, 'pF+',...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
        'FontSize', font_size_large, 'FontWeight', 'bold',...
        'FontName', 'Arial', 'Color', colCodes_bar(1,:));
    text(1.35, 98, 'pF-',...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
        'FontSize', font_size_large, "FontWeight", "bold",...
        'FontName', 'Arial', 'Color', colCodes_bar(2,:));
    % n observations
    text(0, -5, ['n=', num2str(n_total)],...
        "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_large, "FontName", "Arial")

    
%% gather data
pF_start_step = NaN(10,1);
pF_end_step = NaN(10,1);
step_AP_thresh = NaN(length(allData_250ms_2Hz_7APs),...
                     size(allData_250ms_2Hz_7APs(1).analysis.step.stepTimes_ms,1));
step_AP_MVal = NaN(length(allData_250ms_2Hz_7APs),...
                     size(allData_250ms_2Hz_7APs(1).analysis.step.stepTimes_ms,1));
inter_step_APs = NaN(length(allData_250ms_2Hz_7APs),...
                     size(allData_250ms_2Hz_7APs(1).analysis.step.stepTimes_ms,1)-1);
inter_step_Vm = NaN(length(allData_250ms_2Hz_7APs),...
                    size(allData_250ms_2Hz_7APs(1).analysis.step.stepTimes_ms,1)-1);
inter_step_AHP = NaN(length(allData_250ms_2Hz_7APs),...
                     size(allData_250ms_2Hz_7APs(1).analysis.step.stepTimes_ms,1)-1);
step_FRs = NaN(n_total, length(allData_250ms_2Hz_7APs(1).analysis.step.stepFRs_Hz));
% iterate cells
for jter = 1:length(allData_250ms_2Hz_7APs)
    step_FRs(jter,:) = allData_250ms_2Hz_7APs(jter).analysis.step.stepFRs_Hz;
    if ismember(jter, pF_pos)
        pF_start_step(jter) = find(allData_250ms_2Hz_7APs(jter).analysis.step.interStepFRs_Hz>4,1);
        pF_end_step(jter) = max(find(allData_250ms_2Hz_7APs(jter).analysis.step.interStepFRs_Hz>4));
    else
        pF_start_step(jter) = NaN;
        pF_end_step(jter) = NaN;
    end
    % iterate steps
    for iter = 1:size(allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms,1)
        % step stats
        step_inds = allData_250ms_2Hz_7APs(jter).analysis.APs.APts_ms > ...
                    allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms(iter,1) & ...
                    allData_250ms_2Hz_7APs(jter).analysis.APs.APts_ms < ...
                    allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms(iter,2);
        step_AP_thresh(jter, iter) = mean(allData_250ms_2Hz_7APs(jter).analysis.APs.APtsVals_mV(step_inds));
        step_AP_MVal(jter, iter) = mean(allData_250ms_2Hz_7APs(jter).analysis.APs.MVal_mVpms(step_inds));
        % inter-step stats
        if iter < size(allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms,1)  
            inter_step_inds = allData_250ms_2Hz_7APs(jter).analysis.APs.APts_ms > ...
                              allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms(iter,2) & ...
                              allData_250ms_2Hz_7APs(jter).analysis.APs.APts_ms < ...
                              allData_250ms_2Hz_7APs(jter).analysis.step.stepTimes_ms(iter+1,1);
            inter_step_APs(jter, iter) = length(allData_250ms_2Hz_7APs(jter).analysis.APs.APtsVals_mV(inter_step_inds));
            inter_step_Vm(jter, iter) = allData_250ms_2Hz_7APs(jter).analysis.step.interStepVm_mV(iter);
            inter_step_AHP(jter, iter) = allData_250ms_2Hz_7APs(jter).analysis.step.interStepAHP_mV(iter);
            
        end
    end
end


%% boxplot: # ectopic spikes
box_group = ones(size(inter_step_APs,1),1);
box_group(pF_pos) = 2;
interStepAPs_total = sum(inter_step_APs,2);
subplot(6,5,2)
    box_250_ectopic = boxplot(interStepAPs_total, box_group);
    % aesthetics    
    set(box_250_ectopic, 'LineWidth', 1, "Color", "black")
	plotAesthetics(gca, 1, font_size_large);
    set(gca, "XTick", 1:2, 'XTickLabel', {"pF-",'pF+'})
    h = findobj(gca,'Tag','Box');
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_bar(2,:),'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_bar(1,:),'FaceAlpha',.6);   
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    ylabel("Total ectopic AP#")
    % statistics
    ranksum(interStepAPs_total(pF_pos), interStepAPs_total(pF_neg)); % p = 1.7398e-05
    yLims = ylim; yPos = yLims(2) + diff(yLims) * 0.05;
    text(1, yPos, ['n=' num2str(length(pF_neg))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
    text(2, yPos, ['n=' num2str(length(pF_pos))],...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
    yPos = yLims(2) + diff(yLims) * 0.1;
    line([1 2], [yPos yPos], "Color", "black", "LineWidth", width_scalebar);
    text(1.5, yPos, "***",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_stat_star, "FontName", "Arial", "FontWeight", "bold")
    ylim([yLims(1) yPos])


%% --------------------- FR & threshold dynamics --------------------------
% FR
FR_peri_PF = NaN(length(pF_pos)-1, 40);
FR_first_neg = NaN(length(pF_neg), 20);
FR_last_neg = NaN(length(pF_neg), 20);
FR_first_pos = NaN(length(pF_pos)-1, 20);
FR_last_pos = NaN(length(pF_pos)-1, 20);
FR_pre_mean_pos = NaN(length(pF_pos)-1, 1);
FR_post_mean_pos = NaN(length(pF_pos)-1, 1);
% AP threshold
thr_peri_PF = NaN(length(pF_pos)-1, 40);
thr_first_neg = NaN(length(pF_neg), 20);
thr_last_neg = NaN(length(pF_neg), 20);
thr_first_pos = NaN(length(pF_pos)-1, 20);
thr_last_pos = NaN(length(pF_pos)-1, 20);
thr_pre_mean_pos = NaN(length(pF_pos)-1, 1);
thr_post_mean_pos = NaN(length(pF_pos)-1, 1);

for iter = 1:length(pF_pos)-1
    cellNo = pF_pos(iter);
    % pF+ mean FR
    FR_pre_mean_pos(iter) = mean(step_FRs(cellNo,1:pF_start_step(cellNo)),2);
    FR_post_mean_pos(iter) = mean(step_FRs(cellNo,pF_start_step(cellNo):end),2);
    % pF+ FR
    FR_peri_PF(iter,:) = step_FRs(cellNo,pF_start_step(cellNo)-19:pF_start_step(cellNo)+20);
    FR_first_pos(iter,:) = step_FRs(cellNo,1:20);
    FR_last_pos(iter,:) = step_FRs(cellNo,end-19:end);
    % pF+ threshold
	thr_peri_PF(iter,:) = step_AP_thresh(cellNo,pF_start_step(cellNo)-19:pF_start_step(cellNo)+20);
    thr_first_pos(iter,:) = step_AP_thresh(cellNo,1:20);
    thr_last_pos(iter,:) = step_AP_thresh(cellNo,end-19:end);
    % pF+ mean threshold
    thr_pre_mean_pos(iter) = mean(step_AP_thresh(cellNo,1:pF_start_step(cellNo)),2);
    thr_post_mean_pos(iter) = mean(step_AP_thresh(cellNo,pF_start_step(cellNo):end),2);
end

for iter = 1:length(pF_neg)
    % pF- AP
    cellNo = pF_neg(iter);
	FR_first_neg(iter,:) = step_FRs(cellNo,1:20);
	FR_last_neg(iter,:) = step_FRs(cellNo,end-19:end);
    % pF- thresh
	thr_first_neg(iter,:) = step_AP_thresh(cellNo,1:20);
	thr_last_neg(iter,:) = step_AP_thresh(cellNo,end-19:end);
end

pf_FR_pos_dats_long = [FR_first_pos, FR_peri_PF, FR_last_pos];
pf_FR_neg_dats_long = [FR_first_neg, FR_last_neg];
pf_thr_pos_dats_long = [thr_first_pos, thr_peri_PF, thr_last_pos];
pf_thr_neg_dats_long = [thr_first_neg, thr_last_neg];


%% plot firing rates
subplot(6,2,3)
[FR_pos_mean_long, ~, FR_pos_SEM_long, ~] = statistics(pf_FR_pos_dats_long, 1);
[FR_neg_mean_long, ~, FR_neg_SEM_long, ~] = statistics(pf_FR_neg_dats_long, 1);
    rectangle("Position", [21.5 20 20 60], "EdgeColor", "none", "FaceColor", [grey_light+.1 .3]); hold on
	rectangle("Position", [61.5 20 20 60], "EdgeColor", "none", "FaceColor", [grey_light+.1 .3]);
    line([40.5 40.5],[20 70],"Color", "red", "LineStyle", ":", "LineWidth", 1)
    plot([1:80]+.05, FR_pos_mean_long, "o", "MarkerSize", 2, "MarkerEdgeColor", colCodes_bar(1,:), "MarkerFaceColor", colCodes_bar(1,:)); hold on
	plot([1:80;1:80]+.05, [FR_pos_mean_long-FR_pos_SEM_long; FR_pos_mean_long+FR_pos_SEM_long], "Color", colCodes_bar(1,:));
	plot([1:20, 61:80]-.05, FR_neg_mean_long, "o", "MarkerSize", 2, "MarkerEdgeColor", colCodes_bar(2,:), "MarkerFaceColor", colCodes_bar(2,:));
	plot([[1:20, 61:80];[1:20, 61:80]]-.05, [FR_neg_mean_long-FR_neg_SEM_long; FR_neg_mean_long+FR_neg_SEM_long], "Color", colCodes_bar(2,:));
	% aesthetics
	plotAesthetics(gca,1,font_size_large);
    set(gca, "XTick", 0:5:80, "XTickLabel", " ")
	ylabel("Avg. FR [Hz]"); xlabel("Current injection step number");
	xlim([0 81]); ylim([20 80])
    % legend
    text(10, 80, "First 20 steps",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(30, 80, "20 steps pre pF",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(50, 80, "20 steps post pF",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(70, 80, "Last 20 steps",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(40, 35, "pF onset ",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", "red",...
        "Rotation", 90)

    
%% plot AP threshold
subplot(6,2,4)
[thr_pos_mean_long, ~, thr_pos_SEM_long, ~] = statistics(pf_thr_pos_dats_long, 1);
[thr_neg_mean_long, ~, thr_neg_SEM_long, ~] = statistics(pf_thr_neg_dats_long, 1);
    rectangle("Position", [21.5 -47 20 13], "EdgeColor", "none", "FaceColor", [grey_light+.1 .3]); hold on
	rectangle("Position", [61.5 -47 20 13], "EdgeColor", "none", "FaceColor", [grey_light+.1 .3]);
    line([40.5 40.5],[-43 -32],"Color", "red", "LineStyle", ":", "LineWidth", 1)
    plot([1:80]+.05, thr_pos_mean_long, "o", "MarkerSize", 2, "MarkerEdgeColor", colCodes_bar(1,:), "MarkerFaceColor", colCodes_bar(1,:)); hold on
	plot([1:80;1:80]+.05, [thr_pos_mean_long-thr_pos_SEM_long; thr_pos_mean_long+thr_pos_SEM_long], "Color", colCodes_bar(1,:));
	plot([1:20, 61:80]-.05, thr_neg_mean_long, "o", "MarkerSize", 2, "MarkerEdgeColor", colCodes_bar(2,:), "MarkerFaceColor", colCodes_bar(2,:));
	plot([[1:20, 61:80];[1:20, 61:80]]-.05, [thr_neg_mean_long-thr_neg_SEM_long; thr_neg_mean_long+thr_neg_SEM_long], "Color", colCodes_bar(2,:));
	% aesthetics
	plotAesthetics(gca,1,font_size_large);
	set(gca, "XTick", 0:5:80, "XTickLabel", " ")
	ylabel("Avg. AP thresh. [mV]"); xlabel("Current injection step number");
	xlim([0 81]); ylim([-47 -34])
    % legend
    text(10, -34, "First 20 steps",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(30, -34, "20 steps pre pF",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(50, -34, "20 steps post pF",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(70, -34, "Last 20 steps",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_small, "FontName", "Arial")
	text(40, -36, "pF onset ",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", "red",...
        "Rotation", 90)

    
%% include cell 28 
% which was excluded for large plots due to very early pF onset
cellNo = 28;
% pF+ mean FR
FR_pre_mean_pos(18) = mean(step_FRs(cellNo,1:pF_start_step(cellNo)),2);
FR_post_mean_pos(18) = mean(step_FRs(cellNo,pF_start_step(cellNo):end),2);
% pF+ FR
FR_peri_PF(18,:) = NaN;
FR_first_pos(18,:) = NaN;
FR_last_pos(18,:) = NaN;
FR_peri_PF(18,16:25) = step_FRs(cellNo,pF_start_step(cellNo)-4:pF_start_step(cellNo)+5);
FR_first_pos(18,1:4) = step_FRs(cellNo,1:4);
FR_last_pos(18,end-2:end) = step_FRs(cellNo,end-2:end);
% pF+ threshold
thr_peri_PF(18,:) = NaN;
thr_first_pos(18,:) = NaN;
thr_last_pos(18,:) = NaN;
thr_peri_PF(18,16:25) = step_AP_thresh(cellNo,pF_start_step(cellNo)-4:pF_start_step(cellNo)+5);
thr_first_pos(18,1:4) = step_AP_thresh(cellNo,1:4);
thr_last_pos(18,end-2:end) = step_AP_thresh(cellNo,end-2:end);
% pF+ mean threshold
thr_pre_mean_pos(18) = mean(step_AP_thresh(cellNo,1:pF_start_step(cellNo)),2);
thr_post_mean_pos(18) = mean(step_AP_thresh(cellNo,pF_start_step(cellNo):end),2);

    
%% format data short plot
pf_FR_pos_dats_short = [FR_first_pos(:,1:4), FR_pre_mean_pos, FR_peri_PF(:,16:25), FR_post_mean_pos, FR_last_pos(:,end-2:end)];
pf_FR_neg_dats_short = [FR_first_neg(:,1:4), FR_last_neg(:,end-2:end)];
pf_thr_pos_dats_short = [thr_first_pos(:,1:4), thr_pre_mean_pos, thr_peri_PF(:,16:25), thr_post_mean_pos, thr_last_pos(:,end-2:end)];
pf_thr_neg_dats_short = [thr_first_neg(:,1:4), thr_last_neg(:,end-2:end)];


%% FR dynamics mean/sem plot
subplot(6,2,5)
    % legend
    line([0 20], [28 28], "Color", [.1 .1 .1], "LineWidth", width_scalebar); hold on
    rectangle("Position", [4.5 20 1 50], "EdgeColor", [.9 .9 .9 .5], "FaceColor", [.9 .9 .9 .5])
    rectangle("Position", [15.5 20 1 50], "EdgeColor", [.9 .9 .9 .5], "FaceColor", [.9 .9 .9 .5])
    text(2.5, 18, "First 4", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(8, 18, "5 pre pF", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(13, 18, "5 post pF", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(18.5, 18, "Last 3", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(5, 71, "Avg. pre pF ", "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold", "Color", [.6 .6 .6])
	text(16, 71, "Avg. post pF", "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold", "Color", [.6 .6 .6])
	text(10, 21, "    Calibrated FR  (7APs/250ms)",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "Color", [0.1 0.1 0.1])
    line([10.5 10.5], [20 80], "Color", "red", "LineWidth", width_scalebar, "LineStyle", ":");
    text(10.5, 41, "pF onset ",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red",...
        "FontWeight", "bold", "Rotation", 90)

    % data
    [FR_pos_mean_short, ~, FR_pos_SEM_short, ~] = statistics(pf_FR_pos_dats_short, 1);
    [FR_neg_mean_short, ~, FR_neg_SEM_short, ~] = statistics(pf_FR_neg_dats_short, 1);
    plot([1:19]+.05, FR_pos_mean_short, "o", "MarkerSize", 3, "MarkerEdgeColor", colCodes_bar(1,:), "MarkerFaceColor", colCodes_bar(1,:));
    plot([1:19;1:19]+.05, [FR_pos_mean_short-FR_pos_SEM_short; FR_pos_mean_short+FR_pos_SEM_short], "Color", colCodes_bar(1,:));
    plot([1:4, 17:19]-.05, FR_neg_mean_short, "o", "MarkerSize", 3, "MarkerEdgeColor", colCodes_bar(2,:), "MarkerFaceColor", colCodes_bar(2,:));
    plot([[1:4, 17:19];[1:4, 17:19]]-.05, [FR_neg_mean_short-FR_neg_SEM_short; FR_neg_mean_short+FR_neg_SEM_short], "Color", colCodes_bar(2,:));
    % aesthetics
    plotAesthetics(gca,1,font_size_large);
    ylabel("Step FR [Hz]"); xlabel("Step #");
    ylim([20 85])
    set(gca, "XTick", 1:19,...
        "XTickLabel", ["    "],...
        "XTickLabelRotation", 90)

    
%% FR dynamics boxes
FR_first10_neg = reshape(step_FRs(pF_neg,1:10),[],1);
FR_first10_pos = reshape(step_FRs(pF_pos,1:10),[],1);
FR_last10_neg = reshape(step_FRs(pF_neg,end-9:end),[],1);
FR_last10_pos = reshape(step_FRs(pF_pos,end-9:end),[],1);
box_data_FR = [FR_first10_neg;...
               FR_last10_neg;...
               FR_first10_pos;...
               FR_last10_pos]';
box_group_FR = [ones(length(FR_first10_neg),1);...
                ones(length(FR_last10_neg),1)*2;...
                ones(length(FR_first10_pos),1)*3;...
                ones(length(FR_last10_pos),1)*4];

subplot(6,2,7)
	bp_FR = boxplot(box_data_FR, box_group_FR); hold on
    % aesthetics    
    set(bp_FR, "Color", "black")
	plotAesthetics(gca, 1, font_size_large);
    set(gca, 'XTickLabel', {'first 10','last 10', 'first 10','last 10'})
    xlabel('Step #'); ylabel('Avg. step FR [Hz]')
    h = findobj(gca,'Tag','Box');
    patch(get(h(4),'XData'),get(h(4),'YData'),colCodes_bar(2,:),'FaceAlpha',.6);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_bar(2,:),'FaceAlpha',.6); 
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_bar(1,:),'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_bar(1,:),'FaceAlpha',.6);   
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    % statistics
    [p, ~, stats] = kruskalwallis(box_data_FR, box_group_FR, "off");
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    yLims = ylim;
    counter = 1;
    if p < 0.05
        for iter = 1:size(mc,1)
            if mc(iter, end) < 0.05
                yPos = yLims(2) + diff(yLims) * 0.15 * counter;
                line([mc(iter, 1) mc(iter, 2)],[yPos yPos],"Color", "black", "LineWidth", width_scalebar);
                if mc(iter, end) < 0.001; stars = "***";
                elseif mc(iter, end) < 0.01; stars = "**";
                else; stars = "*"; end
                text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                    "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                    "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
                counter = counter + 1;
            end
        end
    end
    ylim([yLims(1) yPos])
    xtickangle(35)

%% threshold dynamics mean/sem plot
subplot(6,2,6)
    % legend
    line([10.5 10.5], [-43 27], "Color", "red", "LineWidth", width_scalebar, "LineStyle", ":"); hold on 
    rectangle("Position", [4.5 -45 1 16], "EdgeColor", [.9 .9 .9 .5], "FaceColor", [.9 .9 .9 .5])
    rectangle("Position", [15.5 -45 1 16], "EdgeColor", [.9 .9 .9 .5], "FaceColor", [.9 .9 .9 .5])
    text(2.5, -46, "First 4", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(8, -46, "5 pre pF", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(13, -46, "5 post pF", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(18.5, -46, "Last 3", "HorizontalAlignment", "center", "VerticalAlignment", "top",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold")
	text(5, -29, "Avg. pre pF ", "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold", "Color", [.6 .6 .6])
	text(16, -29, "Avg. post pF", "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_small, "FontName", "Arial", "FontWeight", "bold", "Color", [.6 .6 .6])
    text(10.5, -32.2, "pF onset ",...
        "HorizontalAlignment", "center", "VerticalAlignment", "bottom",...
        "FontSize", font_size_large, "FontName", "Arial", "Color", "red",...
        "FontWeight", "bold", "Rotation", 90)

    % data
    [thr_pos_mean_short, ~, thr_pos_SEM_short, ~] = statistics(pf_thr_pos_dats_short, 1);
    [thr_neg_mean_short, ~, thr_neg_SEM_short, ~] = statistics(pf_thr_neg_dats_short, 1);
    plot([1:19]+.05, thr_pos_mean_short, "o", "MarkerSize", 3, "MarkerEdgeColor", colCodes_bar(1,:), "MarkerFaceColor", colCodes_bar(1,:)); hold on
    plot([1:19;1:19]+.05, [thr_pos_mean_short-thr_pos_SEM_short; thr_pos_mean_short+thr_pos_SEM_short], "Color", colCodes_bar(1,:));
    plot([1:4, 17:19]-.05, thr_neg_mean_short, "o", "MarkerSize", 3, "MarkerEdgeColor", colCodes_bar(2,:), "MarkerFaceColor", colCodes_bar(2,:));
    plot([[1:4, 17:19];[1:4, 17:19]]-.05, [thr_neg_mean_short-thr_neg_SEM_short; thr_neg_mean_short+thr_neg_SEM_short], "Color", colCodes_bar(2,:));
    % aesthetics
    plotAesthetics(gca,1,font_size_large);
    ylabel("AP threshold [mV]"); xlabel("Step #");
	ylim([-47 -34])
    set(gca, "XTick", 1:19,...
        "XTickLabel", ["     "],...
        "XTickLabelRotation", 90)
    
    
%% AP thresh boxes
step_AP_thresh_norm = 1./(step_AP_thresh./step_AP_thresh(:,1));
thresh_first10_neg = reshape(step_AP_thresh_norm(pF_neg,1:10),[],1);
thresh_first10_pos = reshape(step_AP_thresh_norm(pF_pos,1:10),[],1);
thresh_last10_neg = reshape(step_AP_thresh_norm(pF_neg,end-9:end),[],1);
thresh_last10_pos = reshape(step_AP_thresh_norm(pF_pos,end-9:end),[],1);
box_data_thresh = [thresh_first10_neg;...
                   thresh_last10_neg;...
                   thresh_first10_pos;...
                   thresh_last10_pos]';
box_group_thresh = [ones(length(thresh_first10_neg),1);...
                    ones(length(thresh_last10_neg),1)*2;...
                    ones(length(thresh_first10_pos),1)*3;...
                    ones(length(thresh_last10_pos),1)*4];
subplot(6,2,8)
    bp_thresh = boxplot(box_data_thresh, box_group_thresh); hold on
	% aesthetics    
    set(bp_thresh, "Color", "black")
	plotAesthetics(gca, 1, font_size_large);
    set(gca, 'XTickLabel', {'first 10','last 10', 'first 10','last 10'})
    xlabel('Step #'); ylabel('Norm. avg. AP thresh. [mV]')
    h = findobj(gca,'Tag','Box');
    patch(get(h(4),'XData'),get(h(4),'YData'),colCodes_bar(2,:),'FaceAlpha',.6);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_bar(2,:),'FaceAlpha',.6); 
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_bar(1,:),'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_bar(1,:),'FaceAlpha',.6);   
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    % statistics
    [p, ~, stats] = kruskalwallis(box_data_thresh, box_group_thresh, "off");
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
    yLims = ylim;
    counter = 1;
    if p < 0.05
        for iter = 1:size(mc,1)
            if mc(iter, end) < 0.05
                yPos = yLims(2) + diff(yLims) * 0.15 * counter;
                line([mc(iter, 1) mc(iter, 2)],[yPos yPos],"Color", "black", "LineWidth", width_scalebar);
                if mc(iter, end) < 0.001; stars = "***";
                elseif mc(iter, end) < 0.01; stars = "**";
                else; stars = "*"; end
                text(mean([mc(iter, 1) mc(iter, 2)]), yPos, stars,...
                    "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
                    "FontName", "Arial", "FontSize", font_stat_star, "FontWeight", "bold")
                counter = counter + 1;
            end
        end
    end
    ylim([yLims(1) yPos])
    xtickangle(35)

    
%% persistent firing onset distribution
pf_onset_step = pF_start_step;
subplot(22,1,16)
    bp_onset = boxplot(pf_onset_step, "Orientation", "horizontal"); hold on
    set(bp_onset, "Color", "black");
	h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),blue_std,'FaceAlpha',.8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
	scatter(pf_onset_step(pF_pos), ones(length(pf_onset_step(pF_pos)),1)+.25, 'MarkerEdgeColor', 'black',"LineWidth", 1, "SizeData", 6, 'MarkerFaceColor', 'black');
	xlabel("Current injection step number of peristent firing onset"); 
    set(gca, "YTickLabel", "");
    xlim([0, 191])
    plotAesthetics(gca, 1, font_size_large);
	text(mean(xlim), 1.5, "n=5",...
        "HorizontalAlignment", "center", "VerticalAlignment", "middle",...
        "FontSize", font_size_large, "FontName", "Arial")
    
    
%% example recording
d_plot = squeeze(reshape(d_250ms(:,1,1:10),[],1));
subplot(6,8,[41 48])
    plot(d_plot(3600:52000), "Color", "black"); hold on
    plot(54000:54000+diff([1001440 1304080]), d_plot(1001440:1304080), "Color", "black")
   	plot((54000:54000+diff([4927550 4975950]))+5000+diff([1001440 1304080]), d_plot(4927550:4975950), "Color", "black")
    xlim([-20000 402935+20000]); ylim([-85 10])
    line([46500 50500], [-85 -77], "Color", "black")
    line([46500 50500]+(54000-48500), [-85 -77], "Color", "black")
    line([356565-1000 356565+1000]+1000, [-85 -77], "Color", "black")
    line([358744-1000 358744+1000]+3500, [-82 -74], "Color", "black")
    axis off
    text(20817, -83, 'First steps',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontWeight', 'bold',...
        'FontName', 'Arial');
    text(202054, -83, 'peri-pF steps',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontWeight', 'bold',...
        'FontName', 'Arial');
    text(383680, -83, 'Last steps',...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
        'FontSize', font_size_large, 'FontWeight', 'bold',...
        'FontName', 'Arial');