%% ---------------------- prepare workspace -------------------------------
close all; clc;
cd("/main/")
load('./analysis/misc/aesthetics.mat')

% figure
FS11_VmRest = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");
    

%% ------------------------- example recording ----------------------------
ts = (1:length(reshape(allData(37).recordingData.rawData_mV(:,:),[],1)))*20/1000;
subplot(7,5,1:5);
    % data
    plot(ts, reshape(allData(37).recordingData.rawData_mV(:,:),[],1), "Color", "black")
    % aesthetics
    ylim([-85 30]);
    xLims = xlim; xlim([19084 xLims(2)]);
    axis off;
    % scale
    line([1.105e5 1.105e5]-3*28000, [-20 -10], "Color", "black", "LineWidth", width_scalebar)
    line([1.105e5 1.125e5]-3*28000, [-20 -20], "Color", "black", "LineWidth", width_scalebar)
    text(1.105e5-3*28000, -15, "10mV ", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "middle", "FontWeight", "bold");
    text(mean([1.105e5 1.125e5])-3*28000, -21, "2s", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "center", "VerticalAlignment", "top", "FontWeight", "bold");
    pos = get(gca, "Position");
	set(gca, "Position", [pos(1) pos(2)-.015 pos(3) pos(4)]);

    
%% ------------------------- measure Vm + AHP -----------------------------
Vm_base = []; min_base = [];
Vm_pre = []; min_pre = [];
Vm_post = []; min_post = [];
collectTraces = [];
counter = 1;
isBF = resultTable_currentSteps.isBF;
% iterate data
for iter = 1:length(allData)
    if isBF(iter) == 1  
        % only choose those recordings with min. 3 sweeps after pF
        firstBF = find(allData(iter).analysis.nBarrageSpikes > 3, 1);
        if size(allData(iter).recordingData.rawData_mV,2) - firstBF > 3
            % only chose those recordings with singular pF episode in region of interest
            if sum(allData(iter).analysis.nBarrageSpikes(firstBF+1:firstBF+3)) == 0
                % find first sweep with APs
                supra = allData(iter).analysis.supraThresholdSweep;
                baseTraces(counter,:) =...
                    reshape(allData(iter).recordingData.rawData_mV(:,supra:supra+3),[],1);          
                % find pF sweep
                collectTraces(counter,:) =...
                    reshape(allData(iter).recordingData.rawData_mV(:,firstBF-6:firstBF+4),[],1);


                Vm_pre(counter,:) = [mean(collectTraces(counter,2.2e5:2.5e5)),...
                                     mean(collectTraces(counter,4.2e5:4.5e5)),...
                                     mean(collectTraces(counter,6.2e5:6.5e5)),...
                                     mean(collectTraces(counter,8.2e5:8.5e5)),...
                                     mean(collectTraces(counter,10.2e5:10.5e5)),...
                                     mean(collectTraces(counter,12.2e5:12.5e5))];
                                 
                Vm_post(counter,:) = [mean(collectTraces(counter,16.2e5:16.5e5)),...
                                      mean(collectTraces(counter,18.2e5:18.5e5)),...
                                      mean(collectTraces(counter,20.2e5:20.5e5))];
                                  
                min_pre(counter,:) = [min(collectTraces(counter,1e5:1.2e5)),...
                                      min(collectTraces(counter,3e5:3.2e5)),...
                                      min(collectTraces(counter,5e5:5.2e5)),...
                                      min(collectTraces(counter,7e5:7.2e5)),...
                                      min(collectTraces(counter,9e5:9.2e5)),...
                                      min(collectTraces(counter,11e5:11.2e5))];
                                  
                min_post(counter,:) = [min(collectTraces(counter,15e5:15.2e5)),...
                                       min(collectTraces(counter,17e5:17.2e5)),...
                                       min(collectTraces(counter,19e5:19.2e5))];
                counter = counter + 1;
            end
        end
    end
end

% calculate mean params
mean_Vm_pre = mean(Vm_pre,2);
mean_Vm_post = mean(Vm_post,2);
mean_min_pre = mean(min_pre,2);
mean_min_post = mean(min_post,2);
dPre = mean_Vm_pre - mean_min_pre;
dPost = mean_Vm_post - mean_min_post;
% staistics: paired Wilcoxon signed rank test
p1 = signrank(mean_Vm_pre, mean_Vm_post);
p2 = signrank(dPre, dPost);
% create groups for plotting
group_pre = ones(size(mean_Vm_pre));
group_post = ones(size(mean_Vm_post))*2;
groups = [group_pre; group_post];

% create groups for plotting
groups = [ones(size(mean_Vm_post)),...
          ones(size(mean_Vm_post))*2,...
          ones(size(mean_Vm_post))*3,...
          ones(size(mean_Vm_post))*4,...
          ones(size(mean_Vm_post))*5,...
          ones(size(mean_Vm_post))*6,...
          ones(size(mean_Vm_post))*7,...
          ones(size(mean_Vm_post))*8,...
          ones(size(mean_Vm_post))*9];
% organize data
Vm_min = [(min_pre),...
          (min_post)];
Vm = [(Vm_pre),...
      (Vm_post)];
dVm = Vm - Vm_min;


%% ---------------------- plot mini time course ---------------------------
SEM_Vm = std(Vm,[],1)./sqrt(size(Vm,1));
% Vm
subplot(7,5,6);
    line([6.5 6.5],[-75 -68],"Color", "red", "LineStyle", ":", "LineWidth", width_scalebar); hold on
    plot(1:9, mean(Vm,1), "o", "MarkerSize", 3, "MarkerEdgeColor", "black", "MarkerFaceColor", "black")
    plot([1:9;1:9], [mean(Vm)-SEM_Vm; mean(Vm)+SEM_Vm], "Color", "black")
    % Aesthetics
    plotAesthetics(gca,1,font_size_large);
    set(gca, "XTick", 1:9, "XTickLabel", {"-6","-5","-4", "-3","-2","-1", "+1","+2","+3"});
    ylabel("V_m [mV]"); xlabel("Sweep #");
    % legend
    ylim([-75 -67]); yLims = ylim;
    text(6.5, yLims(2), "pF onset ", "Color", "red", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", "FontWeight", "bold")
	text(mean(xlim), yLims(1)-0.27*diff(yLims), ['n=' num2str(size(Vm,1))],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "VerticalAlignment", "top", "HorizontalAlignment", "center");
    xtickangle(90)
	pos = get(gca, "Position");
	set(gca, "Position", [pos(1)-.05 pos(2) pos(3) pos(4)]);
    
% step AHP 
subplot(7,5,7);
    line([6.5 6.5],[1 4],"Color", "red", "LineStyle", ":", "LineWidth", width_scalebar); hold on
    plot(1:9, mean(dVm,1), "o", "MarkerSize", 3, "MarkerEdgeColor", "black", "MarkerFaceColor", "black")
    SEM_dVm = std(dVm,[],1)./sqrt(size(dVm,1));
    plot([1:9;1:9], [mean(dVm)-SEM_dVm; mean(dVm)+SEM_dVm], "Color", "black");
    plotAesthetics(gca,1,font_size_large);
    set(gca, "XTick", 1:9, "XTickLabel", {"-6","-5","-4", "-3","-2","-1", "+1","+2","+3"});
    xtickangle(90)
    ylabel("Step AHP [mV]"); xlabel("Sweep #");
    ylim([1.5 4.4]); yLims = ylim;
    text(6.5, yLims(2), "pF onset ", "Color", "red", "FontSize", font_size_large, "FontName", "Arial",...
        "HorizontalAlignment", "right", "VerticalAlignment", "top", "FontWeight", "bold")
    text(mean(xlim), yLims(1)-0.27*diff(yLims), ['n=' num2str(size(Vm,1))],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "VerticalAlignment", "top", "HorizontalAlignment", "center");
    pos = get(gca, "Position");
	set(gca, "Position", [pos(1)-.025 pos(2) pos(3) pos(4)]);
    

%% -------------------------- box plots -----------------------------------
% boxplot mean values + statistics
pre_Vm_pool = reshape(Vm_pre(:,4:6),[],1);
post_Vm_pool = reshape(Vm_post,[],1);
dPre_pool = reshape(Vm_pre(:,4:6)-min_pre(:,4:6),[],1);
dPost_pool = reshape(Vm_post-min_post,[],1);
% create groups for plotting
group_pre = ones(size(pre_Vm_pool));
group_post = ones(size(post_Vm_pool))*2;
groups = [group_pre; group_post];

% box + statistics Vm
subplot(7,5,8);
    bp1 = boxplot([pre_Vm_pool; post_Vm_pool], groups); hold on;
    %plot([group_pre group_post]',[pre_Vm_pool post_Vm_pool]', "-", "Color", "black");
    set(bp1, "Color", "black");
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    set(gca, "XTickLabel", ["pre pF", "post PF"]);
    ylabel("V_m [mV]");
    h = findobj(gca,'Tag','Box');
    patch(get(h(1), 'XData'), get(h(1),'YData'),blue_std, 'FaceAlpha', .8);
    patch(get(h(2), 'XData'), get(h(2),'YData'),blue_std, 'FaceAlpha', .8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    % statistics
    p_Vm = ranksum(pre_Vm_pool, post_Vm_pool); % p = 0.0196
    yLims = ylim;   
    line([1 2],[yLims(2) yLims(2)]-0.05*diff(xlim),...
        "Color", "black", "LineWidth", width_scalebar);
    text(1.5, yLims(2)-0.05*diff(xlim), "*",...
        "FontSize", font_stat_star, "FontName", "Arial",...
        "VerticalAlignment", "middle", "HorizontalAlignment", "center");
    % n observations
    text(mean(xlim), yLims(1)-0.12*diff(yLims), ['n=' num2str(size(Vm,1))],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "VerticalAlignment", "top", "HorizontalAlignment", "center");
    pos = get(gca, "Position");
	set(gca, "Position", [pos(1) pos(2) pos(3) pos(4)]);
   
% box + statistics step AHP
subplot(7,5,9);
    bp1 = boxplot([dPre_pool; dPost_pool], groups); hold on;
    %plot([group_pre group_post]',[dPre_pool dPost_pool]', "-", "Color", "black");
    set(bp1, "Color", "black");
    % aesthetics
    plotAesthetics(gca, 1, font_size_large);
    set(gca, "XTickLabel", ["pre pF", "post PF"]);
    ylabel("Step AHP [mV]");
    h = findobj(gca,'Tag','Box');
    patch(get(h(1), 'XData'), get(h(1),'YData'),blue_std, 'FaceAlpha', .8);
    patch(get(h(2), 'XData'), get(h(2),'YData'),blue_std, 'FaceAlpha', .8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    % statistics
    p_AHP = ranksum(dPre_pool, dPost_pool); % p = 0.0146
	yLims = ylim;
    line([1 2],[yLims(2) yLims(2)]-0.05*diff(xlim),"Color", "black", "LineWidth", width_scalebar)
	text(1.5, yLims(2)-0.05*diff(xlim), "*",...
        "FontSize", font_stat_star, "FontName", "Arial",...
        "VerticalAlignment", "middle", "HorizontalAlignment", "center");
	text(mean(xlim), yLims(1)-0.12*diff(yLims), ['n=' num2str(size(Vm,1))],...
        "FontSize", font_size_large, "FontName", "Arial",...
        "VerticalAlignment", "top", "HorizontalAlignment", "center");
    pos = get(gca, "Position");
	set(gca, "Position", [pos(1)+.01 pos(2) pos(3) pos(4)]);
    
    