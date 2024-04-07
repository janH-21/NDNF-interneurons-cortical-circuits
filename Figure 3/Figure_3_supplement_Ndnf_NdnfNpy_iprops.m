%% prepare workspace and load data
close all; clc;
cd('/main/');
load('./analysis/misc/aesthetics.mat')

% figure
FS6 = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");

% ----------------------- collect data ------------------------------------
cellType = resultTable_currentSteps.cellType(cleanDataInd);
ClustIDs =  resultTable_currentSteps.ClustIDs(~isnan(resultTable_currentSteps.ClustIDs));
allNDNF_dats = collectData(cellType == 1,:);
allNDNFNPY_dats = collectData(cellType == 2,:);
clu1_NDNF_dats = allNDNF_dats(ClustIDs == 1,:);
clu2_NDNF_dats = allNDNF_dats(ClustIDs == 2,:);
varNames = {'Spike onset [ms]',...
         'AP amplitude [mV]',...
         'AP threshold [mV]',...
         'AP half-width [ms]',...
         'AP max. slope [mV/ms]',...
         'AHP amplitude [mV]',...
         'AHP half-width [ms]',...
         'ISI9/ISI2',...
         'ISI var',...
         'V_m rest [mV]',...
         'R_{in} [MOhm]',...
         'Tau [ms]',...
         'Depol. hump amplitude [%]',...
         'Sag amplitude [%]'};
group_Clusts = [ones(size(allNDNF_dats,1),1);...
                ones(size(allNDNF_dats(ClustIDs == 1,:),1),1)*2;...
                ones(size(allNDNF_dats(ClustIDs == 2,:),1),1)*3;...
                ones(size(allNDNFNPY_dats,1),1)*4];

% ------------------------- plot figure -----------------------------------
clear sp; close all
fSF22 = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");
orange = [252, 157, 3]./255;
turquois = [42, 162, 199]./255;
green = [103, 194, 58]./255;
for iter = 1:14
    sp{iter} = subplot(7,2,iter); 
    bp = boxplot([allNDNF_dats(:,iter);...
                  clu1_NDNF_dats(:,iter);...
                  clu2_NDNF_dats(:,iter);...
                  allNDNFNPY_dats(:,iter)],...
                  group_Clusts,...
                  'Orientation', 'horizontal',...
                  'PlotStyle','traditional'); hold on;
    set(gca, 'ytick', 1:4,...
             'YTickLabel', {'all Ndnf+','Ndnf+ Clu1','Ndnf+ Clu2','Ndnf+/NPY+'})
    set(bp, "Color", "black")
    xlabel(varNames(iter))
    plotAesthetics(gca, 1, font_size_large); grid off
    h = findobj(gca,'Tag','Box');
    patch(get(h(4),'XData'),get(h(4),'YData'),colCodes_Clu(3,:),'FaceAlpha',.6);
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_Clu(1,:),'FaceAlpha',.8);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_Clu(2,:),'FaceAlpha',.8);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_Clu(4,:),'FaceAlpha',.8);
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    hold on
end

%% MANOVA for group differences for all groups
dats_MANOVA_all = [allNDNF_dats;...
                   allNDNF_dats(ClustIDs == 1,:);...
                   allNDNF_dats(ClustIDs == 2,:);...
                   allNDNFNPY_dats];
           
group_MANOVA_all = [ones(size(allNDNF_dats,1),1);...
                    ones(size(allNDNF_dats(ClustIDs == 1,:),1),1)*2;...
                    ones(size(allNDNF_dats(ClustIDs == 2,:),1),1)*3;...
                    ones(size(allNDNFNPY_dats,1),1)*4];
 
[d_all,p_all,stats_all] = manova1(dats_MANOVA_all,group_MANOVA_all);

%% post-hoc MANOVAs
group_1 = NaN(6,1);
group_2 = NaN(6,1);
d = NaN(6,1);
p = NaN(6,1);
counter = 1;
for iter = 1:4
    for jter = 1:4
        if iter<jter
            group_1(counter) = iter;
            group_2(counter) = jter;
            [d(counter), p(counter), ~] = manova1(...
                dats_MANOVA_all(group_MANOVA_all == iter | group_MANOVA_all == jter, :),...
                group_MANOVA_all(group_MANOVA_all == iter | group_MANOVA_all == jter));
            counter = counter + 1;
        end
    end
end

% multiple test correction
pAdj = p.*length(p); % Bonferroni
[group_1 group_2 d p pAdj]

%% post-hoc kruskal-wallis tests
paramStats = NaN(6*14, 6); % [param group_1 group_2 p padj padj<0.05]
for iter = 1:14
    [p,tbl,stats] = kruskalwallis(dats_MANOVA_all(:, iter), group_MANOVA_all, 'off');
    mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
	axes(sp{iter});
    xLims = xlim;
    counter = 1;
    for jter = 1:size(mc,1)
        if mc(jter,end) < 0.05
            line([xLims(2) xLims(2)]+0.06*(counter-1)*diff(xLims),...
                 [mc(jter,1) mc(jter,2)],"Color", "black", "LineWidth", 1)
            if mc(jter,end) < 0.001; stars = '***';
            elseif mc(jter,end) < 0.01; stars = '**';
            elseif mc(jter,end) < 0.05; stars = '*';
            end
            text(xLims(2)+0.06*(counter-1)*diff(xLims), mean(mc(jter,1:2)), stars,...
                "HorizontalAlignment", "center", "VerticalAlignment", "top",...
                "FontSize", 12, "FontName", "Arial", "FontWeight", "bold", "Rotation", 90)
            counter = counter + 1;
        end
    end
    xlim([xLims(1) xLims(2)+0.05*(counter)*diff(xLims)])
end

%% save figure
set(gcf,'renderer','Painters')
% print('FS6_iprops.pdf','-dpdf')