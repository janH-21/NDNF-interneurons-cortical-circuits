%% prepare workspace and load data
close all; clc;
cd('/main/');
load('./analysis/misc/aesthetics.mat')
load('./analysis/currentSteps/resultTable_ndnf.mat');

% figure
F_S10_BF_iprops = figure("Units", "centimeters", "OuterPosition", [1 1 21 29.7], "Resize", "off");

%% ----------------------------format data---------------------------------
useBF = resultTable_currentSteps.useBFanalysis;
cellType = resultTable_currentSteps.cellType;
isBF = resultTable_currentSteps.isBF(cleanDataInd & cellType == 1 & useBF == 1);
allNDNF_dats = collectData_all(cleanDataInd & cellType == 1 &  useBF == 1,:);
NDNF_dats_pF = allNDNF_dats(isBF==1,:);
NDNF_dats_non_pF = allNDNF_dats(isBF==0,:);
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
                ones(size(allNDNF_dats(isBF == 1,:),1),1)*2;...
                ones(size(allNDNF_dats(isBF == 0,:),1),1)*3];

% -------------------------- plot figure ----------------------------------
clear sp; 
for iter = 1:14
    sp{iter} = subplot(7,2,iter);
    bp = boxplot([allNDNF_dats(:,iter);...
                  NDNF_dats_pF(:,iter);...
                  NDNF_dats_non_pF(:,iter)],...
                  group_Clusts,...
                  'Orientation', 'horizontal',...
                  'PlotStyle','traditional');
    set(gca, 'ytick', 1:4,...
             'YTickLabel', {'Ndnf all ','Ndnf pF+','Ndnf pF-'})
    set(bp, "Color", "black")
    xlabel(varNames(iter))
    plotAesthetics(gca, 1, font_size_large); grid off
    h = findobj(gca,'Tag','Box');
    patch(get(h(3),'XData'),get(h(3),'YData'),colCodes_iprop_box(1,:),'FaceAlpha',.6);
    patch(get(h(2),'XData'),get(h(2),'YData'),colCodes_iprop_box(2,:),'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'),colCodes_iprop_box(3,:),'FaceAlpha',.6);   
	h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    hold on
end

% ------------------------------statistics---------------------------------
% MANOVA for group differences for all groups
dats_MANOVA_all = [allNDNF_dats;...
                   allNDNF_dats(isBF == 1,:);...
                   allNDNF_dats(isBF == 0,:)];
           
group_MANOVA_all = [ones(size(allNDNF_dats,1),1);...
                    ones(size(allNDNF_dats(isBF == 1,:),1),1)*2;...
                    ones(size(allNDNF_dats(isBF == 0,:),1),1)*3];
 
[d_all,p_all,stats_all] = manova1(dats_MANOVA_all,group_MANOVA_all);
p_all

% post-hoc MANOVAs
group_1 = NaN(3,1);
group_2 = NaN(3,1);
d = NaN(3,1);
p = NaN(3,1);
counter = 1;
for iter = 1:3
    for jter = 1:3
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

% post-hoc kruskal-wallis tests
paramStats = NaN(6*14, 6); % [param group_1 group_2 p padj padj<0.05]
for iter = 1:14
    [p,tbl,stats] = kruskalwallis(dats_MANOVA_all(:, iter), group_MANOVA_all, 'off');
    if p < 0.05
        mc = multcompare(stats, "CType", "dunn-sidak", "Display", "off");
        axes(sp{iter});
        xLims = xlim;
        counter = 1;
        for jter = 1:size(mc,1)
            if mc(jter,end) < 0.05
                line([xLims(2) xLims(2)]+0.05*(counter-1)*diff(xLims),...
                     [mc(jter,1) mc(jter,2)],"Color", "black", "LineWidth", width_scalebar)
                if mc(jter,end) < 0.001; stars = '***';
                elseif mc(jter,end) < 0.01; stars = '**';
                elseif mc(jter,end) < 0.05; stars = '*';
                end
                text(xLims(2)+0.05*(counter-1)*diff(xLims), mean(mc(jter,1:2)), stars,...
                    "HorizontalAlignment", "center", "VerticalAlignment", "top",...
                    "FontSize", 12, "FontName", "Arial", "FontWeight", "bold", "Rotation", 90)
                counter = counter + 1;
            end
        end
        xlim([xLims(1) xLims(2)+0.05*(counter)*diff(xLims)])
    end
end


%% save figure
set(gcf,'renderer','Painters')
%print('FS10_BF_iprops.pdf','-dpdf')