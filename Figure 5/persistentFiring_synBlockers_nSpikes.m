%% prepare
clear all; close all; clc
cd("/main/")
load("./analysis/misc/aesthetics.mat")


%% load data
results_pharma = load("./analysis/currentSteps/resultTable_Ndnf_synBlockers.mat");
results_control = load("./analysis/currentSteps/resultTable_ndnf.mat");
pharma = load("./analysis/Figure 5/synaptic_steps_BFstats.mat");
control = load("./analysis/Figure 5/control_steps_BFstats.mat");
% create indices
control_ind = results_control.resultTable_currentSteps.isBF == 1 &...
              results_control.resultTable_currentSteps.cellType == 1;
pharma_ind = results_pharma.resultTable_currentSteps.isBF == 1 &...
             results_pharma.resultTable_currentSteps.cellType == 1;
% assign box plot groups
box_groups = [ones(size(pharma.results,1),1);...
              ones(size(control.results,1),1)*2];
box_data_APs = [pharma.results(:,2);...
                control.results(:,2)];
box_data_BFmaxAPs = [results_pharma.resultTable_currentSteps.nEctopicSupraMax(pharma_ind)
                     results_control.resultTable_currentSteps.nEctopicSupraMax(control_ind)];
% test for normality  
normality = ~[~[adtest(box_data_APs(box_groups==1)), adtest(box_data_APs(box_groups==2));...
                adtest(box_data_BFmaxAPs(box_groups==1)), adtest(box_data_BFmaxAPs(box_groups==2))] & ...
              ~[lillietest(box_data_APs(box_groups==1)), lillietest(box_data_APs(box_groups==2));...
                lillietest(box_data_BFmaxAPs(box_groups==1)), lillietest(box_data_BFmaxAPs(box_groups==2))] &...
              ~[kstest(normalize(box_data_APs(box_groups==1))), kstest(normalize(box_data_APs(box_groups==2)));...
                kstest(normalize(box_data_BFmaxAPs(box_groups==1))), kstest(normalize(box_data_BFmaxAPs(box_groups==2)))]]
          
% statistics: APs until persistent firing     
p_APs =  ranksum(box_data_APs(box_groups==1),...
                 box_data_APs(box_groups==2))
% statistics: max # APs during persistent firing episode
[h,p_BFmaxAPs] = ttest2(box_data_BFmaxAPs(box_groups==1),...
                        box_data_BFmaxAPs(box_groups==2),...
                        "Tail", "both")

                    
%% plot results
% relabel box group for plotting           
box_groups(box_groups==2) = 0;
% plot results               
figure;
subplot(1,2,1)
    bp1 = boxplot(box_data_APs, box_groups);
    set(bp1, "Color", "black", "LineWidth", 1)
	plotAesthetics(gca, 1, 12);
    set(gca, 'XTickLabel', {'Control','Syn. blockers'})
    xlabel('Step #'); ylabel('APs until persistent firing')
    h = findobj(gca,'Tag','Box');
    patch(get(h(2),'XData'),get(h(2),'YData'), blue_std ,'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'), blue_std,'FaceAlpha',.6);
    h = findobj('LineStyle','--'); set(h, 'LineStyle','-');
    subtitle(['Wilcoxon rank-sum test p = ', num2str(p_APs)])
subplot(1,2,2)
    bp1 = boxplot(box_data_BFmaxAPs, box_groups);
    set(bp1, "Color", "black", "LineWidth", 1)
	plotAesthetics(gca, 1, 12);
    set(gca, 'XTickLabel', {'Control','Syn. blockers'})
    xlabel('Step #'); ylabel('Max. AP# during persistent firing')
    h = findobj(gca,'Tag','Box');
    patch(get(h(2),'XData'),get(h(2),'YData'), blue_std ,'FaceAlpha',.6);
    patch(get(h(1),'XData'),get(h(1),'YData'), blue_std,'FaceAlpha',.6);
    h = findobj('LineStyle','--'); set(h, 'LineStyle','-');
    subtitle(['Two-tailed t-test p = ', num2str(p_BFmaxAPs)])
    
    