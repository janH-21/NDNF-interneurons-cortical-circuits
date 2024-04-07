clear all; close all; clc;
cd('/main/')
table_PV = load('./analysis/currentSteps/resultTable_PV.mat');
table_S5E2 = load('./analysis/currentSteps/resultTable_S5E2.mat');
table_INs = load('./analysis/currentSteps/resultTable_INs.mat');
table_VIP = load('./analysis/currentSteps/resultTable_VIP.mat');
table_PC = load('./analysis/currentSteps/resultTable_PC.mat');

resultTable = [table_PV.resultTable_currentSteps;...
               table_S5E2.resultTable_currentSteps;...
               table_INs.resultTable_currentSteps;...
               table_VIP.resultTable_currentSteps;...
               table_PC.resultTable_currentSteps];
               
%%
resultTable.group(resultTable.genetics == "PV-cre") = "PV";
resultTable.group(resultTable.genetics == "S5E2") = "PV";
resultTable.group(resultTable.genetics == "none/mDlx" & resultTable.maxFR_Hz_evoked > 150) = "PV";
resultTable.group(resultTable.genetics == "none/mDlx" & resultTable.maxFR_Hz_evoked <= 150) = "non-FS";
resultTable.group(resultTable.genetics == "VIP-cre") = "non-FS";
resultTable.group(resultTable.genetics == "PC") = "PC";

%%
save("resultTable_pF_groups.mat", "resultTable")