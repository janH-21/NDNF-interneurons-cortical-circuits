%% -------------------- prepare workspace & data ----------------------- %%
clear all; close all; clc;
 
% CHOOSE & LOAD FILE  
cd("/main/")
[file, path] = uigetfile({'*.abf'});  
[d,si,h]=abfload(fullfile(path,file)); size(d)
% d - (datapoint, channel, sweep), si - sampling interval in µs, h - meta

% prepare original data
dats = reshape(d(:,1,:), size(d,1), size(d,3)); % reshape data
t = (1:size(d,1))* si/1000; % time stamps in ms
ind = 1:size(dats,1); % index points
 
% set wd to dir of this script for relative file paths to work
cd(fileparts(which(mfilename('fullpath'))))


%% -----------------------USER EDITABLE BLOCK--------------------------- %%
stepSize = 10; % needs to be adapted for some older recordings of non-NDNF cells

% find zero pA injection step and calculate current steps
[~, minInd] = min(abs(1-mean(d(1300*1000/20:1800*1000/20, 1, :))./mean(d(300*1000/20:800*1000/20, 1, :))));
startValue = (minInd - 1) * -10; 
I = startValue:stepSize:(size(d,3)*stepSize-startValue); % current steps in pA
% times
stepTime = [1062, 2062]; % start and end of current injection
interSweepTime = 0; % [ms]

% Choose cell type. Important for analysis scripts downstream. 
% cellType = "L23_non-FS_IN";
% cellType = "L23_VIP";
% cellType = "L23_PC"; 
% cellType = "L23_FS_IN";
cellType = "NDNF+";
% cellType = "NDNF+/NPY+";

% save result under pathOut
pathOut = path;

%% --------------------END USER EDITABLE BLOCK-------------------------- %%
%% -------------------------- get stats -------------------------------- %%
% detect APs parameters
n = 1:1:size(d,3); % sweep numbers to be analysed
cutOffmVpms = 20;  
tAHPwindow = 8000; 
tSpikeMaxWindow = 1000;
tStepWindow = [1062, 2062+2];
% get APstats
plotXLim = [0 4000];
APstats = barrage_detectAPs(dats, si, t, n, cutOffmVpms, tStepWindow, tAHPwindow, tSpikeMaxWindow, plotXLim);
                                                      
% calculate time lines
[spikeTimeLine, sweepTimeLine, spikeTimes, stepTimes, interSweepTimes] = ...
    extractSpikeAndSweepTimes(stepTime, APstats.threshTs, t, interSweepTime);

% save data
% ID information
saveData = struct();
saveData.recordingFileName = file;
saveData.cellType = cellType;

% raw redording Data
saveData.recordingData.rawData_mV = dats;
saveData.recordingData.si_us = si;
saveData.recordingData.stepSize_pA = stepSize;
saveData.recordingData.Isteps_pA = I;
saveData.recordingData.sweepDuration = max(t);
saveData.recordingData.interSweepTime = interSweepTime;

% analyzed Data
saveData.analysis.stepTimes_ms = stepTimes;
saveData.analysis.interSweepTimes_ms = interSweepTimes;
saveData.analysis.spikeTimes_ms.cell = spikeTimes;
saveData.analysis.spikeTimes_ms.mat = cell2mat(spikeTimes); 
saveData.analysis.FR_Hz_sweep = cell2mat(APstats.nSpikesTotal)./(stepTime(2)-stepTime(1))*1000;
saveData.analysis.FR_Hz_step = cell2mat(APstats.nSpikesDuringStep)./(stepTime(2)-stepTime(1))*1000;
saveData.analysis.APstats = APstats;
saveData.analysis.spikeTimeLine_ms = spikeTimeLine;
saveData.analysis.sweepTimeLine_ms = sweepTimeLine;
saveData.analysis.suitableForBFanalysis = BFanalysis;

% write to file
fOut = strcat(pathOut, file(1:length(file)-4), '_barrageDataSingleCell.mat');

step_FR_trajectory = [APstats.nSpikesDuringStep{:}]
max_step_FR = max([APstats.nSpikesDuringStep{:}])
max_base_FR = max([APstats.nSpikesTotal{:}] - [APstats.nSpikesDuringStep{:}])
plot(d(:,1,end)); hold on; plot(dats(:,1,1)); hold off


%% MANUAL DECISION
% BFanalysis: "true" if cell reached depolarization block, "false" is not.
% Rationale: BFanalysis will be used as inclusion/exclusion criterion to 
% count cells as displaying persistent firing, in order to avoid countin
% cells as non-persistent firing just because they were not stimulated long
% enough. Apply independently of whether cell displays persistent firing in
% order to avoid bias in the opposite direction.
cd(pathOut)
BFanalysis = true; 
saveData.analysis.suitableForBFanalysis = BFanalysis; save(fOut, 'saveData')
close all