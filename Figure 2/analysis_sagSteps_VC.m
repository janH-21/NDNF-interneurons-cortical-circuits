%% ------------------------------ analysis -------------------------------%
close all; clc; clear all;
% choose and load file 
cd("./main/data/sag_VC/")
[file, path] = uigetfile({'*.abf'});  
[d,si,h]=abfload(fullfile(path,file)); 
cd("../../analysis/Figure 2/")

% load sag data
load("sagData.mat")
nCells = length(sagData);
%sagData = struct();
%nCells = 0;

% sweep start and end
startInd = 57423; % startInd = locs(1);
endInd = 207422;  % endInd = locs(2);
startT = startInd / 1000 * 20;
endT = endInd / 1000 * 20;
% 8-pole butterworth with Fc = 1000Hz
G = [0.0038;
     0.0037;
     0.0036;
     0.0035;
     1.0000];
SOS = [1.0000    2.0000    1.0000    1.0000   -1.9369    0.9523;
       1.0000    2.0000    1.0000    1.0000   -1.8551    0.8698;
       1.0000    2.0000    1.0000    1.0000   -1.7970    0.8112;
       1.0000    2.0000    1.0000    1.0000   -1.7670    0.7811];

% voltage data
vBase = -50;
vStep = -55:-5:-105;

% ------------------------------- analysis -------------------------------%
% prepare data
ts = (1:size(d,1))*si/1000;
dats = squeeze(d(:,1,:));
maxDats = NaN(size(dats));
minDats = NaN(size(dats));
figure;
subplot(3,3,[1 8]); hold on
% iterate sweeps
for iter = 1:size(dats,2)
    % normalize sweep baseline
    dats(:,iter) = dats(:,iter) - mean(dats(1:startInd-500,iter));
    plot(ts, dats(:,iter), "Color", "black")
    % calculate sweep baseline
    sweepBase = mean(dats(1:startInd-500,iter));
    % filter and smooth with narrow window for faster dynamics
    maxDats(:,iter) = smooth(filtfilt(SOS, G, dats(:,iter)),111);
    plot(ts, maxDats, "LineWidth", 2, "Color", "blue")
    % search for maximum within 200ms of step start
    [~, maxInd] = max(maxDats(startInd+150-1:startInd+150-1+50*1000/si, iter));
    maxInd = maxInd + startInd + 150 - 2;
    % smooth with borad window for slower dynamics
    minDats = maxDats;
    minDats(startInd+150-1:end,iter)= smooth(minDats(startInd+150-1:end,iter),3001);
    plot(ts, minDats(:,iter), "LineWidth", 2, "Color", "cyan")
    % search for minimum within 400ms after maximum
    [~, minInd] = min(minDats(maxInd:maxInd+400*1000/si, iter));
    minInd = maxInd+minInd-1;
    % calculate step baseline
    stepBase = mean(dats(1.8074e5:2.074e5,iter));
    line([xlim],[stepBase stepBase],"LineWidth",2,"Color","white")
    % calculate max/min values
    maxVal = mean(dats(maxInd-25:maxInd+25, iter));
    maxT = maxInd*20/1000;
    plot(maxT, maxVal, "r*", "LineWidth", 4)
    minVal = mean(minDats(minInd-100:minInd+100, iter));
    minT = minInd*20/1000;
    plot(minT, minVal, "go", "LineWidth", 4)
    % sag
	sagAmpl = maxVal - stepBase;
	sagFraction = (sweepBase-maxVal) / (sweepBase - stepBase);
    % collect Data
	sagData(nCells+1).Vstep_mV(iter) = vStep(iter);
    sagData(nCells+1).sweepBase_pA(iter) = sweepBase;
    sagData(nCells+1).stepBase_pA(iter) = stepBase;
    sagData(nCells+1).maxInd_samples(iter) = maxInd;
    sagData(nCells+1).maxT_ms(iter) = maxT;
    sagData(nCells+1).maxVal_pA(iter) = maxVal;
    sagData(nCells+1).minInd_samples(iter) = minInd;
    sagData(nCells+1).minT_ms(iter) = minT;
    sagData(nCells+1).minVal_pA(iter) = minVal;
    sagData(nCells+1).sagAmpl_pA(iter) = sagAmpl;
    sagData(nCells+1).sagFraction(iter) = sagFraction;
end

% plot aesthetics
xlim([1000 4500])
ylim([minVal*1.2 50])
plotAesthetics(gca, 2, 12);
ylabel("I_{hold}[pA]"); xlabel("Time [ms]")
title("Sag Steps Protocol - Analysis Control Plot")
subtitle(file, 'Interpreter', 'none');

% plot results
subplot(3,3,3)
    plot(sagData(nCells+1).Vstep_mV, sagData(nCells+1).sagFraction, "o-", "LineWidth", 1.5)
    plotAesthetics(gca,2,12); set(gca, 'XDir','reverse');
    title("Sag Fraction"); xlabel("V_{step} [mV]"); ylabel("Fraction")
subplot(3,3,6)
    plot(sagData(nCells+1).Vstep_mV, sagData(nCells+1).sagAmpl_pA, "o-", "LineWidth", 1.5)
    plotAesthetics(gca,2,12); set(gca, 'XDir','reverse');
    title("Sag Amplitude"); xlabel("V_{step} [mV]"); ylabel("Sag Amplitude [pA]")
subplot(3,3,9)
    plot(sagData(nCells+1).Vstep_mV, sagData(nCells+1).minVal_pA - sagData(nCells+1).stepBase_pA, "o-", "LineWidth", 1.5)
    plotAesthetics(gca,2,12); set(gca, 'XDir','reverse');
    title("Dip Amplitude"); xlabel("V_{step} [mV]"); ylabel("Dip Amplitude [pA]")
    
    
%% save data
sagData(nCells+1).fileName = file;
%save("sagData.mat", "sagData")
