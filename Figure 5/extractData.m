%% format pharmacology data (run current steps cross cell master on synaptic blocker data first)
BF_onset = NaN(length(allData),1);
APsToPF = NaN(length(allData),1);
for iter = 1:length(allData)
	if ~isnan(allData(iter).analysis.nBarrageSpikes) &...
       ~isempty(find(allData(iter).analysis.nBarrageSpikes > 10, 1))
   
    	BF_onset(iter) = find(allData(iter).analysis.nBarrageSpikes > 5, 1);
        APsToPF(iter) = sum([allData(iter).analysis.APstats.nSpikesTotal{1:BF_onset(iter)}]);
        
	end
end
% format results --> saved as synaptic_steps_BFstats.mat
results = [BF_onset, APsToPF];
results = results(~isnan(BF_onset),:);


%% control data (run current steps cross cell master on control data first)
useInd = ones(length(allData),1);

for iter = 1:length(allData)
   if allData(iter).cellType ~= "NDNF+"
       useInd(iter) = 0; end
   if allData(iter).analysis.suitableForBFanalysis ~= 1
       useInd(iter) = 0; end
   if isnan(allData(iter).analysis.nBarrageSpikes)
       useInd(iter) = 0; end
   if isempty(allData(iter).analysis.nBarrageSpikes)
       useInd(iter) = 0; end
   if sum(allData(iter).analysis.nBarrageSpikes) <= 5
       useInd(iter) = 0; end
end


BF_onset = NaN(length(allData),1);
APsToPF = NaN(length(allData),1);
% manual encoding for 4 cells escaping automatic processing
for iter = 1:length(allData)
	if useInd(iter) == 1
    	BF_onset(iter) = find(allData(iter).analysis.nBarrageSpikes > 5, 1);
        APsToPF(iter) = sum([allData(iter).analysis.APstats.nSpikesTotal{1:BF_onset(iter)}]);
    end
    if iter == 47
        BF_onset(iter) = 71;
        APsToPF(iter) = sum([allData(iter).analysis.APstats.nSpikesTotal{1:BF_onset(iter)}]);
    end
	if iter == 104
        BF_onset(iter) = 39;
        APsToPF(iter) = sum([allData(iter).analysis.APstats.nSpikesTotal{1:BF_onset(iter)}]);
	end
	if iter == 112
        BF_onset(iter) = 70;
        APsToPF(iter) = sum([allData(iter).analysis.APstats.nSpikesTotal{1:BF_onset(iter)}]);
    end
end

% format data --> saved as control_steps_BFstats.mat
results = [BF_onset, APsToPF];
results = results(~isnan(BF_onset),:)


%% example recording pharmacology
figure
plot(tts, reshape(allData(5).recordingData.rawData_mV,[],1), "Color", "Black")
    xlim([1.24 1.29]*10^5)
    line([1.24 1.24]*10^5,[0 20],"Color", "Black", "LineWidth", 1) % 20mV
    line([1.24 1.24]*10^5+[0 500],[0 0],"Color", "Black", "LineWidth", 1) % 500ms
    axis off