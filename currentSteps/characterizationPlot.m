function sweepInds = characterizationPlot(allData, cellNo)
% INPUTs
%   allData       -   allData structure from barrage master script
%   cellNo        -   # of cell to be analyzed
% OUTPUTs
%   sweepInds     -   [first sweep, subthreshold sweep, suprathreshold sweep, maxFR sweep]

    % loop vars
    notFound = 1;
    jter = 1;
    
    % find max FR sweep
    [~, maxInd] = max(cell2mat(allData(cellNo).analysis.APstats.nSpikesDuringStep));

    % find sweep with max. FR
    while notFound == 1 
        if allData(cellNo).analysis.APstats.nSpikesDuringStep{jter} > 0 &&...
           allData(cellNo).analysis.APstats.nSpikesDuringStep{jter+1} > 0 &&... 
           allData(cellNo).analysis.APstats.nSpikesDuringStep{jter+2} > 0
            
            subThreshInd = jter - 1;
            supraThreshInd = jter;
            notFound = 0;

        else
            jter = jter + 1;
        end
    end
    
    % return
    sweepInds = [1  subThreshInd supraThreshInd maxInd];
    
end