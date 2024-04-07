function [spikeTimeLine, sweepTimeLine, spikeTimes, stepTimes, interSweepTimes] = extractSpikeAndSweepTimes(stepTime, APtimes, t, interSweepTime)
% INPUTS:
%   stepTime - 2-value vector with step start and end within one sweep
%   APtimes - cell of vectors with ts of all action potentials of all
%             sweeps; APtimes{x} = [y z ...], x = sweep#, [y z ...] = AP ts
%             in sweep x           
%   t - vector of ts
%   interSweepTime - time between sweeps
% OUPUTS explained below

    % prepare data structures
    sweepLength = zeros(length(APtimes)*max(t) + (length(APtimes)-1)*interSweepTime, 1); % vector of all sweeps ts in ms resolution
    spikeTimeLine = sweepLength; % vector of all sweeps ts in ms resolution to hold spike ts
    sweepTimeLine = sweepLength; % vector of all sweeps ts in ms resolution to hold sweep component identity; 0 = baseline, 1 = current step, 2 = inter-sweep time
    spikeTimes = {};
    stepTimes = NaN(length(APtimes), 2); % vector of step start and end times
    interSweepTimes = NaN(length(APtimes), 2); % vector of inter sweep interval start and end times

    % loop sweeps
    timeCounter = 0;
    for iter = 1:length(APtimes)

        % if spikes in current sweep, add to spikeTimes
        if ~isempty(APtimes{iter})
            spkTmp = int64(APtimes{iter}); % NOTE: THIS CONVERSION REDUCES PRECISION TO MILISECOND LEVELS (which is not of concern for this plot)
            spkTmp = spkTmp + timeCounter; % indent spike times by current time
            spikeTimeLine(spkTmp) = 1; % mark spike times in time line
            spikeTimes(iter) = {APtimes{iter} + timeCounter};
        end

        % get step times in two different formats #1
        sweepTimeLine(stepTime(1):stepTime(2)) = 1; % mark steps in time line
        stepTimes(iter, :) = stepTime;

        % progress time
        timeCounter = timeCounter + max(t) + interSweepTime;
        stepTime = stepTime + + max(t) + interSweepTime;

        % get inter-sweep times in two different formats
        if ~(interSweepTime == 0)
            interSweepTimes(iter, :) = [timeCounter-interSweepTime timeCounter];
            sweepTimeLine(timeCounter-interSweepTime : timeCounter) = 2;
        end
    end % end for loop
    
end