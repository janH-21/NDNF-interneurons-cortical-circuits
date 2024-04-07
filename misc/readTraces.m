function [traces, ts] = readTraces(dirPath, regexPattern, nSamples, si)
% INPUTS:
% dirPath          path to directory to search
% regexPattern     regex pattern for file names to include
% nSamples         recording length in samples
% si               sampling interval
% tsLims           start and end bin to extract from trace

    % list of files in directory
    lDir = getFilesFromDir(dirPath, regexPattern, 'names');
    % preallocate trace array
    traces = NaN(length(lDir), nSamples);
    % read files
    for iter = 1:length(lDir)
        load(lDir(iter));
        if exist('d1_mean') == 1
            traces(iter,:) = d1_mean;
        elseif exist('filtDats2') == 1
            traces(iter,:) = filtDats2;
        end
    end
    % calculate time stamps
    ts = (1:length(d1_mean)) * si / 1000;
end