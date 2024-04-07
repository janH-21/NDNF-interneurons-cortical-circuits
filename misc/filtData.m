function [filtDats1, filtDats2, filtDats3] = filtData(dats, movAvgWin1, movAvgWin2)
% filtDats1 for plots and charge
% filtDats2 for high frequency domains, as peak and rise time
% filtDats3 for low frequency domains as decay time
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
    filtDats1 = filtfilt(SOS, G, dats);
    filtDats2 = movmean(filtDats1, movAvgWin1);
    filtDats3 = movmean(filtDats1, movAvgWin2);
end