% ---------------------------PREPARE WORKSPACE ---------------------------
clearvars -except dirPath lDir fName uter fID

% ------------------------ USER EDITABLE BLOCK ---------------------------
% directory for .abf files to read
searchPath = dirPath; 
% directory for result files to save (can be the same or different)
savePath = searchPath; 
% output file identifier
% fID = 'Ndnf_nonFS';
% LED duration in ms
pulseDuration = 0.5;
% save data boolean
saveData = 0;

% --------------------------- MAIN SCRIPT --------------------------------
% load current file
cd(searchPath)
%[fName, path] = uigetfile({'*.abf'});
%[d,si,h]=abfload(fullfile(path,fName)); 
[d,si,h]=abfload(fullfile(fName));
if size(d,2) == 4 
    d(:,2,:) = d(:,4,:);
    d = d(:,1:2,:);
end
t = (1:size(d,1))*si/1000; % time stamps in ms

% sanity check: returns 1 if timing always exactly the same
[LEDsane, LEDtimes] = IPSC_checkLEDtimes(d, si);
% LEDtimes = [2.7188 3.7188 4.7188 5.7188]*1000;
% LEDtimes = [2.6568 3.6568 4.6568 5.6568]*1000;

% NOTE: ARTEFACT REMOVAL: The aim is to get as close as possible to 15
% sweeps per cell. For most cells, this is no problem. In some cells,
% single sweeps had to be excluded if they contained electric artefacts or
% spiking activity in the analysis time window. In this cases, the sweep(s)
% immediately following sweep 15 were added into the analysis to reach a
% total of 15 sweeps.
% Example: If sweeps 2 and 10 had to be skipped due to artifacts, sweeps 1,
% 3 to 9, and 11 to 17 were used.


% ------------------------- pre-process data ------------------------------
d1 = squeeze(d(:,1,sweepNo));
d1_mean = mean(d1,2);
d1_mean = d1_mean - mean(d1_mean(2300*1000/si:2700*1000/si));


[filtDats1, filtDats2, filtDats3] = filtData(d1_mean, 100, 500);
% NOTE:
    % filtDats1 for plots and charge
    % filtDats2 for high frequency domains, as peak and rise time
    % filtDats3 for low frequency domains as decay time
   

% plot filter results
fFilt = plotFilterResults(t, d1_mean,...
                          filtDats1, filtDats2, filtDats3,...
                          LEDtimes, pulseDuration);

% ---------------------------- analyse data -------------------------------                      
% z-Score data
SD_rise = std(filtDats2(2200*1000/si:2700*1000/si));
zDats_rise = (filtDats2 - mean(filtDats2(2200*1000/si:2700*1000/si)))./SD_rise;
SD_fall = std(filtDats3(2200*1000/si:2700*1000/si));
zDats_fall = (filtDats3 - mean(filtDats3(2200*1000/si:2700*1000/si)))./SD_fall;

% extract IPSC stats
[IPSCchar, fStats] = extractIPSCstats(t, LEDtimes, filtDats2, filtDats3, zDats_rise, zDats_fall, si);
IPSCresults = array2table([convertCharsToStrings(fName) IPSCchar], 'VariableNames',...
{'fileNo',...
 'IPSCno_1', 'amplitude_1_pA', 'tmaxAmpl_1_ms', 'tZStart_1_ms', 'tZEnd_1_ms', 'dtZ_1_ms', 't20rise_1_ms', 't80rise_1_ms', 't80decay_1_ms', 't20decay_1_ms', 'dtRise_1_ms', 'dtDecay_1_ms', 'charge_1_z3_pC', 'charge_1_1000ms_pC',...
 'IPSCno_2', 'amplitude_2_pA', 'tmaxAmpl_2_ms', 'tZStart_2_ms', 'tZEnd_2_ms', 'dtZ_2_ms', 't20rise_2_ms', 't80rise_2_ms', 't80decay_2_ms', 't20decay_2_ms', 'dtRise_2_ms', 'dtDecay_2_ms', 'charge_2_z3_pC', 'charge_2_1000ms_pC',...
 'IPSCno_3', 'amplitude_3_pA', 'tmaxAmpl_3_ms', 'tZStart_3_ms', 'tZEnd_3_ms', 'dtZ_3_ms', 't20rise_3_ms', 't80rise_3_ms', 't80decay_3_ms', 't20decay_3_ms', 'dtRise_3_ms', 'dtDecay_3_ms', 'charge_3_z3_pC', 'charge_3_1000ms_pC',...
 'IPSCno_4', 'amplitude_4_pA', 'tmaxAmpl_4_ms', 'tZStart_4_ms', 'tZEnd_4_ms', 'dtZ_4_ms', 't20rise_4_ms', 't80rise_4_ms', 't80decay_4_ms', 't20decay_4_ms', 'dtRise_4_ms', 'dtDecay_4_ms', 'charge_4_z3_pC', 'charge_4_1000ms_pC'}) %#ok<NOPTS>

pause
saveFiles(savePath, fID, fName, IPSCresults, d1_mean, fFilt, fStats);
close all;


%% --------------------------- SUBROUTINES --------------------------------
function [filtDats1, filtDats2, filtDats3] = filtData(dats, movAvgWin1, movAvgWin2)
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

function fFilt = plotFilterResults(t, d1_mean, filtDats1, filtDats2, filtDats3, LEDtimes, pulseDuration)
    fFilt = figure('units','normalized','outerposition',[0 0 1 1]);
    plot(NaN); hold on
    for iter = 1:length(LEDtimes)
        rectangle('Position', [LEDtimes(iter) -200 pulseDuration 500],...
            'FaceColor', 'cyan', 'EdgeColor', 'cyan');
    end
    plot(t, d1_mean, "Color", [.85 .85 .85])
    plot(t, filtDats1, "Color", "black")
    plot(t, filtDats2, "Color", "blue")
    plot(t, filtDats3, "Color", "red")
    plotAesthetics(gca, 1, 14);
    ylabel("Amplitude [pA]"); xlabel("Time [ms]");
    ylim([-10 150]); xlim([2500 7000]);
    legend('LED',...
           'raw',...
           '8-pole Butterworth Fc = 1000Hz',...
           'Butterworth + MovAvg, WinSz = 100',...
           'Butterworth + MovAvg, WinSz = 500')
end

function [IPSCchar, fStats] = extractIPSCstats(t, LEDtimes, filtDats1, filtDats2, zDats_rise, zDats_fall, si)  
    IPSCchar = NaN(1,length(LEDtimes)*14);
    fStats = figure('units','normalized','outerposition',[0 0 1 1]);
    for iter = 1:length(LEDtimes)
        base_ms = 50;
        % ---------- ANALYSIS ----------
        ts_start = (LEDtimes(iter)-(base_ms+1))*1000/si;
        ts_end = (LEDtimes(iter)+1000)*1000/si;
        % extract IPSC max. amplitude & time stamp
        dats_pA = filtDats2(ts_start:ts_end);
        dats_pA = dats_pA - mean(dats_pA(1:base_ms*1000/si));
        [maxAmpl_pA, maxAmpl_Ind] = max(dats_pA(1:(base_ms+100)*1000/si));
        % charge
        t_charge = t(ts_start:ts_end);
        dats_pA_charge  = filtDats1(ts_start:ts_end);
        dats_pA_charge = dats_pA_charge - mean(dats_pA_charge(1:base_ms*1000/si));
        charge_pC_1000ms = trapz(t_charge((base_ms+1)*1000/si:(base_ms+1001)*1000/si)./1000,...
                                 dats_pA_charge((base_ms+1)*1000/si:(base_ms+1001)*1000/si)); % charge for complete trace between LEDtimes
        IPSCchar(1+14*(iter-1)) = iter;
        IPSCchar(2+14*(iter-1)) = maxAmpl_pA;
        IPSCchar(3+14*(iter-1)) = (maxAmpl_Ind)*si/1000 - base_ms - 1;
        IPSCchar(14+14*(iter-1)) = charge_pC_1000ms;
        try % try-catch to avoid crushing the script is no connectivity
           % extract z-scored data for current IPSC
            dats_z_rise = zDats_rise(ts_start:ts_end);
            dats_z_fall = zDats_fall(ts_start:ts_end);
            z3_start = find(dats_z_rise((base_ms+1)*1000/si:end) > 3, 1);
            z3_start = z3_start + (base_ms+1)*1000/si;
            z3_end = find(dats_z_fall(maxAmpl_Ind:end) < 3, 1) - 1 + maxAmpl_Ind;
            if isempty(z3_end); z3_end = (base_ms+1001)*1000/si; end
           % extract 20->80% rise times
           	dats_pA_rise = dats_pA(1:maxAmpl_Ind);
            dats_pA_rise = dats_pA_rise ./ max(dats_pA_rise) * 100; 
            rise20 = find(dats_pA_rise > 20, 1);
            rise80 = find(dats_pA_rise > 80, 1);  
           % extract  80->20% decay times 
            dats_pA_fall = dats_pA(maxAmpl_Ind:end);
            dats_pA_fall = dats_pA_fall ./ max(dats_pA_fall) * 100;
            decay20 = find(dats_pA_fall < 20, 1);
            decay80 = find(dats_pA_fall < 80, 1);
           % calculate charge
            charge_pC_z3 = trapz(t_charge(z3_start:z3_end)./1000,...
                                 dats_pA_charge(z3_start:z3_end)); % charge for area where z > 3; ./1000 to convert to seconds
           % ------------ PLOTs ------------
           % plot pA data, max. ampl & charge
            subplot(4,9,[1:3]+(iter-1)*9)  
                plot(t_charge, dats_pA_charge, 'Color', [0.5 0.5 0.9]); hold on; grid on;
                area(t_charge((base_ms+1)*1000/si:(base_ms+1001)*1000/si), dats_pA((base_ms+1)*1000/si:(base_ms+1001)*1000/si)',...
                     'EdgeColor', 'none', 'FaceColor', 'blue', 'FaceAlpha', 0.05)
                area(t_charge(z3_start:z3_end), dats_pA(z3_start:z3_end)',...
                     'EdgeColor', 'none', 'FaceColor', 'blue', 'FaceAlpha', 0.1)
                plot(t_charge, dats_pA, 'Color', 'black');
                line(xlim, [maxAmpl_pA maxAmpl_pA],...
                     "LineStyle", ":", "Color", "black", "LineWidth", 1.5)
                line([maxAmpl_Ind maxAmpl_Ind]./1000*si+ts_start/1000*si, ylim,...
                     "LineStyle", ":", "Color", "black", "LineWidth", 1.5)
                legend('Butterworth', 'Butterworth + MovAvg100',...
                       'Charge by > 3 y-Score', 'Charge over 1000ms',...
                       'Peak Amplitude', 'Location', 'best')
                xlabel('Time [ms]'); ylabel('Amplitude [pA]')
                title('IPSC')
            % plot z-scored data
             subplot(4,9,[4:6]+(iter-1)*9)
                plot(dats_z_rise, 'LineWidth', 1); hold on; grid on;
                plot(dats_z_fall, 'LineWidth', 1);
                line(xlim, [3 3], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1); 
                xlabel('Samples'); ylabel('z-Score'); title('z-Scored IPSC')
                line([maxAmpl_Ind maxAmpl_Ind], ylim,'Color', 'red', 'LineStyle', '-.', 'LineWidth', 1)
                line([z3_start z3_start], ylim,'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5)
                line([z3_end z3_end], ylim,'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5); 
                legend('rise data', 'fall data', 'threshold: z-Score = 3 ',...
                       'peak amplitude', 'z-Score threshold crossing', '', 'Location', 'best')
            % plot rising component
            subplot(4,9,7+(iter-1)*9)
                plot(dats_pA_rise, 'LineWidth', 1); grid on; ylim([-10 110])
                xlabel('Samples'); ylabel('Amplitude [%]'); title('IPSC rise')
                line([xlim], [20 20], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1); 
                line([xlim], [80 80], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1); 
                line([rise20 rise20],[ylim],'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5)
                line([rise80 rise80],[ylim],'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5)
            % plot decaying component
            subplot(4,9,[8:9]+(iter-1)*9)
                plot(dats_pA_fall, 'LineWidth', 1); grid on; ylim([-10 110])
                xlabel('Samples'); ylabel('Amplitude [%]'); title('IPSC decay')
                line([xlim], [20 20], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1); 
                line([xlim], [80 80], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1); 
                line([decay20 decay20],[ylim],'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5)
                line([decay80 decay80],[ylim],'Color', 'red', 'LineStyle', ':', 'LineWidth', 1.5)
                legend('data', '%-threshold: amplitude 20/80', '', '%-threshold crossing', '', 'Location', 'best')
            % -------- CoLLECT DATA --------
            IPSCchar(4+14*(iter-1)) = (ts_start + z3_start - 1)/1000*si - LEDtimes(iter);
            IPSCchar(5+14*(iter-1)) = (ts_start + z3_end - 1)/1000*si - LEDtimes(iter);
            IPSCchar(6+14*(iter-1)) = abs(z3_end-z3_start)./1000*si;
            IPSCchar(7+14*(iter-1)) = (ts_start + rise20 - 1)/1000*si - LEDtimes(iter);
            IPSCchar(8+14*(iter-1)) = (ts_start + rise80 - 1)/1000*si - LEDtimes(iter);
            IPSCchar(9+14*(iter-1)) = (ts_start + maxAmpl_Ind + decay20 - 2)/1000*si - LEDtimes(iter);
            IPSCchar(10+14*(iter-1)) = (ts_start + maxAmpl_Ind + decay80 - 2)/1000*si - LEDtimes(iter);
            IPSCchar(11+14*(iter-1)) = abs(rise20-rise80)/1000*si;
            IPSCchar(12+14*(iter-1)) = abs(decay20-decay80)/1000*si;
            IPSCchar(13+14*(iter-1)) = charge_pC_z3;
        catch
            IPSCchar(1+14*(iter-1)) = iter;
            strjoin({'Unable to extract stats for IPSC #', num2str(iter)},'')
        end
    end
end

function saveFiles(savePath, fID, fName, IPSCresults, d1_mean, fFilt, fStats)
    % save IPSC statistics
    fName_IPSCs = convertCharsToStrings(strjoin({fID, '_results.csv'},''));
    fName_IPSCs = fullfile(savePath, fName_IPSCs);
    if isfile(fName_IPSCs)
        saveFile = readtable(fName_IPSCs);
        saveFile = [saveFile; IPSCresults] %#ok<NOPRT>
        writetable(saveFile, fName_IPSCs);
    else
        writetable(IPSCresults, fName_IPSCs);
    end
    % get file name
    fName_split = strsplit(fName, '.');
    fName_trunk = convertStringsToChars(fName_split{1});
    fName_trunk = fullfile(savePath, fName_trunk);
    % save mean trace & figures
    save(strjoin({fName_trunk, '_meanTrace.mat'}, ''), 'd1_mean')
    savefig(fFilt, strjoin({fName_trunk, '_analysis_filt.fig'}, ''))
    savefig(fStats, strjoin({fName_trunk, '_analysis_stats.fig'}, ''))
end

function inpt = repeatInput(cellOfChars, prompt)
    repeat = true;
    while repeat
        inpt = input(prompt,"s");
        if find(strcmp(inpt, cellOfChars))
            repeat = false;
        end
    end
end