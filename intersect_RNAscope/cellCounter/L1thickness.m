%%========================== user editable block ==========================
% input file -> export DAPI channel from .csv to .tif
fIn = 'fileName.tif';

% meta - adapt according to used image
pxls_to_um = 0.42;
rotate_by = 300; % angle in degree
bin_size_pxls = 100;

% current CAVEATS:
% - pia detection needs to be manually double checked
% - rotation needs to be manually double checked
% - cropping needs to be manually double checked
% - L1 border detection needs to be manually double checked

%========================== end user editable block =======================
% make dir
% get file names
fparts = strsplit(fIn, '_');
dirName = strjoin(fparts(1:2), '\_');
% make dir
mkdir(dirName)

% read .tof
fDapi = Tiff(fIn,'r');
iDapi  = read(fDapi);
% select DAPI channel
dapi = iDapi(:,:,3);


%% if necessary, rotate image: imshow(imrotate(dapi, rotate_by))
dapi = imrotate(dapi, rotate_by);

%% if necessary, crop image
%dapi = dapi(0:end,:);

%% display & save 
fTif = figure('units','normalized','outerposition',[0 0 1 1]);
imshow(dapi)
saveas(fTif, fullfile(dirName, strjoin({dirName, 'DAPI.png'},'_')))


% bin figure to 100pxl bins, align and plot average intensity
fAlign = figure('units','normalized','outerposition',[0 0 1 1]);
% empty matrix to hold aligned values
dapiAligned = zeros(floor(size(dapi,1)/100),size(dapi,2));
% index array for subplots
index = reshape(1:floor(size(dapi,1)/100)*2, 2, floor(size(dapi,1)/100)).';

% iterate bins
nBins = floor(size(dapi,1)/bin_size_pxls);
for iter = 1:nBins
    
    % binned image
    subplot(nBins,2,index(iter,1))
    dapiBin = mean(dapi(1+((iter-1)*100):iter*100,:),1);
    plot(dapiBin); hold on; set(gca,'visible','off');
    % find and plot pia
    med = mean(dapiBin);
    [row, col] = find(dapiBin > med, 1);
    plot(col, dapiBin(col), 'Color', 'green', 'Marker','x')
    
    % aligned bins
    subplot(nBins,2,index(iter,2))
    dapiAligned(iter,1:size(dapiBin,2)-col+1) = dapiBin(col:end);
    plot(dapiAligned(iter,:)); set(gca,'visible','off');
    
end

saveas(fAlign, fullfile(dirName, strjoin({dirName, 'binAlign.png'},'_')))


% calculate L1 thickness
fResult = figure('units','normalized','outerposition',[0 0 1 1]);

% mean DAPI distribution over aligned bins
meanDapi = mean(dapiAligned,1);

% savgol filter to smooth for easy min/max detection
filtDapi = lowpass(meanDapi,0.001);
filtDapi = sgolayfilt(filtDapi,1,63);
filtDapi = sgolayfilt(filtDapi,1,63);
filtDapi = sgolayfilt(filtDapi,1,63);


% plot traces
subplot(2, 1, 1)
plot(meanDapi, 'Color', 'black'); hold on
plot(filtDapi, 'Color', 'blue', 'LineWidth',2);

% % find min/max
% [pksP,locsP] = findpeaks(filtDapi);
% [pksN,locsN] = findpeaks(-filtDapi);
% 
% % calculate average of L1 = bordermin/max, plot min, max and average
% L1_border_pxls = (locsP(1) + locsN(1))/2; 
% line([locsP(1) locsP(1)],[ylim], 'Color', 'red', 'LineStyle', '-.')
% line([locsN(1) locsN(1)],[ylim], 'Color', 'red', 'LineStyle', '-.')
% line([L1_border_pxls L1_border_pxls], [ylim], 'Color', 'red')

% calculate L1 border via derivative and plot
[pks,locs] = findpeaks(diff(filtDapi), 'MinPeakDistance', 100/pxls_to_um*2/3);
L1_border_pxls = locs(2);
line([L1_border_pxls L1_border_pxls], [ylim], 'Color', 'red')

% plot derivative
subplot(2, 1, 2)
plot(diff(filtDapi)); hold on
line([L1_border_pxls L1_border_pxls], [ylim], 'Color', 'red')

% calculate L1 border in um and save to .txt
L1_border_um = L1_border_pxls * pxls_to_um;
writematrix(L1_border_um, strjoin({dirName, 'result.txt'},'_'))
saveas(fResult, fullfile(dirName, strjoin({dirName, 'result.png'},'_')))

