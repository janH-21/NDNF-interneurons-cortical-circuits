% The program imports different channels as separate single channel .tif 
% files, the paths to which are provided by the user in the first section 
% below. The channel histograms need to be adjusted prior to loading the 
% files (e.g. using ZEN or Fiji/ImageJ), and cannot be adjusted further in 
% this program. The channel .tif which is used to identify the pia should 
% be loaded in fName_Ch1.
%
% Workflow:
% 1. Create single channel, histogram-adjusted 8-bit .tif files 
%   (ZEN, Fiji,...)
% 2. Provide the the file paths in the user editable section below.
% 3. Channel 1 will be displayed, and you can demarcate the pia with
%    left-mouse clicks
% 4. When finished, press Enter
% 5. Now, cells can be marked using the keys 1-4, channels can be changed
%    using keys q-r.
% 6. When finished, press space
% 7. Results will be saved to the current working directory.
%
% NOTE that all measurements and calculations use pixles as units.
%
% Writen in Matlab R2020b. Jan Hartung, March 2021
%% -------------------------- user editable block -------------------------
% provide chanel file names
fName_Ch1 = "Ndnf168-9-LeftSection-LeftHemi-20x-Orthogonal Projection-Adjust-splitchannels_DAPI(lightblue).tif";
fName_Ch2 = "Ndnf168-9-LeftSection-LeftHemi-20x-Orthogonal Projection-Adjust-splitchannels_B2(red).tif";
fName_Ch3 = "Ndnf168-9-LeftSection-LeftHemi-20x-Orthogonal Projection-Adjust-splitchannels_Ndnf(green).tif";
fName_Ch4 = "Ndnf168-9-LeftSection-LeftHemi-20x-Orthogonal Projection-Adjust-splitchannels_VIP(blue).tif";

% file prefix for output file
fOut_pre = "originalIntersect_slice03_20x_allChannels_orthogonalProjection(z11-19)_stich";

% prefered colors for the channel markers can be chosen, purely aesthetic
markerColor{1} = "red";
markerColor{2} = "green";
markerColor{3} = "blue";
markerColor{4} = "yellow";
% -------------------------- end user editable block ----------------------


%% read files
fTif_Ch1 = Tiff(fName_Ch1,'r'); vTif_Ch{1} = read(fTif_Ch1);
fTif_Ch2 = Tiff(fName_Ch2,'r'); vTif_Ch{2} = read(fTif_Ch2);
fTif_Ch3 = Tiff(fName_Ch3,'r'); vTif_Ch{3} = read(fTif_Ch3);
fTif_Ch4 = Tiff(fName_Ch4,'r'); vTif_Ch{4} = read(fTif_Ch4);


%% create figure
fig = figure;
imshow(vTif_Ch{1}); hold on
set(fig,'units','normalized','outerposition',[0 0 1 1],'color','black');


%% draw pia
pia=[];
[pia(1,:), pia(2,:)] = ginput;
int_pia = min(pia(1,:)):0.1:max(pia(1,:));
spline_pia = spline( pia(1,:) , pia(2,:) , int_pia);
pia = [int_pia; spline_pia];
plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white"); hold on


%% main counter
% matrix to collect click history
clickHistory = []; 

% counting loop
repeat = true;
counter = 1;
while repeat
    
    [x,y,button] = ginput(1);
    switch button
        
        % add marker commands
        case 49  % 1
            [fig, clickHistory, counter] = addCounter(fig, clickHistory, x, y, 1, counter, markerColor);
        case 50  % 2
            [fig, clickHistory, counter] = addCounter(fig, clickHistory, x, y, 2, counter, markerColor);
        case 51  % 3
            [fig, clickHistory, counter] = addCounter(fig, clickHistory, x, y, 3, counter, markerColor);
        case 52  % 4
            [fig, clickHistory, counter] = addCounter(fig, clickHistory, x, y, 4, counter, markerColor);
            
        % switch channel commands
        case 113 % q
            fig = plotImage(fig, vTif_Ch{1}, clickHistory, markerColor); 
            plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white");
        case 119 % w
            fig = plotImage(fig, vTif_Ch{2}, clickHistory, markerColor);
            plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white");
        case 101 % e
            fig = plotImage(fig, vTif_Ch{3}, clickHistory, markerColor); 
            plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white");
        case 114 % r
            fig = plotImage(fig, vTif_Ch{4}, clickHistory, markerColor); 
            plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white");
            
        % end selection + void commands   
        case 32  % space     
            repeat = false;
        otherwise
            
    end
end


%% calculate results
distPia = shortestDistance(clickHistory, pia);
CellCounterResults = [clickHistory distPia];


%% save to file
CellCounterResults_export = array2table(CellCounterResults, 'VariableNames', {'x_coord_pxls', 'y_coord_pxls', 'marker_ID', 'distance_to_pia_pxls'});
fOut = join([fOut_pre, 'cellCounting_results.csv'], '_');
writetable(CellCounterResults_export, fOut)  

%% functions
function [fig, clickHistory, counter] = addCounter(fig, clickHistory, x, y, markerN, counter, markerColor)
        clickHistory(counter,:) = [x y markerN];
        counter = counter + 1;
        text(x, y, int2str(markerN), "Color", markerColor{markerN}, 'FontWeight', 'bold')
end

function fig = plotImage(fig, vTif_Ch, clickHistory, markerColor)  
    % plot image
    imshow(vTif_Ch); hold on 
    set(fig,'units','normalized','outerposition',[0 0 1 1],'color','black');   
    % replot marker
    for iter = 1:size(clickHistory,1)
        clickHistory(iter,1)
        text(clickHistory(iter,1), clickHistory(iter,2), int2str(clickHistory(iter,3)), "Color", markerColor{clickHistory(iter,3)}, 'FontWeight', 'bold')
    end  
end

function distPia = shortestDistance(clickHistory, pia)
% Calculates euclidean distance to nearset point in pia. Due to the dense
% spline interpolation, this can be considered equal to the distance to the
% continuous pia for all practical purposes
    distPia = nan(size(clickHistory,1),1);
    for iter = 1:length(distPia)
        allDist = abs(sqrt( (clickHistory(iter,1)-pia(1,:)).^2 + (clickHistory(iter,2)-pia(2,:)).^2 ));
        distPia(iter) = min(allDist);
    end
end