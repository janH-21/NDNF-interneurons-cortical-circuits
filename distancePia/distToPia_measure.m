% Writen in Matlab R2020b. Jan Hartung, November 2022. Based on cellCounter.m from March 2021.
% NOTE that all measurements and calculations use pixles as units.

%% main loop
clear all; close all; clc
mainPath = "";
path = mainPath;
recNo = {};
imgNo = {}; 
loop = true;
counter = 1;
while loop
    % get file
    cd(path) % for convenience
    [fName, path] = uigetfile({'*.tif'}); % choose file 
    fFull = fullfile(path, fName);
    [X,map] = imread(fFull); % read file
    
    % repeat measurement loop
    measureLoop = true;
    while measureLoop
        % figure
        close all
        fig = figure;
        imshow(rot90(X), map); hold on
        set(fig,'units','normalized','outerposition',[0 0 1 1],'color','black');

        % draw pia
        pia=[];
        [pia(1,:), pia(2,:)] = ginput;
        int_pia = min(pia(1,:)):0.1:max(pia(1,:));
        spline_pia = spline( pia(1,:) , pia(2,:) , int_pia);
        pia = [int_pia; spline_pia];
        plot(pia(1,:),pia(2,:),"LineWidth",3,"Color","white"); hold on

        % mark cell
        [x,y,button] = ginput(1);
        text(x, y, "X", "Color", "white", 'FontWeight', 'bold')

        % store data
        pathSplit = strsplit(path, "/");
        recNo{counter} = pathSplit{end-1};
        imgNo{counter} = fName;
        distPia(counter) = shortestDistance([x, y], pia);

        % continue?
        inputLoop = true;
        while inputLoop
            prompt = "Proceed (p) / Repeat (r) / Finish (q)"; 
            userIn = input(prompt, "s");
            if userIn == "p"        % exit input decision and measurement loops, procced to next image
                inputLoop = false; 
                measureLoop = false; 
                counter = counter + 1; 
            elseif userIn == "r"    % exit input decision loop but repeat measurement
                inputLoop = false;
            elseif userIn == "q"    % exit all loops and proceed to analysis
                inputLoop = false;      % exit input decision tree
                loop = false;           % exit main loop
                measureLoop = false;    % exit current image measure loop
            end
        end % input loop
    end % measure loop
end % main loop

%%
resultTable = table(recNo', imgNo', distPia', 'VariableNames', {'recording', 'imgNo',  'distance_to_pia_pxls'})

%% save to file
writetable(resultTable, fullfile(mainPath, "NDNFdistToPia.csv"), "Delimiter", ",")  

%% functions
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