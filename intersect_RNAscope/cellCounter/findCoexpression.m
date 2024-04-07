function [tbl, tbl_intersect, fig] = findCoexpression(inputFile, clustDist, res, experimentID)
% This function takes the exported data from CellCounter in .csv format and
% calculates coexpression of markers in cells based on proximity of marks.
% The results are automatically plotted and saved to .csv files
%


    
    % check inputs
    if isempty(inputFile), [file, path] = uigetfile({'*.csv'}); inputFile = fullfile(path, file); end
    if isempty(clustDist), clustDist = input('Accepted marker distance [um]:'); end
    if isempty(res), res = input('Data resolution [um/pxl]:'); end
    if isempty(experimentID), experimentID = input('Experiment ID:','s'); end

    % load file
    tbl = readtable(inputFile);
    
    % convert pixles to um
    tbl.x_um = tbl{:,2}*res;
    tbl.y_um = tbl{:,3}*res;
    tbl.dPia_um = tbl{:,4}*res;
    
    % rename pxl columns
    tbl.Properties.VariableNames{2} = 'x_pxl';
    tbl.Properties.VariableNames{3} = 'y_pxl';
    tbl.Properties.VariableNames{4} = 'dPia_pxl';
    
    % add metadata
    tbl.cellNo = (1:size(tbl,1))';
    tbl.experimentID = repmat(experimentID,size(tbl,1),1);
    tbl.cellID = [tbl.experimentID repmat('_cell',size(tbl,1),1) cell2mat(compose('%05d',tbl.cellNo))]; 
    
    % calculate y diffs and x diffs & find cells
    sets = {};
    setCounter = 1;
    for iter = 1:size(tbl,1)
        for jter = 1:size(tbl,1)
            
            % if close enough to each other
            if abs(tbl.x_um(iter) - tbl.x_um(jter)) < clustDist & abs(tbl.y_um(iter) - tbl.y_um(jter)) < clustDist & iter ~= jter
                
                % pull out cell IDs
                cell01 = convertCharsToStrings(tbl.cellID(iter,:));
                cell02 = convertCharsToStrings(tbl.cellID(jter,:));
                
                if setCounter == 1
                    sets(setCounter) = {[cell01, cell02]}; % add first item and 
                    setCounter = setCounter + 1; % increase counter
                elseif setCounter > 1
                    % check if set overlapping with one set on list, if yes, union, if no, add
                    included = false;
                    for kter = 1:length(sets)
                        if ~isempty(intersect([cell01, cell02], sets{kter})) % if intersection (i.e. overlap) with other set
                            sets(kter) = {union([cell01, cell02], sets{kter})}; % union and
                            included = true; % count as included
                        end
                    end
                    if ~(included == true) % if not yet included
                        sets(setCounter) = {[cell01, cell02]}; % add as new item and
                        setCounter = setCounter + 1; % increase counter
                    end                  
                end
      
            end % end if close enough
        end % end jter
    end % end iter
   
    % add single marker cells
    entryCount = length(sets)+1;
    for iter = 1:size(tbl,1)
        included = false;
        cell_i = convertCharsToStrings(tbl.cellID(iter,:));
        for jter = 1:length(sets)
            if ismember(cell_i, sets{jter})
                included = true;
            end
        end
        if included == false
            sets(entryCount) = {cell_i};
            entryCount = entryCount + 1;
        end
    end
    
    % make second table with intersect IDs
    tbl_intersect = table(); % create empty table
    counter = zeros(size(sets,2),length(unique(tbl.counter))); % create zero matrix to hold counters
    x_um = nan(size(sets,2),length(unique(tbl.counter)));
    y_um = nan(size(sets,2),length(unique(tbl.counter)));
    dPia_um = nan(size(sets,2),length(unique(tbl.counter)));
    
    % fill table
    for iter = 1:size(sets,2)
        % fill metadata
        tbl_intersect.cellNo(iter) = iter;
        tbl_intersect.cellID(iter) = {[experimentID, '_intersect_cell', cell2mat(compose('%05d',iter))]}; 
        tbl_intersect.experimentID(iter) = {experimentID};
        
        % extract counters and distances
        for jter = 1:size(tbl,1)
            if ismember(convertCharsToStrings(tbl.cellID(jter,:)),sets{iter})
                
                % check counter
                counter(iter,tbl.counter(jter)) = 1;
                % 1 = EYFP
                % 2 = Npy
                % 3 = Ndnf
                
                % insert distances
                x_um(iter,tbl.counter(jter)) = tbl.x_um(jter);
                y_um(iter,tbl.counter(jter)) = tbl.y_um(jter);
                dPia_um(iter,tbl.counter(jter)) = tbl.dPia_um(jter);
                
            end
        end    
    end
    
    % --> edit this part to enable more or less counters <-- 
    tbl_intersect.counter_1 = counter(:,1);
    tbl_intersect.counter_2 = counter(:,2);
    tbl_intersect.counter_3 = counter(:,3);
    
    % average distances
    tbl_intersect.x_um = mean(x_um,2,'omitnan');
    tbl_intersect.y_um = mean(y_um,2,'omitnan');
    tbl_intersect.dPia_um = mean(dPia_um,2,'omitnan');
    
    % plot results
    fig = figure('Name', experimentID); ax = axis; hold on; 
    set(gca, 'color', [0.45 0.45 0.45]); set(gcf, 'color', [0.7 0.7 0.7]); set(gca, 'FontSize', 20)
    % plot table 
    scatter(tbl.x_um(tbl.counter == 1), tbl.y_um(tbl.counter == 1), 'o', 'filled', 'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'MarkerFaceAlpha', 0.5, 'LineWidth',3)
    scatter(tbl.x_um(tbl.counter == 2), tbl.y_um(tbl.counter == 2), 'o', 'filled', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'MarkerFaceAlpha', 0.5, 'LineWidth',3)
    scatter(tbl.x_um(tbl.counter == 3), tbl.y_um(tbl.counter == 3), 'o', 'filled', 'MarkerEdgeColor', 'yellow', 'MarkerFaceColor', 'yellow', 'MarkerFaceAlpha', 0.5, 'LineWidth',3)
    % plot intersect table
    scatter(tbl_intersect.x_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 1),...
        tbl_intersect.y_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 1), 'p', 'MarkerEdgeColor', 'black', 'LineWidth',2)
    scatter(tbl_intersect.x_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 0),...
        tbl_intersect.y_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 0), '+', 'MarkerEdgeColor', 'black', 'LineWidth',2)
    scatter(tbl_intersect.x_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 0 & tbl_intersect.counter_3 == 1),...
        tbl_intersect.y_um(tbl_intersect.counter_1 == 1 & tbl_intersect.counter_2 == 0 & tbl_intersect.counter_3 == 1), '*', 'MarkerEdgeColor', 'black', 'LineWidth',2)
    scatter(tbl_intersect.x_um(tbl_intersect.counter_1 == 0 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 1),...
        tbl_intersect.y_um(tbl_intersect.counter_1 == 0 & tbl_intersect.counter_2 == 1 & tbl_intersect.counter_3 == 1), 'x', 'MarkerEdgeColor', 'black', 'LineWidth',2)
    % cosmetics
    title(strrep(experimentID, '_', ' ')); ylabel('um'); xlabel('um'); grid on
    legend('cOn/fOn EYFP', 'npy mRNA', 'ndnf mRNA', 'EYFP + npy + ndnf', 'EYFP + npy','EYFP + ndnf','npy + ndnf',...
        'Location','eastoutside', 'Color', [0.45 0.45 0.45])
    
    % save results
    [savePath, fileName, ~] = fileparts(inputFile);
    writetable(tbl,strjoin({savePath, '\', fileName, '_extended.csv'}, ''), 'Delimiter', 'comma');
    writetable(tbl_intersect, strjoin({savePath, '\', fileName, '_intersect.csv'}, ''), 'Delimiter', 'comma');
    savefig(fig, strjoin({savePath, '\', fileName, '.fig'}, ''))

end