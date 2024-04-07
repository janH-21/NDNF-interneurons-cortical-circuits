function lDir = getFilesFromDir(sDir, rxPattern, rxMode)
% INPUT:
% sDir          -   specify directory as string
% rxPattern     -   specify regex search pattern as string
% rxMode        -   spexify regex mode as string
% OUTPUT:
% lDir          -   list of filenames in directory matching pattern

% get files in directory
sDir = dir(sDir);
% perform regex search
rxDir = regexp({sDir.name}, rxPattern, rxMode);
% create string array
lDir = strings;

% loop sDir.names and extract those files matching the search
    for nfile = 1:length(sDir)
        if ~isempty(rxDir{nfile})
            % lDir(1) = sDir(nfile).name;
            lDir(end+1,1) = string(sDir(nfile).name);
        end
    end

% remove 1st element, which is empty and was created creating the array
lDir = lDir(2:end);

end