function [PSD_data, PSD_slopes] = UVP5_PSD(dataFolder, projectNames, stationUse, depthRange, ESDmin, ESDmax, numBins, saveDir)
%% UVP5_PSD - Calculate particle size distributions (PSD) from UVP data
% Inputs:
%   dataFolder   - folder containing .mat particle data
%   projectNames - cell array of project names to include
%   stationUse   - station numbers to include
%   depthRange   - [startDepth, endDepth] in meters
%   ESDmin/ESDmax- min/max particle size in mm
%   numBins      - number of logarithmic bins
%   saveDir      - folder to save output tables
%
% Outputs:
%   PSD_data     - table of PSD counts per liter per mm for each cast
%   PSD_slopes   - table of PSD slopes for each cast

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

PSD_data   = table();
PSD_slopes = table();

% Define bin edges and centers
binEdges   = logspace(log10(ESDmin), log10(ESDmax), numBins+1);
binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
binCenters = binCenters(:);
binWidths  = diff(binEdges);

% List all .mat files for the given projects
allFiles = dir(fullfile(dataFolder,'*.mat'));
if isempty(allFiles)
    error('No .mat files in %s', dataFolder);
end

fileList = {};
for p = 1:numel(projectNames)
    matches = contains({allFiles.name}, projectNames{p}, 'IgnoreCase', true);
    fileList = [fileList; {allFiles(matches).name}']; %#ok<AGROW>
end
if isempty(fileList)
    error('No files found for projects %s', strjoin(projectNames, ', '));
end

% --- Loop through files ---
for f = 1:numel(fileList)
    filePath = fullfile(dataFolder, fileList{f});
    S = load(filePath);
    fn = fieldnames(S);
    T = S.(fn{1});

    if ~istable(T)
        warning('File %s has no table, skipping.', fileList{f});
        continue;
    end

    % Lowercase column names
    T.Properties.VariableNames = lower(T.Properties.VariableNames);
    varNamesLower = T.Properties.VariableNames;
    findVar = @(keyword) varNamesLower(contains(varNamesLower, keyword,'IgnoreCase',true));

    % Detect required columns
    projectVar  = findVar('project');
    cruiseVar   = findVar('cruise');
    stationVar  = findVar('station');
    castVar     = findVar('cast');
    literVar    = findVar('liter');
    depthVar    = findVar('depth');
    esdVar      = findVar('esd');
    datetimeVar = findVar('datetime');

    if any(cellfun(@isempty, {projectVar, cruiseVar, stationVar, castVar, literVar, depthVar, esdVar, datetimeVar}))
        warning('Missing required columns in %s, skipping.', fileList{f});
        continue;
    end

    % Filter by project
    isProj = ismember(lower(string(T.(projectVar{1}))), lower(projectNames));

    % Filter by station (allow all if stationUse is empty)
    if isempty(stationUse)
        isStn = true(height(T),1);
    else
        isStn = ismember(T.(stationVar{1}), stationUse);
    end

    % Filter by depth and ESD
    isDepth = T.(depthVar{1}) >= depthRange(1) & T.(depthVar{1}) <= depthRange(2);
    isESD   = T.(esdVar{1}) >= ESDmin & T.(esdVar{1}) <= ESDmax;

    % Apply all filters
    T = T(isProj & isStn & isDepth & isESD, :);
    if isempty(T)
        continue;
    end

    % Convert station and cast to numeric if needed
    if iscell(T.(stationVar{1})) || isstring(T.(stationVar{1}))
        T.(stationVar{1}) = double(string(T.(stationVar{1})));
    end
    if iscell(T.(castVar{1})) || isstring(T.(castVar{1}))
        T.(castVar{1}) = double(string(T.(castVar{1})));
    end

    % --- Loop over all unique (cruise, station, cast) combinations ---
    uniqueComb = unique([T.(cruiseVar{1}), T.(stationVar{1}), T.(castVar{1})], 'rows');

    for idx = 1:size(uniqueComb,1)
        cr = uniqueComb(idx,1); % cruise
        st = uniqueComb(idx,2); % station
        ca = uniqueComb(idx,3); % cast

        Tc = T(T.(cruiseVar{1})==cr & T.(stationVar{1})==st & T.(castVar{1})==ca, :);
        if isempty(Tc), continue; end

        % Histogram counts
        counts = histcounts(Tc.(esdVar{1}), binEdges);
        counts = counts(:);

        % Total volume
        totalVolumeL = numel(unique(Tc.(literVar{1}))) * 1.22;

        % Normalize counts
        countsPerL = counts ./ binWidths(:) / totalVolumeL;

        % Extract metadata
        cruiseName  = string(Tc.(cruiseVar{1})(1));
        projectName = string(Tc.(projectVar{1})(1));
        castVal     = ca;
        stationVal  = st;
        timeStart   = Tc.(datetimeVar{1})(1);
        timeEnd     = Tc.(datetimeVar{1})(end);

        % Build PSD table
        nRows = numBins;
        PSDtable = table();
        PSDtable.bin_center_mm          = binCenters;
        PSDtable.bin_width_mm           = binWidths(:);
        PSDtable.counts_per_liter_per_mm = countsPerL;
        PSDtable.project  = repmat(projectName, nRows, 1);
        PSDtable.cruise   = repmat(cruiseName, nRows, 1);
        PSDtable.cast     = repmat(castVal, nRows, 1);
        PSDtable.station  = repmat(stationVal, nRows, 1);
        PSDtable.datetime_start = repmat(timeStart, nRows, 1);
        PSDtable.datetime_end   = repmat(timeEnd, nRows, 1);

        PSD_data = [PSD_data; PSDtable];

        % Compute PSD slope
        valid = countsPerL > 0;
        x = log10(binCenters(valid));
        y = log10(countsPerL(valid));
        if numel(x) > 1
            slope = polyfit(x,y,1);
            slope = slope(1);
        else
            slope = NaN;
        end

        PSD_slopes = [PSD_slopes; table(projectName, cruiseName, castVal, stationVal, slope, timeStart, timeEnd, ...
            'VariableNames', {'Project','Cruise','Cast','Station','Slope','DatetimeStart','DatetimeEnd'})];
    end
end

% Add depth info
PSD_data.depthStart   = repmat(depthRange(1), height(PSD_data), 1);
PSD_data.depthEnd     = repmat(depthRange(2), height(PSD_data), 1);
PSD_slopes.depthStart = repmat(depthRange(1), height(PSD_slopes), 1);
PSD_slopes.depthEnd   = repmat(depthRange(2), height(PSD_slopes), 1);

% Save PSD tables
depthStartStr = regexprep(num2str(depthRange(1)), '[^\w]', '_');
depthEndStr   = regexprep(num2str(depthRange(2)), '[^\w]', '_');

save(fullfile(saveDir, sprintf('PSD_data_%s_%s_%.0f.mat', depthStartStr, depthEndStr, ESDmax)), 'PSD_data');
save(fullfile(saveDir, sprintf('PSD_slopes_%s_%s_%.0f.mat', depthStartStr, depthEndStr, ESDmax)), 'PSD_slopes');

disp('PSD computation completed.');
end
