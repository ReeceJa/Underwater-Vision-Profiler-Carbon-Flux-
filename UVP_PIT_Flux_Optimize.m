function [optParams, jackknifeResults, fluxStats, ...
    uvpCastFlux, uvpBinFlux, uvpCruiseFlux] = ...
UVP_PIT_Flux_Optimize( ...
psdFile, pitFile, excludeCruises, initParams, POC_factor, ...
minESD, maxESD)

%% ========================================================================
% UVP_PIT_Flux_Optimize_Jackknife_HOT
%
% Optimize:
%
%       mass = A * d^b
%
% against PIT carbon flux.
%
% The UVP flux is calculated at three levels:
%
%   1. Particle-size-bin × cast × cruise
%   2. Cast × cruise
%   3. Cruise
%
% The optimization is performed at the CRUISE level.
%
% For each cruise:
%
%   - UVP flux is calculated for every cast.
%   - Cast UVP fluxes are averaged to obtain one cruise-level UVP flux.
%   - The cruise-level UVP flux is compared with the cruise-level PIT flux.
%
% A classical leave-one-cruise-out jackknife is then performed.
%
% Each jackknife iteration removes ONE cruise and re-optimizes A and b.
%
% Outputs:
%
%   optParams
%       [A b] from the full dataset
%
%   jackknifeResults
%       Table containing jackknife A and b estimates
%
%   fluxStats
%       Summary statistics including jackknife standard errors
%
%   uvpCastFlux
%       UVP flux for every cast using the full-data optimized parameters
%
%   uvpBinFlux
%       UVP flux contribution from every particle-size bin for every cast
%
%   uvpCruiseFlux
%       Cruise-level UVP flux and corresponding PIT flux
%
%% ========================================================================

%% ========================================================================
% INITIALIZE ALL OUTPUTS
% ========================================================================

optParams = [NaN NaN];

jackknifeResults = table();

fluxStats = table();

% Empty cast-level table
uvpCastFlux = table( ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', {'cruise','cast','uvp_flux'});

% Empty bin-level table
uvpBinFlux = table( ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', { ...
        'cruise', ...
        'cast', ...
        'bin_center_mm', ...
        'bin_width_mm', ...
        'number_per_L_per_mm', ...
        'particle_mass', ...
        'uvp_bin_flux'});

% Empty cruise-level table
uvpCruiseFlux = table( ...
    strings(0,1), ...
    zeros(0,1), ...
    zeros(0,1), ...
    'VariableNames', {'cruise','uvp_flux','pit_flux'});
%% ------------------------------------------------------------------------
% LOAD PSD DATA
% -------------------------------------------------------------------------

S = load(psdFile);

vars = fieldnames(S);

allResults = S.(vars{1});

allResults.cruise = string(allResults.cruise);


%% ------------------------------------------------------------------------
% REMOVE EXCLUDED CRUISES
% -------------------------------------------------------------------------

if ~isempty(excludeCruises)

    excludeCruises = string(excludeCruises);

    isExcluded = ismember(allResults.cruise, excludeCruises);

    allResults = allResults(~isExcluded,:);

end


%% ------------------------------------------------------------------------
% LOAD PIT DATA
% -------------------------------------------------------------------------

PIT_data = readtable(pitFile);

if ismember('crn', PIT_data.Properties.VariableNames)

    PIT_data.Properties.VariableNames{'crn'} = 'cruise';

end

if ismember('carb', PIT_data.Properties.VariableNames)

    PIT_data.Properties.VariableNames{'carb'} = 'carb_flux';

end

PIT_data.cruise = string(PIT_data.cruise);


%% ------------------------------------------------------------------------
% FIND COMMON CRUISES
% -------------------------------------------------------------------------

uvpCruises = unique(allResults.cruise);

pitCruises = unique(PIT_data.cruise);

commonCruises = intersect(uvpCruises, pitCruises);

commonCruises = setdiff(commonCruises, string(excludeCruises));

fprintf('\nCommon UVP/PIT cruises used for optimization:\n');

disp(commonCruises);

fprintf('Number of cruises used: %d\n', numel(commonCruises));

if numel(commonCruises) < 3

    error(['Fewer than 3 common UVP/PIT cruises are available. ' ...
        'Optimization/jackknife cannot be performed reliably.']);

end


%% ------------------------------------------------------------------------
% OPTIMIZATION FUNCTION
% -------------------------------------------------------------------------

objFun = @(params, cruisesToUse) ...
    cruiseFluxError( ...
    params, ...
    allResults, ...
    PIT_data, ...
    cruisesToUse, ...
    POC_factor, ...
    maxESD, ...
    minESD);


%% ------------------------------------------------------------------------
% FULL-DATA OPTIMIZATION
% -------------------------------------------------------------------------

fprintf('\n========================================\n');
fprintf('FULL-DATA UVP-PIT OPTIMIZATION\n');
fprintf('========================================\n');

options = optimset('Display','iter');

[optParams, ~] = fminsearch( ...
    @(params) objFun(params, commonCruises), ...
    initParams, ...
    options);

fprintf('\nOptimized parameters:\n');

fprintf('A = %.10g\n', optParams(1));

fprintf('b = %.10g\n', optParams(2));


%% ------------------------------------------------------------------------
% JACKKNIFE
%
% Leave one cruise out at a time.
% -------------------------------------------------------------------------

nCruises = numel(commonCruises);

jackA = nan(nCruises,1);

jackB = nan(nCruises,1);

jackCruise = strings(nCruises,1);

fprintf('\n========================================\n');
fprintf('LEAVE-ONE-CRUISE-OUT JACKKNIFE\n');
fprintf('========================================\n');

for j = 1:nCruises

    omittedCruise = commonCruises(j);

    jackCruise(j) = omittedCruise;

    cruisesJack = commonCruises;

    cruisesJack(j) = [];

    fprintf('\nJackknife %d of %d\n', j, nCruises);

    fprintf('Leaving out cruise %s\n', omittedCruise);

    [jackParams, ~] = fminsearch( ...
        @(params) objFun(params, cruisesJack), ...
        optParams, ...
        optimset('Display','off'));

    jackA(j) = jackParams(1);

    jackB(j) = jackParams(2);

    fprintf('A = %.10g\n', jackA(j));

    fprintf('b = %.10g\n', jackB(j));

end


%% ------------------------------------------------------------------------
% JACKKNIFE STATISTICS
%
% Classical delete-one jackknife standard error:
%
% SE = sqrt((n-1)/n * sum((theta_i - theta_bar)^2))
% -------------------------------------------------------------------------

meanJackA = mean(jackA,'omitnan');

meanJackB = mean(jackB,'omitnan');

jackSE_A = sqrt( ...
    (nCruises-1)/nCruises * ...
    sum((jackA - meanJackA).^2,'omitnan'));

jackSE_b = sqrt( ...
    (nCruises-1)/nCruises * ...
    sum((jackB - meanJackB).^2,'omitnan'));


%% ------------------------------------------------------------------------
% JACKKNIFE BIAS ESTIMATE
% -------------------------------------------------------------------------

biasA = ...
    (nCruises-1) * ...
    (meanJackA - optParams(1));

biasB = ...
    (nCruises-1) * ...
    (meanJackB - optParams(2));


%% ------------------------------------------------------------------------
% JACKKNIFE RESULTS TABLE
% -------------------------------------------------------------------------

jackknifeResults = table( ...
    jackCruise, ...
    jackA, ...
    jackB, ...
    'VariableNames', ...
    {'OmittedCruise','A','b'});


%% ========================================================================
% CALCULATE UVP FLUX AT BIN, CAST, AND CRUISE LEVELS
% ========================================================================
%
% IMPORTANT:
%
% uvpBinFlux:
%   One row for every particle-size bin in every cast.
%
% uvpCastFlux:
%   One row for every cast.
%
% uvpCruiseFlux:
%   One row for every cruise.
%
% The optimized A and b from the FULL dataset are used here.
%
%% ========================================================================

%% ------------------------------------------------------------------------
% CALCULATE BIN AND CAST FLUX
% -------------------------------------------------------------------------

for i = 1:nCruises

    cruiseName = commonCruises(i);

    cruiseData = allResults( ...
        allResults.cruise == cruiseName,:);

    casts = unique(cruiseData.cast);


    for c = 1:numel(casts)

        castData = cruiseData( ...
            cruiseData.cast == casts(c),:);

        castData = sortrows( ...
            castData,'bin_center_mm');


        %% ---------------------------------------------------------------
        % PARTICLE-SIZE INFORMATION
        % ---------------------------------------------------------------

        d = castData.bin_center_mm(:);

        delta_d = castData.bin_width_mm(:);

        n_d = ...
            castData.counts_per_liter_per_mm(:) * 1000;


        %% ---------------------------------------------------------------
        % APPLY FLUX SIZE RANGE
        % ---------------------------------------------------------------

        validSize = ...
            d >= minESD & ...
            d <= maxESD;

        d = d(validSize);

        delta_d = delta_d(validSize);

        n_d = n_d(validSize);


        %% ---------------------------------------------------------------
        % CALCULATE PARTICLE MASS
        %
        % mass = A * d^b
        % ---------------------------------------------------------------

        mw_d = ...
            optParams(1) .* ...
            (d .^ optParams(2));


        %% ---------------------------------------------------------------
        % CALCULATE FLUX FOR EACH PARTICLE-SIZE BIN
        % ---------------------------------------------------------------
        %
        % n_d:
        %       number / L / mm
        %
        % mw_d:
        %       particle mass
        %
        % delta_d:
        %       bin width
        %
        % POC_factor:
        %       converts particle mass to carbon mass
        %
        % Therefore:
        %
        % binFlux = n_d * mass * bin width * POC_factor
        %
        % The existing code uses:
        %
        %       n_d = counts_per_liter_per_mm * 1000
        %
        % which is retained here unchanged.
        % ---------------------------------------------------------------

        binFlux = ...
            n_d .* ...
            mw_d .* ...
            delta_d .* ...
            POC_factor;


        %% ---------------------------------------------------------------
        % STORE BIN-LEVEL FLUX
        % ---------------------------------------------------------------

        binTable = table( ...
            repmat(cruiseName,numel(d),1), ...
            repmat(casts(c),numel(d),1), ...
            d, ...
            delta_d, ...
            n_d, ...
            mw_d, ...
            binFlux, ...
            'VariableNames', ...
            { ...
            'cruise', ...
            'cast', ...
            'bin_center_mm', ...
            'bin_width_mm', ...
            'number_per_L_per_mm', ...
            'particle_mass', ...
            'uvp_bin_flux'});


        uvpBinFlux = ...
            [uvpBinFlux; binTable];


        %% ---------------------------------------------------------------
        % SUM BIN FLUXES TO OBTAIN TOTAL CAST FLUX
        % ---------------------------------------------------------------

        flux = nansum(binFlux);


        castTable = table( ...
            cruiseName, ...
            casts(c), ...
            flux, ...
            'VariableNames', ...
            {'cruise','cast','uvp_flux'});


        uvpCastFlux = ...
            [uvpCastFlux; castTable];

    end

end


%% ========================================================================
% CALCULATE CRUISE-LEVEL UVP FLUX
% ========================================================================
%
% Each cruise receives ONE UVP flux value.
%
% This is calculated as the mean of the cast-level UVP fluxes.
%
% Therefore:
%
%   Cruise 321 with 5 casts = one UVP value
%
%   Cruise 322 with 10 casts = one UVP value
%
% Neither cruise receives extra weighting simply because it has more casts.
%
%% ========================================================================

for i = 1:nCruises

    cruiseName = commonCruises(i);


    %% --------------------------------------------------------------------
    % CAST FLUXES FOR THIS CRUISE
    % ---------------------------------------------------------------------

    castFluxThisCruise = ...
        uvpCastFlux.uvp_flux( ...
        uvpCastFlux.cruise == cruiseName);


    %% --------------------------------------------------------------------
    % MEAN CAST FLUX = CRUISE FLUX
    % ---------------------------------------------------------------------

    cruiseUVP = ...
        mean(castFluxThisCruise,'omitnan');


    %% --------------------------------------------------------------------
    % PIT FLUX FOR THIS CRUISE
    % ---------------------------------------------------------------------

    cruisePIT = ...
        mean( ...
        PIT_data.carb_flux( ...
        PIT_data.cruise == cruiseName), ...
        'omitnan');


    %% --------------------------------------------------------------------
    % STORE CRUISE RESULTS
    % ---------------------------------------------------------------------

    newRow = table( ...
        cruiseName, ...
        cruiseUVP, ...
        cruisePIT, ...
        'VariableNames', ...
        {'cruise','uvp_flux','pit_flux'});


    uvpCruiseFlux = ...
        [uvpCruiseFlux; newRow];

end


%% ------------------------------------------------------------------------
% EXTRACT CRUISE-LEVEL FLUX FOR REGRESSION
% -------------------------------------------------------------------------

uvpFlux = uvpCruiseFlux.uvp_flux;

pitFlux = uvpCruiseFlux.pit_flux;


%% ========================================================================
% LOG-LOG REGRESSION
% ========================================================================

valid = ...
    uvpFlux > 0 & ...
    pitFlux > 0 & ...
    isfinite(uvpFlux) & ...
    isfinite(pitFlux);

uvpPlot = uvpFlux(valid);

pitPlot = pitFlux(valid);

logUVP = log10(uvpPlot);

logPIT = log10(pitPlot);

lm_log = fitlm(logUVP,logPIT);


%% ========================================================================
% CREATE SUMMARY STATISTICS
% ========================================================================

fluxStats = table();

fluxStats.A = optParams(1);

fluxStats.b = optParams(2);

fluxStats.jackknife_mean_A = meanJackA;

fluxStats.jackknife_mean_b = meanJackB;

fluxStats.jackknife_SE_A = jackSE_A;

fluxStats.jackknife_SE_b = jackSE_b;

fluxStats.jackknife_bias_A = biasA;

fluxStats.jackknife_bias_b = biasB;

fluxStats.R2_log_fit = ...
    lm_log.Rsquared.Ordinary;

fluxStats.p_log = ...
    lm_log.Coefficients.pValue(2);

fluxStats.nCruises = nCruises;

fluxStats.minESD = minESD;

fluxStats.maxESD = maxESD;

fluxStats.POC_factor = POC_factor;


%% ========================================================================
% DISPLAY RESULTS
% ========================================================================

fprintf('\n========================================\n');
fprintf('OPTIMIZATION RESULTS\n');
fprintf('========================================\n');

fprintf('\nFull-data parameters:\n');

fprintf('A = %.10g\n', ...
    optParams(1));

fprintf('b = %.10g\n', ...
    optParams(2));


fprintf('\nJackknife results:\n');

fprintf('A = %.10g +/- %.10g\n', ...
    optParams(1), ...
    jackSE_A);

fprintf('b = %.10g +/- %.10g\n', ...
    optParams(2), ...
    jackSE_b);


fprintf('\nJackknife mean estimates:\n');

fprintf('A = %.10g\n', ...
    meanJackA);

fprintf('b = %.10g\n', ...
    meanJackB);


fprintf('\nJackknife bias estimates:\n');

fprintf('A bias = %.10g\n', ...
    biasA);

fprintf('b bias = %.10g\n', ...
    biasB);


fprintf('\nFlux outputs:\n');

fprintf('Number of bin-level flux rows: %d\n', ...
    height(uvpBinFlux));

fprintf('Number of cast-level flux rows: %d\n', ...
    height(uvpCastFlux));

fprintf('Number of cruise-level flux rows: %d\n', ...
    height(uvpCruiseFlux));


end


%% ========================================================================
% HELPER FUNCTION
% ========================================================================

function err = cruiseFluxError( ...
    params, ...
    allResults, ...
    PIT_data, ...
    cruisesToUse, ...
    POC_factor, ...
    maxESD, ...
    minESD)


%% ------------------------------------------------------------------------
% EXTRACT PARAMETERS
% -------------------------------------------------------------------------

A = params(1);

b = params(2);


%% ------------------------------------------------------------------------
% INITIALIZE
% -------------------------------------------------------------------------

uvpMeanFlux = ...
    nan(numel(cruisesToUse),1);

pitFlux = ...
    nan(numel(cruisesToUse),1);


%% ------------------------------------------------------------------------
% PREVENT NONPHYSICAL PARAMETER COMBINATIONS
% -------------------------------------------------------------------------

if A <= 0 || b <= 0

    err = 1e100;

    return

end


%% ------------------------------------------------------------------------
% LOOP THROUGH CRUISES
% -------------------------------------------------------------------------

for i = 1:numel(cruisesToUse)

    cruiseName = cruisesToUse(i);


    %% --------------------------------------------------------------------
    % GET UVP DATA FOR THIS CRUISE
    % ---------------------------------------------------------------------

    cruiseData = allResults( ...
        allResults.cruise == cruiseName,:);


    if isempty(cruiseData)

        continue

    end


    %% --------------------------------------------------------------------
    % GET CASTS
    % ---------------------------------------------------------------------

    casts = unique(cruiseData.cast);

    fluxes = ...
        nan(numel(casts),1);


    %% --------------------------------------------------------------------
    % CALCULATE FLUX FOR EVERY CAST
    % ---------------------------------------------------------------------

    for c = 1:numel(casts)

        castData = cruiseData( ...
            cruiseData.cast == casts(c),:);


        castData = sortrows( ...
            castData,'bin_center_mm');


        %% ---------------------------------------------------------------
        % PARTICLE SIZE
        % ---------------------------------------------------------------

        d = ...
            castData.bin_center_mm(:);

        delta_d = ...
            castData.bin_width_mm(:);

        n_d = ...
            castData.counts_per_liter_per_mm(:) * 1000;


        %% ---------------------------------------------------------------
        % SIZE RANGE
        % ---------------------------------------------------------------

        validSize = ...
            d >= minESD & ...
            d <= maxESD;

        d = d(validSize);

        delta_d = delta_d(validSize);

        n_d = n_d(validSize);


        %% ---------------------------------------------------------------
        % PARTICLE MASS
        % ---------------------------------------------------------------

        mw_d = ...
            A .* ...
            (d .^ b);


        %% ---------------------------------------------------------------
        % TOTAL CAST FLUX
        % ---------------------------------------------------------------

        fluxes(c) = ...
            nansum( ...
            n_d .* ...
            mw_d .* ...
            delta_d);


        fluxes(c) = ...
            fluxes(c) * POC_factor;

    end


    %% --------------------------------------------------------------------
    % CRUISE UVP FLUX
    %
    % Mean across casts.
    % ---------------------------------------------------------------------

    uvpMeanFlux(i) = ...
        mean(fluxes,'omitnan');


    %% --------------------------------------------------------------------
    % CRUISE PIT FLUX
    % ---------------------------------------------------------------------

    pitFlux(i) = ...
        mean( ...
        PIT_data.carb_flux( ...
        PIT_data.cruise == cruiseName), ...
        'omitnan');

end


%% ========================================================================
% ONLY VALID POSITIVE FLUXES
% ========================================================================

valid = ...
    isfinite(uvpMeanFlux) & ...
    isfinite(pitFlux) & ...
    uvpMeanFlux > 0 & ...
    pitFlux > 0;


if sum(valid) < 2

    err = 1e100;

    return

end


%% ========================================================================
% LOG-SPACE OBJECTIVE
% ========================================================================

log_diff = ...
    log(pitFlux(valid)) - ...
    log(uvpMeanFlux(valid));


err = ...
    sum(log_diff.^2);


end