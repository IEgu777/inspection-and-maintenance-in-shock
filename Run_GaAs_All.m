function Run_GaAs_All()
% RUN_GAAS_ALL -- Single-click launcher for Example 2 (GaAs laser).
%
% 1. Unzip the complete GaAs_Cm20_Cc90_Cd40_RI_Dense_Full.zip archive.
% 2. Open this file in MATLAB and click Run (F5).
% 3. All four policies and the comparison report execute automatically.
%
% SEARCH RESOLUTION: 'quick' to test; 'full' for paper grid.
% PHYSICAL/COST PARAMETERS: change GaAs_case_parameters.m ONLY.
RUN_MODE = 'full';  % full mode; Hybrid and Policy 1 use 117x99 grids
%
% All results are written to the same extracted folder.
% A MATLAB installation is required; parallel toolbox is optional.
% No original GaAs CSV is necessary for this policy optimization stage.

assert(any(strcmpi(RUN_MODE, {'quick','full'})), ...
    'RUN_MODE must be quick or full.');

rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir)
    error('Cannot locate this file. Save Run_GaAs_All.m before running.');
end

scriptNames = { ...
    'GaAs_hybrid_fast.m', ...
    'GaAs_policy1_fast.m', ...
    'GaAs_policy2_fast.m', ...
    'GaAs_periodic_fast.m', ...
    'GaAs_compare_results.m'};

expectedFiles = { ...
    'GaAs_hybrid_results.mat', ...
    'GaAs_policy1_results.mat', ...
    'GaAs_policy2_results.mat', ...
    'GaAs_periodic_results.mat', ...
    'GaAs_policy_comparison.mat'};

for k = 1:numel(scriptNames)
    if ~isfile(fullfile(rootDir, scriptNames{k}))
        error('Missing %s in %s. Extract the COMPLETE ZIP archive.', ...
            scriptNames{k}, rootDir);
    end
end

originalFolder = pwd;
oldMode = getenv('GAAS_POLICY_RUN_MODE');
folderRestore = onCleanup(@() cd(originalFolder)); %#ok<NASGU>
modeRestore = onCleanup(@() setenv('GAAS_POLICY_RUN_MODE', oldMode)); %#ok<NASGU>
cd(rootDir);
par = GaAs_case_parameters();
fprintf('\nSHARED GaAs SETTINGS: lambda=%.3g, p=%g, q=%g; costs (Ci,Cm,Cp,Cc,Cd)=(%g,%g,%g,%g,%g)\n', ...
    par.lambda,par.p,par.q,par.Ci,par.Cm,par.Cp,par.Cc,par.Cd);
setenv('GAAS_POLICY_RUN_MODE', lower(RUN_MODE));

logName = fullfile(rootDir, sprintf('GaAs_run_log_%s.txt', lower(RUN_MODE)));
if isfile(logName)
    delete(logName);
end
diary(logName);
diary('on');
diaryStop = onCleanup(@() diary('off')); %#ok<NASGU>

fprintf('\n============================================================\n');
fprintf(' GaAs laser Example 2: ONE-CLICK RUN\n');
fprintf(' Mode: %-6s | 1 time unit = 250 hours\n', upper(RUN_MODE));
fprintf(' Output folder: %s\n', rootDir);
fprintf('============================================================\n');

wholeTimer = tic;
for k = 1:numel(scriptNames)
    fprintf('\n\n>>>>> STAGE %d / %d: %s <<<<<\n', ...
        k, numel(scriptNames), scriptNames{k});
    stageTimer = tic;

    % Never read results from an older run if a stage fails midway.
    previousResult = fullfile(rootDir, expectedFiles{k});
    if isfile(previousResult)
        delete(previousResult);
    end

    try
        % The original scripts call clearvars/clc/close all. Run them in
        % MATLAB's BASE workspace so they cannot clear the launcher's
        % function-local state and for-loop variables.
        scriptPath = fullfile(rootDir, scriptNames{k});
        escapedPath = strrep(scriptPath, '''', '''''');
        evalin('base', ['run(''', escapedPath, ''');']);
    catch ME
        fprintf(2, '\nFAILED at stage %d: %s\n', k, scriptNames{k});
        fprintf(2, 'Error: %s\n', ME.message);
        fprintf(2, 'Details and earlier outputs are in %s\n', rootDir);
        rethrow(ME);
    end

    if ~isfile(previousResult)
        error('%s finished without creating %s.', ...
            scriptNames{k}, expectedFiles{k});
    end

    % Check the ACTUAL MAT output, not only the source text.
    if k <= 2 && strcmpi(RUN_MODE,'full')
        Q = load(previousResult,'results');
        gRI = Q.results.coarseGrid.RIGrid(:)';
        expectedRI = [0.01:0.01:0.98,0.981:0.001:0.999];
        assert(numel(gRI)==117 && max(abs(gRI-expectedRI))<1e-11, ...
            'Actual results do not contain the 117-point dense RI grid: %s',previousResult);
        costs = Q.results.parameters;
        assert(isequal([costs.Ci,costs.Cm,costs.Cp,costs.Cc,costs.Cd], ...
            [3,20,50,90,40]),'Cost parameters mismatch in %s',previousResult);
        fprintf('ACTUAL RESULT VERIFIED: %s | 117 RI points, 99 M points\n', ...
            scriptNames{k});
    end
    fprintf('<<<<< Completed stage %d in %.1f s >>>>>\n', ...
        k, toc(stageTimer));
end

comparisonFile = fullfile(rootDir, 'GaAs_policy_comparison.mat');
S = load(comparisonFile, 'comparison');
comparison = S.comparison;
try
    sortedComparison = sortrows(comparison, 'Cinf_per_250h');
catch
    sortedComparison = comparison;
end

fprintf('\n============================================================\n');
fprintf(' ALL FOUR GAAS POLICIES COMPLETED SUCCESSFULLY\n');
fprintf(' Total elapsed time: %.1f s\n', toc(wholeTimer));
fprintf(' Summary (cost rate in cost units per 250 hours):\n');
fprintf('============================================================\n');
disp(sortedComparison);
fprintf('Comparison CSV : %s\n', fullfile(rootDir,'GaAs_policy_comparison.csv'));
fprintf('Comparison MAT : %s\n', comparisonFile);
fprintf('Execution log  : %s\n', logName);
fprintf('Note: Grid optima are numerical approximations.\n');
if any(~isfinite(comparison.MTTF_250h)) || any(~isfinite(comparison.Availability))
    warning('Some metrics are NaN/Inf. Results require numerical investigation.');
end
if any(abs(comparison.RI-0.01)<1e-12 | abs(comparison.RI-0.999)<1e-12)
    warning('At least one optimum is on an RI grid boundary; interpret carefully.');
end
end
