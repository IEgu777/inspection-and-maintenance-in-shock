function Run_GaAs_Compare_True_Estimated()
% Compare TRUE injected vs EM-ESTIMATED shock MAGNITUDE parameters.
% In BOTH groups, the shock arrival intensity is set to the SAME lambda_hat
% estimated from observed eligible grid intervals; lambda=0.25 is NOT
% passed to the TRUE-shock-magnitude group.
% First run Run_GaAs_Grid_EM.m to create GaAs_Grid_EM_Summary.csv.
% To perform both tasks, run Run_GaAs_Grid_All.m instead.
%
% NOTE: Policy cost from TRUE and ESTIMATED optimization is computed under
% DIFFERENT assumed stochastic models. The difference between the two minima
% is NOT the true-model regret of using an estimated-parameter policy.

rootDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
oldSource = getenv('GAAS_POLICY_PARAMETER_SOURCE');
restoreDir = onCleanup(@()cd(oldDir)); %#ok<NASGU>
restoreEnv = onCleanup(@()setenv('GAAS_POLICY_PARAMETER_SOURCE',oldSource)); %#ok<NASGU>
cd(rootDir);
if ~isfile('GaAs_Grid_EM_Summary.csv')
    error('Run Run_GaAs_Grid_EM.m first to generate parameter estimates.');
end

runFiles = {'GaAs_hybrid_results.mat', 'GaAs_policy1_results.mat', ...
    'GaAs_policy2_results.mat', 'GaAs_periodic_results.mat', ...
    'GaAs_policy_comparison.mat', 'GaAs_policy_comparison.csv'};
for k=1:numel(runFiles)
    if isfile(runFiles{k})
        delete(runFiles{k});
    end
end

scenarios = {'true','estimated'};
allTables = cell(2,1);
summary = readtable('GaAs_Grid_EM_Summary.csv');
for is=1:2
    source = scenarios{is};
    fprintf('\n=======================================================\n');
    fprintf(' Shock-magnitude parameters %d/2: %s (COMMON estimated lambda)\n',is,upper(source));
    fprintf('=======================================================\n');
    setenv('GAAS_POLICY_PARAMETER_SOURCE',source);
    P = GaAs_case_parameters();
    fprintf('lambda=%.8f, shock shape=%.8f, shock RATE=%.8f\n', ...
        P.lambda,P.alphaS,P.betaS);
    Run_GaAs_All();

    folder = fullfile(rootDir,['Results_',source]);
    if ~isfolder(folder)
        mkdir(folder);
    end
    for k=1:numel(runFiles)
        f = runFiles{k};
        assert(isfile(fullfile(rootDir,f)), 'Missing output: %s',f);
        copyfile(fullfile(rootDir,f),fullfile(folder,f),'f');
    end
    logFiles = dir(fullfile(rootDir,'GaAs_run_log_*.txt'));
    for k=1:numel(logFiles)
        copyfile(fullfile(rootDir,logFiles(k).name), ...
            fullfile(folder,logFiles(k).name),'f');
    end

    T = readtable('GaAs_policy_comparison.csv');
    assert(height(T)==4,'Expected exactly four policy results.');
    T.ParameterSet = repmat({source},height(T),1);
    T.Lambda = repmat(P.lambda,height(T),1);
    T.ShockShape = repmat(P.alphaS,height(T),1);
    T.ShockRate = repmat(P.betaS,height(T),1);
    % Label computed costs as MODEL-PREDICTED for each parameter set.
    allTables{is}=T;
    fprintf('Saved four-policy %s comparison in %s\n',source,folder);
end

assert(abs(allTables{1}.Lambda(1)-allTables{2}.Lambda(1)) < 1e-12, ...
    'Both optimization groups must use the SAME estimated lambda.');
combined=[allTables{1};allTables{2}];
combined=movevars(combined, ...
    {'ParameterSet','Lambda','ShockShape','ShockRate'},'Before','Policy');
writetable(combined,'GaAs_True_vs_Estimated_Comparison.csv');
save('GaAs_True_vs_Estimated_Comparison.mat','combined','summary');

% Pairwise differences in decisions and predicted minimum costs (not regret).
A=allTables{1}; B=allTables{2};
paired=table(A.Policy,A.RI,B.RI,A.M,B.M,A.T_in_250h,B.T_in_250h, ...
    A.Cinf_per_250h,B.Cinf_per_250h, ...
    B.Cinf_per_250h-A.Cinf_per_250h, ...
    'VariableNames',{'Policy','RI_True','RI_Estimated', ...
    'M_True','M_Estimated','T_True','T_Estimated', ...
    'Cost_TrueModel','Cost_EstimatedModel','CostModelDifference'});
writetable(paired,'GaAs_True_vs_Estimated_Paired.csv');

fprintf('\n=======================================================\n');
fprintf(' Finished 8 optimizations: 4 TRUE + 4 ESTIMATED\n');
fprintf(' COMMON estimated lambda used by BOTH groups = %.6f\n', ...
    summary.Lambda_Estimated_Frequency(1));
fprintf(' Bernoulli injection probability (data-design only) = %.6f\n', ...
    summary.ShockProbability_Set(1));
fprintf(' HPP implied from grid event indicators = %.6f (NOT used)\n', ...
    summary.Lambda_HPP_Implied_From_Indicators(1));
fprintf(' Main comparison: GaAs_True_vs_Estimated_Comparison.csv\n');
fprintf(' Paired table   : GaAs_True_vs_Estimated_Paired.csv\n');
fprintf(' Subfolders: Results_true and Results_estimated\n');
fprintf(' WARNING: differences in minimum model-predicted costs are NOT regret.\n');
fprintf('=======================================================\n');
disp(combined);
end
