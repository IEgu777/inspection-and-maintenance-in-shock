function Run_GaAs_Grid_All()
% SINGLE CLICK: uncensored 4000h grid shocks + EM, then 8 policy optimizations.
% Gamma shock shape=5.0, rate=7.5 (mean=0.6667), fixed preset dataset.
% All 15*16 monitoring intervals are retained for EM regardless of L=10.
% L=10 REMAINS the failure threshold in the maintenance policy model.
% All eight optimizations share a single estimated shock-event frequency.
% MATLAB: open this file and press F5.
% Use Run_GaAs_Grid_EM.m alone to obtain just the shock+EM dataset.
rootDir = fileparts(mfilename('fullpath'));
oldDir = pwd;
folderCleanup = onCleanup(@()cd(oldDir)); %#ok<NASGU>
cd(rootDir);
fprintf('\n=== STEP 1: FULL 4000-h trajectories, 25%% random grid shocks; EM ===\n');
Run_GaAs_Grid_EM();
Check_RI_Grid_Before_Run();
fprintf('\n=== STEP 2: Eight policy optimizations (COMMON lambda_hat; TRUE/ESTIMATED shock size) ===\n');
Run_GaAs_Compare_True_Estimated();
fprintf('\n=== STEP 3: MTTF DIAGNOSTICS (existing results only) ===\n');
Run_GaAs_MTTF_Diagnostics();
Audit_GaAs_Dense_Results();
fprintf('\n=== COMPLETE; see GaAs_True_vs_Estimated_Comparison.csv ===\n');
end
