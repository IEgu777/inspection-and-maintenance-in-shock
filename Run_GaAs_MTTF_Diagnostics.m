function Run_GaAs_MTTF_Diagnostics()
% Print key diagnostics of MTTF at each CURRENT optimized decision.
% Does not recompute MTTF or costs.
rootDir=fileparts(mfilename('fullpath'));
files={'GaAs_hybrid_results.mat','GaAs_policy1_results.mat', ...
       'GaAs_policy2_results.mat','GaAs_periodic_results.mat'};
policies={'Hybrid','Policy 1','Policy 2','Periodic'};
folders={fullfile(rootDir,'Results_true'),fullfile(rootDir,'Results_estimated')};
labels={'true_magnitude','estimated_magnitude'};
if ~isfile(fullfile(folders{1},files{1})) || ...
   ~isfile(fullfile(folders{2},files{1}))
    folders={rootDir};labels={'single_run'};
end
for a=1:numel(folders)
    fprintf('\n========== MTTF DIAGNOSTICS: %s ==========\n',labels{a});
    for k=1:4
        fn=fullfile(folders{a},files{k});
        if ~isfile(fn), fprintf('%s: file missing\n',policies{k});continue;end
        S=load(fn,'results');
        r=S.results.finalBest;
        d=r.diagnostics;
        rc=NaN;
        if isfield(d,'minMTTFLinearMatrixReciprocalCondition')
            rc=d.minMTTFLinearMatrixReciprocalCondition;
        end
        quad=NaN;
        if isfield(d,'maxUncorrectedComponentError')
            quad=d.maxUncorrectedComponentError;
        elseif isfield(d,'maxRawComponentError')
            quad=d.maxRawComponentError;
        end
        fprintf('%-9s: MTTF_new=%9.3f h  MTTF_stationary=%9.3f h  Avail=%.6f\n', ...
            policies{k},r.MTTF_fromNew*250,r.MTTF_stationary*250,r.averageAvailability);
        fprintf('            cycle=%8.4f x 250 h  down/cycle=%9.5f  rcond=%9.2e  rawMassErr=%9.2e\n', ...
            r.expectedCycleLength,r.expectedDowntimePerCycle,rc,quad);
        if ~isfinite(r.MTTF_fromNew) || r.MTTF_fromNew<=0 || ...
           ~isfinite(rc) || rc < 1e-10
            warning('%s: MTTF numerical validity needs investigation.',policies{k});
        end
    end
end
fprintf('\nNOTE: Different policies have different first-failure distributions.\n');
fprintf('MTTF can be much lower for Policy 2 because only shocks trigger inspections.\n');
end
