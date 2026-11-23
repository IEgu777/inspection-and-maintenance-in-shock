%GAAS_COMPARE_RESULTS Compare four GaAs inspection/maintenance policies.
% Run after all four GaAs_*.m optimization programs have completed.
% One time unit = 250 hours; C_inf is in cost units per 250 hours.
clearvars; clc;

names  = {'Hybrid'; 'Policy 1'; 'Policy 2'; 'Periodic'};
files  = {'GaAs_hybrid_results.mat'; 'GaAs_policy1_results.mat'; ...
          'GaAs_policy2_results.mat'; 'GaAs_periodic_results.mat'};
RI = nan(4,1); M=nan(4,1); T=nan(4,1);
CostRate=nan(4,1); Availability=nan(4,1); MTTFnew=nan(4,1);
GridMode=cell(4,1);
for k=1:4
    if ~isfile(files{k})
        error('Missing %s; first run the corresponding GaAs_*.m script.',files{k});
    end
    loaded = load(files{k},'results');
    r = loaded.results;
    assert(strcmp(r.caseName,'GaAs laser, Example 2'), ...
        'Loaded results are not from this GaAs case.');
    assert(r.timeUnitHours==250,'Inconsistent time units.');
    expected = GaAs_case_parameters();
    fields = fieldnames(expected);
    for j = 1:numel(fields)
        f = fields{j};
        assert(isfield(r.parameters,f) && ...
            abs(r.parameters.(f)-expected.(f)) <= 1e-9 * max(1,abs(expected.(f))), ...
            'Different scenarios/parameters detected: %s',f);
    end
    GridMode{k}=r.runMode;
    CostRate(k)=r.finalBest.costRate;
    Availability(k)=r.finalBest.averageAvailability;
    MTTFnew(k)=r.finalBest.MTTF_fromNew;
    M(k)=r.finalBest.policy.M;
    if isfield(r.finalBest.policy,'RI')
        RI(k)=r.finalBest.policy.RI;
    end
    if isfield(r.finalBest.policy,'T')
        T(k)=r.finalBest.policy.T;
    end
end
if numel(unique(GridMode))>1
    warning('Results mix quick and full grid modes; rerun in consistent modes.');
end
M_over_L=M/10;
MTTF_hours=MTTFnew*250;
Cinf_per_hour=CostRate/250;
comparison=table(names,GridMode,RI,M,M_over_L,T,CostRate, ...
    Cinf_per_hour,Availability,MTTFnew,MTTF_hours, ...
    'VariableNames',{'Policy','GridMode','RI','M','M_over_L', ...
    'T_in_250h','Cinf_per_250h','Cinf_per_hour', ...
    'Availability','MTTF_250h','MTTF_hours'});
disp(comparison);
if any(~isfinite(CostRate)) || any(~isfinite(Availability)) || any(~isfinite(MTTFnew))
    warning(['Some metrics are NaN/Inf; this run is NOT fully validated. ', ...
        'Inspect GaAs_hybrid_results.mat diagnostics before reporting.']);
end
writetable(comparison,'GaAs_policy_comparison.csv');
save('GaAs_policy_comparison.mat','comparison');
fprintf('Saved GaAs_policy_comparison.csv and .mat\n');
