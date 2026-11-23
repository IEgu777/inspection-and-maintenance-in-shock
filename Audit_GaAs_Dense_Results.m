function Audit_GaAs_Dense_Results()
% Validate both shock-parameter groups AFTER full optimization.
rootDir=fileparts(mfilename('fullpath'));
expectedRI=[0.01:0.01:0.98,0.981:0.001:0.999];
expectedM=(0.01:0.01:0.99)*10;
scenarios={'true','estimated'};
policies={'GaAs_hybrid_results.mat','GaAs_policy1_results.mat'};
rows=cell(4,8);
idx=0;
for a=1:numel(scenarios)
    for b=1:numel(policies)
        idx=idx+1;
        f=fullfile(rootDir,['Results_',scenarios{a}],policies{b});
        assert(isfile(f),'Missing result: %s',f);
        S=load(f,'results');
        assert(isfield(S,'results'),'No results struct in %s',f);
        r=S.results;
        assert(strcmpi(r.runMode,'full'),'Not FULL mode: %s',f);
        ri=r.coarseGrid.RIGrid(:)';
        m=r.coarseGrid.MGrid(:)';
        assert(numel(ri)==117 && max(abs(ri-expectedRI))<1e-11, ...
            'WRONG RI GRID in %s',f);
        assert(numel(m)==99 && max(abs(m-expectedM))<1e-11, ...
            'WRONG M GRID in %s',f);
        P=r.parameters;
        assert(isequal([P.Ci,P.Cm,P.Cp,P.Cc,P.Cd],[3,20,50,90,40]), ...
            'WRONG COSTS in %s',f);
        rows(idx,:)={scenarios{a},policies{b},r.runMode, ...
            numel(ri),numel(m),numel(ri)*numel(m), ...
            r.bestPolicy.RI,r.bestPolicy.M};
        fprintf('ACTUAL RESULTS VERIFIED: %s / %s | FULL 117 x 99 | RI*=%.4f M*=%.4f\n', ...
            scenarios{a},policies{b},r.bestPolicy.RI,r.bestPolicy.M);
    end
end
summary=cell2table(rows,'VariableNames',{'ParameterSet','File', ...
    'GridMode','RI_Count','M_Count','Candidates','Best_RI','Best_M'});
writetable(summary,fullfile(rootDir,'GaAs_RI_Grid_Audit.csv'));
fprintf('PASSED: both groups and both policies use dense full RI grid, original costs.\n');
end
