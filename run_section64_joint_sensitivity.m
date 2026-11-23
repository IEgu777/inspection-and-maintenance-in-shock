function [T, results] = run_section64_joint_sensitivity(cfg)
% RUN_SECTION64_JOINT_SENSITIVITY
% Continuous one-at-a-time sensitivity study for Sections 6.2/6.4.
% Uses the IDENTICAL analytical policy evaluator from the original
% run_shock_cbm_case_v6_fast.m (copied unchanged below, except this file's
% separate optimization driver). The underlying model is not modified.
%
% Outputs: optimal long-run cost, optimal RI and M, elapsed computing time,
% cost-rate surfaces for selected scenarios, and figures/CSV/MAT files.
%
% Usage:
%   [T, results] = run_section64_joint_sensitivity();
%   % Quick smoke test of the workflow (NOT publication-ready):
%   cfg = struct('mode','smoke');
%   [T, results] = run_section64_joint_sensitivity(cfg);
%   % Original 0.01-by-0.01L exhaustive grid (VERY expensive):
%   cfg = struct('searchMode','full');
%   [T, results] = run_section64_joint_sensitivity(cfg);
%
% Default adaptive search: 0.05 global grid, then 0.01 local refinement
% around up to three coarse-grid local minima. This DOES NOT assume
% convexity. To check a reported optimum against the full global grid,
% use searchMode='full' for selected scenarios.
%
% WARNING: 12 parameters x 11 ratios can be a long run. Every completed
% setting is immediately checkpointed to CSV and MAT (including time).
% Different runs use different output directories by search mode.
%
% MATLAB R2019b or later recommended (table, writetable, tiledlayout).
% No Parallel Computing Toolbox is required for this self-contained file.

if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = sensitivityDefaults(cfg);

par = struct('alpha',1.02, 'beta',1.20, 'lambda',0.25, ...
    'alphaS',2.5, 'betaS',1.0, 'p',6, 'q',5, 'L',12.5, ...
    'Ci',10, 'Cm',100, 'Cp',280, 'Cc',350, 'Cd',50);
opts = struct('Nz',60,'Nv',80,'Ny',18,'NtInt',80,'Nconv',60, ...
    'NuRel',18,'NtRel',60,'tRelMaxFactor',2.5, ...
    'computeReliabilityCurve',false,'useParallel',false, ...
    'showProgress',false,'progressEverySeconds',10, ...
    'minMassWarning',5e-2,'warnOnMass',false);

if strcmpi(cfg.mode,'smoke')
    cfg.useParallel=false;
    % Use only for checking code/data export, NEVER for manuscript results.
    cfg.parameterNames = {'lambda','Ci'};
    cfg.ratios = [0.5, 1, 1.5];
    cfg.coarseStep = 0.20;
    cfg.fineStep = 0.05;
    cfg.refineRadius = 0.20;
    cfg.maxBasins = 2;
    opts.Nz=16; opts.Nv=24; opts.Ny=6; opts.NtInt=24; opts.Nconv=25;
end

if ~isempty(cfg.overrideOptions)
    fieldNames = fieldnames(cfg.overrideOptions);
    for i=1:numel(fieldNames)
        opts.(fieldNames{i}) = cfg.overrideOptions.(fieldNames{i});
    end
end

[~, ~] = mkdir(cfg.outputDir);
root = cfg.outputDir;
% Optional parallel worker pool. No dependency if the toolbox is absent.
cfg.parallelActive=false;
if cfg.useParallel && license('test','Distrib_Computing_Toolbox')
    try
        pool=gcp('nocreate');
        if isempty(pool), pool=parpool; end
        cfg.parallelActive=true;
        fprintf(' Parallel pool active: %d workers.\n',pool.NumWorkers);
    catch ME
        warning('Parallel pool unavailable, switching to serial: %s',ME.message);
    end
end
csvFile = fullfile(root,'Table_section64_continuous_sensitivity.csv');
matFile = fullfile(root,'section64_joint_sensitivity_results.mat');

fprintf('\n=========================================================\n');
fprintf(' Section 6.4 continuous sensitivity (original V6 model)\n');
fprintf(' Mode: %-10s  Search: %s\n',cfg.mode,cfg.searchMode);
fprintf(' Ratios: %s\n',mat2str(cfg.ratios,4));
fprintf(' Parameters: %s\n',strjoin(cfg.parameterNames,', '));
fprintf(' Output folder: %s\n',root);
fprintf('=========================================================\n');

records = repmat(emptyRecord(),0,1);
surfaces = struct('Parameter',{},'Ratio',{},'RIGrid',{},'MGrid',{}, ...
    'CostMatrix',{},'RIStar',{},'MStar',{},'CostStar',{});
startAll = tic;

% Resume a long run from the last completed scenario if requested.
resuming=false;
if cfg.resume && isfile(matFile)
    old=load(matFile,'T','surfaces','par','opts','cfg','baseline','baseCache');
    if ~isfield(old,'baseCache') || ...
       ~isequal(old.par,par) || ~isequal(old.opts,opts) || ...
       ~isequal(old.cfg.ratios,cfg.ratios) || ...
       ~isequal(old.cfg.parameterNames,cfg.parameterNames) || ...
       ~strcmpi(old.cfg.searchMode,cfg.searchMode) || ...
       old.cfg.coarseStep~=cfg.coarseStep || ...
       old.cfg.fineStep~=cfg.fineStep || ...
       old.cfg.maxBasins~=cfg.maxBasins || ...
       old.cfg.refineRadius~=cfg.refineRadius
        error(['Existing checkpoint uses different parameters/grid settings. ' ...
            'Choose a new cfg.outputDir or set cfg.resume=false to restart.']);
    end
    records=table2struct(old.T);
    surfaces=old.surfaces;
    baseCache=old.baseCache;
    sol0=old.baseline;
    if isempty(surfaces)
        error('Existing checkpoint has no baseline surface. Restart with cfg.resume=false.');
    end
    surf0=struct('RIGrid',surfaces(1).RIGrid, ...
        'MGrid',surfaces(1).MGrid,'CostMatrix',surfaces(1).CostMatrix);
    resuming=true;
    fprintf(' Resuming from %d saved rows. Baseline will not be recomputed.\n',numel(records));
end
if ~resuming
    % Compute the baseline cost surface/basis ONCE. For cost-only sensitivity,
    % the probability law and expected cycle quantities do NOT depend on costs.
    % Therefore their per-policy numerical evaluations can be reused exactly.
    fprintf('\n[Baseline model] computing reusable policy grid ...\n');
    [sol0,baseCache,surf0] = solveOneSetting(par,opts,cfg,[]);
    if strcmp(sol0.Status,'OK')
        surfaces(end+1,1) = makeSurface('Baseline',1,surf0,sol0); %#ok<AGROW>
        fprintf(' Baseline: RI*=%.4f, M*=%.4f, C*=%.6f, %.1f s\n', ...
            sol0.RIStar,sol0.MStar,sol0.CostStar,sol0.SearchSeconds);
        if strcmpi(cfg.mode,'publication') && ~strcmpi(cfg.searchMode,'full')
            if abs(sol0.CostStar-37.62366)>0.15 || ...
                    abs(sol0.RIStar-0.92)>0.04 || abs(sol0.MStar-9.625)>0.5
                warning(['Baseline result differs from the original Section 6.2 result. ' ...
                    'Check search precision and input parameters before publication.']);
            end
        end
    else
        warning('Baseline failed: %s',sol0.Message);
    end

    % Persist the baseline timing separately. Ratio-1 rows below reuse this
    % evaluation and therefore correctly report zero *additional* search time.
    baseRec=emptyRecord();
    baseRec.Parameter='Baseline';
    baseRec.Multiplier=1;
    baseRec.RIStar=sol0.RIStar;
    baseRec.MStar=sol0.MStar;
    baseRec.MStarOverL=sol0.MStar/par.L;
    baseRec.CostStar=sol0.CostStar;
    baseRec.CostChangePct=0;
    baseRec.NumRequestedEvaluations=sol0.NumRequestedEvaluations;
    baseRec.NumNewEvaluations=sol0.NumNewEvaluations;
    baseRec.NumCacheHits=sol0.NumCacheHits;
    baseRec.SearchSeconds=sol0.SearchSeconds;
    baseRec.TotalSeconds=sol0.SearchSeconds;
    baseRec.ReusedModelCache=false;
    baseRec.Status=sol0.Status;
    baseRec.Message=sol0.Message;
    records(end+1,1)=baseRec;
    T=struct2table(records);
    saveCurrent(csvFile,matFile,T,surfaces,par,opts,cfg,sol0,baseCache,toc(startAll));

end

% First: all non-cost parameters (cost-only cases reuse the base cache).
% These seven parameters correspond to the left panel of existing Fig. 8.
physicalNames = {'beta','betaS','alpha','lambda','alphaS','p','q'};
costNames = {'Cc','Cm','Ci','Cp','Cd'};
orderedNames = [physicalNames, costNames];
orderedNames = orderedNames(ismember(orderedNames,cfg.parameterNames));

for ip = 1:numel(orderedNames)
    pname = orderedNames{ip};
    fprintf('\n====== Sensitivity of %s (baseline = %.6g) ======\n',pname,par.(pname));
    isCost = ismember(pname,costNames);
    for ir = 1:numel(cfg.ratios)
        ratio = cfg.ratios(ir);
        alreadyDone=strcmp({records.Parameter},pname) & ...
                    abs([records.Multiplier]-ratio)<1e-10;
        if any(alreadyDone)
            fprintf(' %s ratio %.3f : already checkpointed, skipping.\n',pname,ratio);
            continue;
        end
        pa = par;
        pa.(pname) = pa.(pname)*ratio;
        fprintf(' %s ratio %.3f : ',pname,ratio);
        wall = tic;
        if abs(ratio-1) <= 1e-10
            sol = sol0;
            sol.SearchSeconds=0;
            sol.NumNewEvaluations=0;
            sol.NumCacheHits=sol0.NumRequestedEvaluations;
            sol.NumRequestedEvaluations=0;
            sol.Message='Reused baseline (elapsed time recorded in baseline row).';
            extraTime=0;
            srf=surf0;
            cacheUsed=true;
        elseif isCost
            [sol,baseCache,srf] = solveOneSetting(pa,opts,cfg,baseCache);
            extraTime=toc(wall);
            cacheUsed=true;
        else
            [sol,~,srf] = solveOneSetting(pa,opts,cfg,[]);
            extraTime=toc(wall);
            cacheUsed=false;
        end
        if abs(ratio-1)<=1e-10, extraTime=0; end
        rec=emptyRecord();
        rec.Parameter=pname;
        rec.Multiplier=ratio;
        rec.ParameterValue=pa.(pname);
        rec.RIStar=sol.RIStar;
        rec.MStar=sol.MStar;
        rec.MStarOverL=sol.MStar/par.L;
        rec.CostStar=sol.CostStar;
        rec.CostChangePct=100*(sol.CostStar/sol0.CostStar-1);
        rec.NumRequestedEvaluations=sol.NumRequestedEvaluations;
        rec.NumNewEvaluations=sol.NumNewEvaluations;
        rec.NumCacheHits=sol.NumCacheHits;
        rec.SearchSeconds=sol.SearchSeconds;
        rec.TotalSeconds=extraTime;
        rec.ReusedModelCache=cacheUsed;
        rec.Status=sol.Status;
        rec.Message=sol.Message;
        records(end+1,1)=rec; %#ok<AGROW>
        if strcmp(sol.Status,'OK')
            fprintf('RI*=%.4f, M*=%.4f, C*=%.6f, %.1f s [new %d; cached %d]\n', ...
                sol.RIStar,sol.MStar,sol.CostStar,rec.TotalSeconds, ...
                rec.NumNewEvaluations,rec.NumCacheHits);
        else
            fprintf('FAILED: %s\n',sol.Message);
        end
        if isSelectedSurface(pname,ratio,cfg)
            surfaces(end+1,1)=makeSurface(pname,ratio,srf,sol); %#ok<AGROW>
        end
        T = struct2table(records);
        saveCurrent(csvFile,matFile,T,surfaces,par,opts,cfg,sol0,baseCache,toc(startAll));
    end
end

T = struct2table(records);
results=struct('baseline',sol0,'parameters',par,'options',opts, ...
    'config',cfg,'surfaces',surfaces,'totalSeconds',toc(startAll));
saveCurrent(csvFile,matFile,T,surfaces,par,opts,cfg,sol0,baseCache,results.totalSeconds);
if ~isempty(surfaces)
    shapeStats=surfaceShapeDiagnostics(surfaces,par);
    writetable(shapeStats,fullfile(root,'Table_section64_surface_checks.csv'));
    results.surfaceShapeDiagnostics=shapeStats;
end

if cfg.makePlots
    try
        plotContinuousResults(T,par,cfg,root);
    catch ME
        warning('Unable to generate continuous-sensitivity figures: %s',ME.message);
    end
    try
        plotCostSurfaces(surfaces,par,root);
    catch ME
        warning('Unable to generate cost-surface figure: %s',ME.message);
    end
end

fprintf('\n=========================================================\n');
fprintf(' Completed. Total elapsed time: %.3f seconds (%.2f hours).\n', ...
    results.totalSeconds,results.totalSeconds/3600);
fprintf(' Results: %s\n',csvFile);
fprintf(' Full results: %s\n',matFile);
fprintf(' IMPORTANT: adaptive results are NOT mathematical proof of convexity.\n');
fprintf('=========================================================\n');
end

function cfg = sensitivityDefaults(cfg)
if ~isfield(cfg,'mode'), cfg.mode='publication'; end
if ~isfield(cfg,'searchMode'), cfg.searchMode='adaptive'; end
if ~isfield(cfg,'ratios'), cfg.ratios=linspace(0.5,1.5,11); end
if ~isfield(cfg,'parameterNames')
    cfg.parameterNames={'beta','betaS','alpha','lambda','alphaS','p','q', ...
                        'Cc','Cm','Ci','Cp','Cd'};
end
if ~isfield(cfg,'coarseStep'), cfg.coarseStep=0.05; end
if ~isfield(cfg,'fineStep'), cfg.fineStep=0.01; end
if ~isfield(cfg,'refineRadius'), cfg.refineRadius=0.05; end
if ~isfield(cfg,'maxBasins'), cfg.maxBasins=3; end
if ~isfield(cfg,'makePlots'), cfg.makePlots=true; end
if ~isfield(cfg,'useParallel'), cfg.useParallel=true; end
if ~isfield(cfg,'resume'), cfg.resume=true; end
if ~isfield(cfg,'overrideOptions'), cfg.overrideOptions=struct(); end
if ~isfield(cfg,'outputDir')
    root=fileparts(mfilename('fullpath'));
    cfg.outputDir=fullfile(root,sprintf('section64_results_%s_%s_%dpts', ...
        lower(cfg.mode),lower(cfg.searchMode),numel(cfg.ratios)));
end
if ~ismember(lower(cfg.searchMode),{'adaptive','full'})
    error('cfg.searchMode must be either adaptive or full.');
end
if ~all(cfg.ratios>0), error('Every parameter multiplier must be positive.'); end
if any(strcmpi(cfg.parameterNames,'L'))
    error('The failure threshold L is fixed in the original Section 6.4 figure.');
end
end

function r=emptyRecord()
r=struct('Parameter','','Multiplier',NaN,'ParameterValue',NaN, ...
    'RIStar',NaN,'MStar',NaN,'MStarOverL',NaN,'CostStar',NaN, ...
    'CostChangePct',NaN,'NumRequestedEvaluations',0, ...
    'NumNewEvaluations',0,'NumCacheHits',0, ...
    'SearchSeconds',NaN,'TotalSeconds',NaN, ...
    'ReusedModelCache',false,'Status','','Message','');
end

function cache=emptyCache()
cache=struct();
cache.keys=containers.Map('KeyType','char','ValueType','double');
cache.basis=zeros(0,5);
cache.err=cell(0,1);
cache.allRI=zeros(0,1);
cache.allM=zeros(0,1);
end

function [sol,cache,surface] = solveOneSetting(par,opts,cfg,cache)
% Search without assuming convexity. Cost-only scenarios may pass a cache
% from the baseline physical model: basis=[1,EPM,EPR,ECR,Edowntime]/Ecycle.
if isempty(cache), cache=emptyCache(); end
start=tic;
newCalls=0; cacheHits=0; totalRequests=0;
sol=struct('RIStar',NaN,'MStar',NaN,'CostStar',NaN, ...
    'NumRequestedEvaluations',0,'NumNewEvaluations',0, ...
    'NumCacheHits',0,'SearchSeconds',NaN,'Status','ERROR','Message','');

% The exact original grid is (RI,M/L) = 0.01,0.02,...,0.99.
if strcmpi(cfg.searchMode,'full')
    rgrid=(0.01:0.01:0.99)';
    mgrid=(0.01:0.01:0.99)'*par.L;
else
    rgrid=(cfg.coarseStep:cfg.coarseStep:1-cfg.coarseStep)';
    mgrid=(cfg.coarseStep:cfg.coarseStep:1-cfg.coarseStep)'*par.L;
end
[rr,mm]=ndgrid(rgrid,mgrid);
[coarseCost,~]=scorePairs(rr(:),mm(:));
C=reshape(coarseCost,numel(rgrid),numel(mgrid));
surface=struct('RIGrid',rgrid,'MGrid',mgrid,'CostMatrix',C);

valid=find(isfinite(C));
if isempty(valid)
    sol.Message='All coarse-grid policy evaluations failed.';
    finish();
    return;
end

[rFine,mFine,Cfine]=deal([]);
if ~strcmpi(cfg.searchMode,'full')
    % Find ALL discrete 8-neighbor local minima, then refine the best few.
    % Do not use a convexity assumption or a single-start local optimizer.
    basinIdx=gridLocalMinima(C);
    [~,ord]=sort(C(basinIdx),'ascend');
    basinIdx=basinIdx(ord(1:min(numel(ord),cfg.maxBasins)));
    [bi,bj]=ind2sub(size(C),basinIdx);
    for ib=1:numel(bi)
        rlo=max(0.01,rgrid(bi(ib))-cfg.refineRadius);
        rhi=min(0.99,rgrid(bi(ib))+cfg.refineRadius);
        mlo=max(0.01,mgrid(bj(ib))/par.L-cfg.refineRadius);
        mhi=min(0.99,mgrid(bj(ib))/par.L+cfg.refineRadius);
        rlocal=fineAxis(rlo,rhi,cfg.fineStep);
        mlocal=fineAxis(mlo,mhi,cfg.fineStep)*par.L;
        [rmat,mmat]=ndgrid(rlocal,mlocal);
        rFine=[rFine;rmat(:)]; %#ok<AGROW>
        mFine=[mFine;mmat(:)]; %#ok<AGROW>
    end
    if ~isempty(rFine)
        points=unique(round([rFine,mFine/par.L]*1e8)/1e8,'rows','stable');
        rFine=points(:,1); mFine=points(:,2)*par.L;
        [Cfine,~]=scorePairs(rFine,mFine);
    end
end

allRI=[rr(:);rFine(:)];
allM=[mm(:);mFine(:)];
allC=[coarseCost(:);Cfine(:)];
val=find(isfinite(allC));
if isempty(val)
    sol.Message='No finite cost rate after refinement.';
else
    [sol.CostStar,k]=min(allC(val));
    k=val(k);
    sol.RIStar=allRI(k);
    sol.MStar=allM(k);
    sol.Status='OK';
    if numel(val)<numel(allC)
        sol.Message=sprintf('%d invalid/failed policies; see MATLAB warnings.',numel(allC)-numel(val));
    else
        sol.Message='OK';
    end
    if abs(sol.RIStar-0.01)<1e-10 || abs(sol.RIStar-0.99)<1e-10 || ...
       abs(sol.MStar/par.L-0.01)<1e-10 || abs(sol.MStar/par.L-0.99)<1e-10
        warning('Optimal setting (RI=%.3f,M/L=%.3f) reached the search boundary.', ...
            sol.RIStar,sol.MStar/par.L);
    end
end
finish();

    function [costs,isNew] = scorePairs(riValues,mValues)
        n=numel(riValues);
        costs=NaN(n,1); isNew=false(n,1);
        co=[par.Ci,par.Cm,par.Cp,par.Cc,par.Cd];
        indexes=zeros(n,1);
        need=[];
        for j=1:n
            totalRequests=totalRequests+1;
            % Every scenario keeps L=12.5. Normalize M by L for stable keys.
            ky=sprintf('%0.8f|%0.8f',riValues(j),mValues(j)/par.L);
            if isKey(cache.keys,ky)
                indexes(j)=cache.keys(ky);
                cacheHits=cacheHits+1;
            else
                newCalls=newCalls+1;
                isNew(j)=true;
                index=size(cache.basis,1)+1;
                cache.keys(ky)=index;
                cache.allRI(index,1)=riValues(j);
                cache.allM(index,1)=mValues(j);
                cache.basis(index,:)=NaN(1,5);
                cache.err{index,1}='Not evaluated';
                indexes(j)=index;
                need(end+1,1)=j; %#ok<AGROW>
            end
        end
        if ~isempty(need)
            nNew=numel(need);
            addedBasis=NaN(nNew,5);
            addedErr=cell(nNew,1);
            newR=riValues(need); newM=mValues(need);
            if cfg.parallelActive && nNew>1
                parfor k=1:nNew
                    [addedBasis(k,:),addedErr{k}]=costBasisAtPolicy(newR(k),newM(k),par,opts);
                end
            else
                for k=1:nNew
                    [addedBasis(k,:),addedErr{k}]=costBasisAtPolicy(newR(k),newM(k),par,opts);
                end
            end
            for k=1:nNew
                idx=indexes(need(k));
                cache.basis(idx,:)=addedBasis(k,:);
                cache.err{idx,1}=addedErr{k};
            end
        end
        costs=cache.basis(indexes,:)*co(:);
    end
    function finish()
        sol.NumRequestedEvaluations=totalRequests;
        sol.NumNewEvaluations=newCalls;
        sol.NumCacheHits=cacheHits;
        sol.SearchSeconds=toc(start);
    end
end

function [basis,message]=costBasisAtPolicy(RI,M,par,opts)
% This is the EXACT original V6 policy objective, decomposed linearly by cost.
basis=NaN(1,5);
message='';
try
    rep=evaluatePolicy(RI,M,par,opts,false,false);
    v=rep.expectedCountsPerCycle;
    basis=[1,v.PM,v.preventiveReplacement,v.correctiveReplacement, ...
        rep.expectedDowntimePerCycle]/rep.expectedCycleLength;
    if any(~isfinite(basis)) || any(basis<0)
        error('Invalid expected cycle/cost basis.');
    end
catch ME
    message=ME.message;
end
end

function a=fineAxis(lo,hi,step)
% Explicitly align fine points to a COMMON 0.01 grid, not to moving centers.
% This makes cached cost vectors reusable and policy precision comparable.
first=ceil((lo-1e-12)/step);
last=floor((hi+1e-12)/step);
a=(first:last)'*step;
a=a(a>=0.01-1e-12 & a<=0.99+1e-12);
end

function idx=gridLocalMinima(C)
[nr,nc]=size(C);
idx=[];
for i=1:nr
    for j=1:nc
        if ~isfinite(C(i,j)), continue; end
        near=C(max(1,i-1):min(nr,i+1),max(1,j-1):min(nc,j+1));
        finiteNear=near(isfinite(near));
        if isempty(finiteNear), continue; end
        if C(i,j)<=min(finiteNear)+1e-12
            idx(end+1,1)=sub2ind(size(C),i,j); %#ok<AGROW>
        end
    end
end
if isempty(idx)
    good=find(isfinite(C));
    if ~isempty(good)
        [~,j]=min(C(good)); idx=good(j);
    end
end
end

function tf=isSelectedSurface(name,ratio,cfg)
% Use the SAME scenarios as the continuous sensitivity study.
if strcmpi(cfg.mode,'smoke')
    tf=strcmp(name,'lambda');
else
    tf=ismember(name,{'lambda','beta'}) && ...
       (abs(ratio-0.5)<1e-10 || abs(ratio-1.5)<1e-10);
end
end

function s=makeSurface(name,ratio,grid,sol)
s=struct('Parameter',name,'Ratio',ratio,'RIGrid',grid.RIGrid, ...
    'MGrid',grid.MGrid,'CostMatrix',grid.CostMatrix, ...
    'RIStar',sol.RIStar,'MStar',sol.MStar,'CostStar',sol.CostStar);
end

function out=surfaceShapeDiagnostics(surfaces,par)
% DISCRETE diagnostics on the sampled grid ONLY. They cannot prove
% general convexity. Negative second differences can indicate nonconvex
% behavior, but small negatives may also be quadrature/grid artifacts.
rows=cell(numel(surfaces),1);
for k=1:numel(surfaces)
    s=surfaces(k);
    C=s.CostMatrix;
    finiteC=C(isfinite(C));
    nMin=numel(gridLocalMinima(C));
    r2=diff(C,2,1);
    m2=diff(C,2,2);
    r2=r2(isfinite(r2));
    m2=m2(isfinite(m2));
    if isempty(r2), minR=NaN; negR=NaN; else, minR=min(r2); negR=sum(r2 < -1e-4); end
    if isempty(m2), minM=NaN; negM=NaN; else, minM=min(m2); negM=sum(m2 < -1e-4); end
    if isempty(finiteC), cmin=NaN; else, cmin=min(finiteC); end
    rows{k}=struct('Parameter',s.Parameter,'Multiplier',s.Ratio, ...
        'GridPointsRI',numel(s.RIGrid),'GridPointsM',numel(s.MGrid), ...
        'DiscreteLocalMinima',nMin,'MinSecondDiffRI',minR, ...
        'MinSecondDiffM',minM,'NegativeSecondDiffCountRI',negR, ...
        'NegativeSecondDiffCountM',negM,'GridMinimumCost',cmin);
end
out=struct2table(vertcat(rows{:}));
end

function saveCurrent(csvFile,matFile,T,surfaces,par,opts,cfg,baseline,baseCache,elapsed)
writetable(T,csvFile);
% Keep a complete checkpoint after every finished setting. Use -v7.3
% because all surface matrices may exceed the standard MAT size limit.
save(matFile,'T','surfaces','par','opts','cfg','baseline','baseCache','elapsed','-v7.3');
end

function plotContinuousResults(T,par,cfg,root)
% Each panel shows all parameter paths; original paper splits physical and
% cost parameters. Baseline ratios and outputs are NOT normalized in y.
model={'beta','betaS','alpha','lambda','alphaS','p','q'};
price={'Cc','Cm','Ci','Cp','Cd'};
plotOne(T,model,price,par,cfg,root,'CostStar','Optimal cost rate', ...
    'Figure_section64_cost_rate');
plotOne(T,model,price,par,cfg,root,'RIStar','Optimal inspection threshold R_I^*', ...
    'Figure_section64_RIstar');
plotOne(T,model,price,par,cfg,root,'MStar','Optimal maintenance threshold M^*', ...
    'Figure_section64_Mstar');
end

function plotOne(T,model,price,par,cfg,root,varName,yName,fileName) %#ok<INUSD>
f=figure('Color','w','Position',[100,100,1150,425], 'Visible','off');
cols={model,price}; panels={'Degradation, shock and PM parameters','Cost parameters'};
for side=1:2
    subplot(1,2,side); hold on; grid on; box on;
    names=cols{side};
    for k=1:numel(names)
        mask=strcmp(T.Parameter,names{k}) & strcmp(T.Status,'OK');
        if any(mask)
            xx=T.Multiplier(mask); yy=T.(varName)(mask);
            [xx,order]=sort(xx);
            plot(xx,yy(order),'-o','LineWidth',1.2,'MarkerSize',3, ...
                'DisplayName',names{k});
        end
    end
    xlabel('Parameter value / baseline value');
    ylabel(yName,'Interpreter','tex');
    title(panels{side});
    legend('show','Location','best','Interpreter','none');
    xlim([min(cfg.ratios),max(cfg.ratios)]);
end
try
    exportgraphics(f,fullfile(root,[fileName,'.png']),'Resolution',200);
catch
    saveas(f,fullfile(root,[fileName,'.png']));
end
savefig(f,fullfile(root,[fileName,'.fig']));
close(f);
end

function plotCostSurfaces(surfaces,par,root)
if isempty(surfaces), return; end
n=numel(surfaces);
ncol=min(3,n); nrow=ceil(n/ncol);
f=figure('Color','w','Position',[80,80,1150,400*nrow],'Visible','off');
for i=1:n
    s=surfaces(i);
    subplot(nrow,ncol,i);
    contourf(s.RIGrid,s.MGrid/par.L,s.CostMatrix',16,'LineStyle','none');
    colorbar; hold on;
    if isfinite(s.RIStar), plot(s.RIStar,s.MStar/par.L,'rx','LineWidth',1.5,'MarkerSize',8); end
    xlabel('R_I'); ylabel('M/L');
    title(sprintf('%s, ratio %.3f',s.Parameter,s.Ratio),'Interpreter','none');
    xlim([min(s.RIGrid),max(s.RIGrid)]);
    ylim([min(s.MGrid),max(s.MGrid)]/par.L);
end
try
    exportgraphics(f,fullfile(root,'Figure_section64_cost_surfaces.png'),'Resolution',200);
catch
    saveas(f,fullfile(root,'Figure_section64_cost_surfaces.png'));
end
savefig(f,fullfile(root,'Figure_section64_cost_surfaces.fig'));
close(f);
end

%% ======================================================================
% Original analytical evaluator and mathematical numerical helpers below.
% Copied from run_shock_cbm_case_v6_fast.m, with NO algebraic changes.
%% ======================================================================
function rep = evaluatePolicy(RI, M, par, opts, doReliabilityCurve, doMTTF)
if nargin < 5
    doReliabilityCurve = false;
end
if nargin < 6
    doMTTF = true;
end
if doReliabilityCurve && ~doMTTF
    error('MTTF is needed to set the automatic reliability-curve horizon.');
end
validatePolicy(RI, M, par);

% Midpoint grid on (0,M), representing the continuous post-maintenance state.
Nz = opts.Nz;
z  = ((1:Nz)' - 0.5) / Nz * M;
wz = ones(Nz,1) * (M / Nz);
states = [0; z];
Ns = numel(states);

% Midpoint grid on [M,L] for maintenance-state integrals.
Nv = opts.Nv;
v  = M + ((1:Nv)' - 0.5) / Nv * (par.L - M);
wv = ones(Nv,1) * ((par.L - M) / Nv);
pmSurvVec = betaSurvRate(M ./ v, par.p, par.q);  % P(Beta >= M/v)
% This matrix does not depend on the initial state s or shock time y.
% It was previously recomputed hundreds of times per policy evaluation.
pmPdfMat = pmPdfMatrix(z, v, par.p, par.q);

% Build kernels and conditional quantities for all starting states.
Kdens      = zeros(Ns, Nz);     % full transition density to interior states
Q1full     = zeros(Ns, 1);      % full transition mass to state 0
Psurv      = zeros(Ns, Ns);     % first-cycle survival transition matrix
PrelSurv   = zeros(Ns, Ns);     % scheduled reliability-trigger survival matrix only

% Diagnostics before and after component-wise mass correction.
PrawMassError       = zeros(Ns,1);
PcorrectedMassError = zeros(Ns,1);
PdoMassError        = zeros(Ns,1);
PpmSuccessMassError = zeros(Ns,1);
PpmEventMassError   = zeros(Ns,1);

cond = initConditionalStruct(Ns);
for i = 1:Ns
    row = buildOneStepRow(states(i), RI, M, par, opts, z, wz, v, wv, pmSurvVec, pmPdfMat, doMTTF);
    Kdens(i,:)       = row.fullInteriorDensity;
    Q1full(i)        = row.fullMass0;
    Psurv(i,:)       = row.PsurvRow;
    PrelSurv(i,:)    = row.PrelSurvRow;

    PrawMassError(i)       = row.rawFullMassError;
    PcorrectedMassError(i) = row.correctedFullMassError;
    PdoMassError(i)        = row.doNothingMassError;
    PpmSuccessMassError(i) = row.pmSuccessMassError;
    PpmEventMassError(i)   = row.pmEventMassError;

    cond.tau(i)      = row.tau;
    cond.mu(i)       = row.mu;
    cond.downtime(i) = row.downtime;
    cond.Ppm(i)      = row.Ppm;
    cond.Ppr(i)      = row.Ppr;
    cond.Pcr(i)      = row.Pcr;
    cond.mttfSource(i) = row.mttfSource;
end

% Solve the stationary distribution in the form pi = a0*delta0 + a0*H(z)dz.
% H(z_i) = K(0,z_i) + integral K(s,z_i)H(s)ds.
A = zeros(Nz,Nz);
for iz = 1:Nz
    for js = 1:Nz
        A(iz,js) = Kdens(js+1,iz) * wz(js);  % row js+1 corresponds to s=z_js
    end
end
b = Kdens(1,:)';
H = (eye(Nz) - A) \ b;
% Small negative values can occur because of numerical quadrature/linear solve.
H(H < 0 & H > -1e-10) = 0;
if any(H < -1e-8)
    warning('Stationary density H contains noticeable negative values. Increase grid sizes.');
end
H = max(H,0);
a0 = 1 / (1 + sum(H .* wz));
piVec = [a0; a0 * H .* wz];
piVec = piVec / sum(piVec);

% Expected cycle quantities under stationary distribution.
Ecycle = sum(piVec .* cond.mu);
Epm    = sum(piVec .* cond.Ppm);
Epr    = sum(piVec .* cond.Ppr);
Ecr    = sum(piVec .* cond.Pcr);
Edown  = sum(piVec .* cond.downtime);

costInspection = par.Ci / Ecycle;
costPM         = par.Cm * Epm / Ecycle;
costPR         = par.Cp * Epr / Ecycle;
costCR         = par.Cc * Ecr / Ecycle;
costDowntime   = par.Cd * Edown / Ecycle;
Cinf = costInspection + costPM + costPR + costCR + costDowntime;

share = [costInspection, costPM, costPR, costCR, costDowntime] ./ Cinf;

availability = 1 - Edown / Ecycle;

% MTTF does not enter the cost objective; skip for grid screening.
MTTF_stationary = NaN;
MTTF_fromNew = NaN;
rcondMTTF = NaN;
if doMTTF
    I = eye(Ns);
    mVec = (I - Psurv) \ cond.mttfSource;
    MTTF_stationary = sum(piVec .* mVec);
    MTTF_fromNew    = mVec(1);
    rcondMTTF = rcond(I-Psurv);
end

% Optional reliability curve from Eq. (59).
rel = struct('t', [], 'Rstationary', [], 'RfromNew', [], 'RbyState', []);
MTTFcurve = NaN;
relativeMTTFcurveError = NaN;
reliabilityTail = NaN;

if doReliabilityCurve
    % Keep the original automatic horizon rule unchanged.
    tMax = opts.tRelMaxFactor * max(MTTF_fromNew, Ecycle);
    if ~isfinite(tMax) || tMax <= 0
        tMax = 2e5;
    end
    tMax = min(tMax, 1.0e6);

    NtRel = max(2, round(opts.NtRel));
    tGrid = linspace(0, tMax, NtRel);

    Rstate = reliabilityCurveByRecursion(tGrid, states, z, wz, v, wv, pmSurvVec, ...
        RI, M, par, opts, cond.tau, PrelSurv, pmPdfMat);

    rel.t = tGrid(:);
    rel.RbyState = Rstate;
    rel.Rstationary = Rstate * piVec;
    rel.RfromNew = Rstate(:,1);

    % Independent consistency check: integral of the reliability curve
    % should approach the stationary-initial MTTF when the time horizon is long enough.
    MTTFcurve = trapz(rel.t, rel.Rstationary);
    relativeMTTFcurveError = abs(MTTFcurve - MTTF_stationary) / max(MTTF_stationary, eps);
    reliabilityTail = rel.Rstationary(end);
end

rep = struct();
rep.policy.RI = RI;
rep.policy.M  = M;
rep.costRate  = Cinf;
rep.expectedCycleLength = Ecycle;
rep.expectedCountsPerCycle = struct('inspection', 1, 'PM', Epm, ...
    'preventiveReplacement', Epr, 'correctiveReplacement', Ecr);
rep.expectedDowntimePerCycle = Edown;
rep.costRateComponents = struct( ...
    'inspection', costInspection, ...
    'PM', costPM, ...
    'preventiveReplacement', costPR, ...
    'correctiveReplacement', costCR, ...
    'downtime', costDowntime);
rep.costShare = struct( ...
    'inspection', share(1), ...
    'PM', share(2), ...
    'preventiveReplacement', share(3), ...
    'correctiveReplacement', share(4), ...
    'downtime', share(5));
rep.averageAvailability = availability;
rep.MTTF_stationary = MTTF_stationary;
rep.MTTF_fromNew = MTTF_fromNew;
rep.reliability = rel;
rep.parameters = par;
rep.stationary = struct('a0', a0, 'z', z, 'wz', wz, 'H', H, 'piVec', piVec);
rep.conditional = cond;

rep.diagnostics = struct( ...
    'maxRawTransitionMassError', max(abs(PrawMassError)), ...
    'meanRawTransitionMassError', mean(abs(PrawMassError)), ...
    'maxCorrectedTransitionMassError', max(abs(PcorrectedMassError)), ...
    'meanCorrectedTransitionMassError', mean(abs(PcorrectedMassError)), ...
    'maxDoNothingMassError', max(abs(PdoMassError)), ...
    'maxPMSuccessMassError', max(abs(PpmSuccessMassError)), ...
    'maxPMEventMassError', max(abs(PpmEventMassError)), ...
    'maxFullTransitionMassError', max(abs(PcorrectedMassError)), ...
    'meanFullTransitionMassError', mean(abs(PcorrectedMassError)), ...
    'minMTTFLinearMatrixReciprocalCondition', rcondMTTF, ...
    'MTTFfromReliabilityCurve', MTTFcurve, ...
    'relativeMTTFcurveError', relativeMTTFcurveError, ...
    'reliabilityTail', reliabilityTail);

maxRawComponentError = max([ ...
    rep.diagnostics.maxRawTransitionMassError, ...
    rep.diagnostics.maxDoNothingMassError, ...
    rep.diagnostics.maxPMSuccessMassError, ...
    rep.diagnostics.maxPMEventMassError]);

if maxRawComponentError > opts.minMassWarning && opts.warnOnMass
    warning(['Max raw transition/component mass error is %.3g. ', ...
        'Consider increasing opts.Nz/Nv/Ny.'], maxRawComponentError);
end
end
function row = buildOneStepRow(s, RI, M, par, opts, z, wz, v, wv, pmSurvVec, pmPdfMat, doMTTF)
L = par.L;
lam = par.lambda;
Nz = numel(z);

% Reliability-triggered scheduled inspection time.
tau = tauFromReliabilityThreshold(s, RI, par);
aTau = par.alpha * tau;
eNoShock = exp(-lam * tau);


%% ---------------- Scheduled reliability-triggered branch ----------------
% This branch occurs only if no shock arrives before tau.
FL_rel = gammaCdfRate(L - s, aTau, par.beta);
FM_rel = gammaCdfRate(M - s, aTau, par.beta);

pr_cr_mass   = eNoShock * max(0, 1 - FL_rel);
p_pm_rel     = eNoShock * max(0, FL_rel - FM_rel);
p_do_rel     = eNoShock * max(0, FM_rel);

fR_v = gammaPdfRate(v - s, aTau, par.beta);
fR_v(v <= s) = 0;

% Raw PM split computed on the v-grid.
pr_pr_raw = eNoShock * sum(wv .* pmSurvVec .* fR_v);
pr_success_prob_raw = eNoShock * sum(wv .* (1 - pmSurvVec) .* fR_v);
pr_pm_raw = pr_pr_raw + pr_success_prob_raw;

% Preserve the raw PR/successful-PM ratio, but force their sum to the
% theoretically computed PM-event probability p_pm_rel.
if p_pm_rel <= 0
    pr_pr = 0;
    pr_success_target = 0;
elseif pr_pm_raw > 0
    pmScale = p_pm_rel / pr_pm_raw;
    pr_pr = pr_pr_raw * pmScale;
    pr_success_target = pr_success_prob_raw * pmScale;
else
    error('Scheduled PM target is positive but the v-grid gives zero PM mass.');
end

% Raw successful-PM density and do-nothing density.
pr_success_density_raw = eNoShock * (pmPdfMat * (wv .* fR_v));

x_do = z - s;
pr_do_density_raw = eNoShock * gammaPdfRate(x_do, aTau, par.beta);
pr_do_density_raw(z <= s) = 0;

% Correct the two interior components separately. This avoids changing the
% PM/do-nothing balance when only one quadrature component has error.
[pr_success_density, pr_success_density_raw_mass] = ...
    scaleDensityToTarget(pr_success_density_raw, wz, pr_success_target);

[pr_do_density, pr_do_density_raw_mass] = ...
    scaleDensityToTarget(pr_do_density_raw, wz, p_do_rel, find(z >= s, 1, 'first'));

% Scheduled survival row used later by the reliability recursion.
PrelRow = zeros(1, Nz+1);
PrelRow(1) = pr_pr;
PrelRow(2:end) = ((pr_success_density + pr_do_density) .* wz)';

%% ----------------------- Shock-triggered branch --------------------------
qs_cr = 0;
qs_pr = 0;
qs_success_density = zeros(Nz,1);
qs_do_density = zeros(Nz,1);
p_pm_shock = 0;
p_do_shock = 0;

% Raw quantities retained only for diagnostics.
qs_pr_raw = 0;
qs_success_prob_raw = 0;
qs_success_density_raw_mass = 0;
qs_do_density_raw_mass = 0;

if tau > 0 && lam > 0
    Ny = opts.Ny;
    yGrid = ((1:Ny)' - 0.5) / Ny * tau;
    wy = tau / Ny;

    % Midpoint quadrature weights are rescaled only so that their total
    % exactly equals P(Y<tau)=1-exp(-lambda*tau). The model and grid are unchanged.
    coeffVec = lam .* exp(-lam .* yGrid) .* wy;
    shockTriggerProb = 1 - eNoShock;
    coeffSum = sum(coeffVec);
    if shockTriggerProb > 0 && coeffSum > 0
        coeffVec = coeffVec * (shockTriggerProb / coeffSum);
    end

    for iy = 1:Ny
        y = yGrid(iy);
        coeff = coeffVec(iy);

        comp = shockTransitionComponentsAtU( ...
            s, y, z, wz, v, wv, pmSurvVec, pmPdfMat, M, par, opts.Nconv);

        qs_cr = qs_cr + coeff * comp.Pcr;
        qs_pr = qs_pr + coeff * comp.Ppr;
        qs_success_density = qs_success_density + coeff * comp.pmSuccessDensity;
        qs_do_density = qs_do_density + coeff * comp.doNothingDensity;
        p_pm_shock = p_pm_shock + coeff * comp.Ppm;
        p_do_shock = p_do_shock + coeff * comp.Pdo;

        qs_pr_raw = qs_pr_raw + coeff * comp.PprRaw;
        qs_success_prob_raw = qs_success_prob_raw + coeff * comp.pmSuccessProbRaw;
        qs_success_density_raw_mass = qs_success_density_raw_mass + ...
            coeff * comp.pmSuccessDensityRawMass;
        qs_do_density_raw_mass = qs_do_density_raw_mass + ...
            coeff * comp.doNothingDensityRawMass;
    end
end

%% -------------------------- Combine branches -----------------------------
Pcr = pr_cr_mass + qs_cr;
Ppr = pr_pr + qs_pr;
Ppm = p_pm_rel + p_pm_shock;
Pdo = p_do_rel + p_do_shock;

pmSuccessDensity = pr_success_density + qs_success_density;
doNothingDensity = pr_do_density + qs_do_density;
interiorDensity = pmSuccessDensity + doNothingDensity;

mass0_full = Pcr + Ppr;

% Diagnostics BEFORE correction. These reveal the actual quadrature errors.
PprRawTotal = pr_pr_raw + qs_pr_raw;
pmSuccessProbRawTotal = pr_success_prob_raw + qs_success_prob_raw;
pmSuccessDensityRawMass = pr_success_density_raw_mass + qs_success_density_raw_mass;
doNothingDensityRawMass = pr_do_density_raw_mass + qs_do_density_raw_mass;

rawFullMassError = ...
    Pcr + PprRawTotal + pmSuccessDensityRawMass + doNothingDensityRawMass - 1;

doNothingMassError = doNothingDensityRawMass - Pdo;
pmSuccessMassError = pmSuccessDensityRawMass - pmSuccessProbRawTotal;
pmEventMassError = (PprRawTotal + pmSuccessProbRawTotal) - Ppm;

% Corrected full row. Because the PM and do-nothing components are corrected
% separately, the cost probabilities and the stationary transition kernel
% now use the same event masses.
correctedFullMassError = ...
    mass0_full + sum(interiorDensity .* wz) - 1;

% Survival transition row for MTTF: excludes corrective replacement, keeps
% preventive replacement mass at state 0 and all nonfailed interior transitions.
PsurvRow = zeros(1, Nz+1);
PsurvRow(1) = Ppr;
PsurvRow(2:end) = (interiorDensity .* wz)';

% Conditional cycle length.
if lam > 0
    mu = (1 - exp(-lam * tau)) / lam;
else
    mu = tau;
end

% Conditional downtime and source term in MTTF equation.
downtime = downtimeIntegral(s, tau, par, opts.NtInt);
if doMTTF
    mttfSource = mttfSourceIntegral(s, tau, par, opts.NtInt);
else
    mttfSource = 0;
end

row = struct();
row.tau = tau;
row.mu = mu;
row.Ppm = min(max(Ppm,0),1);
row.Ppr = min(max(Ppr,0),1);
row.Pcr = min(max(Pcr,0),1);
row.downtime = downtime;
row.mttfSource = mttfSource;
row.fullMass0 = min(max(mass0_full,0),1);
row.fullInteriorDensity = interiorDensity(:)';
row.PsurvRow = PsurvRow;
row.PrelSurvRow = PrelRow;

% New diagnostics; keep fullMassError as an alias for backward compatibility.
row.rawFullMassError = rawFullMassError;
row.correctedFullMassError = correctedFullMassError;
row.doNothingMassError = doNothingMassError;
row.pmSuccessMassError = pmSuccessMassError;
row.pmEventMassError = pmEventMassError;
row.fullMassError = correctedFullMassError;
end
function cond = initConditionalStruct(Ns)
cond = struct();
cond.tau = zeros(Ns,1);
cond.mu = zeros(Ns,1);
cond.downtime = zeros(Ns,1);
cond.Ppm = zeros(Ns,1);
cond.Ppr = zeros(Ns,1);
cond.Pcr = zeros(Ns,1);
cond.mttfSource = zeros(Ns,1);
end

%% ========================================================================
%  Reliability curve recursion, Eq. (59)
% ========================================================================
function Rstate = reliabilityCurveByRecursion(tGrid, states, z, wz, v, wv, pmSurvVec, RI, M, par, opts, tauVec, PrelSurv, pmPdfMat)
% Hybrid-policy first-failure reliability.
%
% The original code used R(t_{n-1}) whenever t-u was between t_{n-1} and
% the current time t_n. That "last-value hold" can create an artificial
% reliability plateau. Here the same recursion is kept, but the unresolved
% part of R(t-u) is represented by linear interpolation between R(t_{n-1})
% and the unknown R(t_n), yielding a small implicit linear system at each
% time step. The shock survival row also uses the same component-wise mass
% correction as the one-step kernel used for cost/MTTF calculations.

Ns = numel(states);
Nt = numel(tGrid);
if Nt < 2
    error('Reliability recursion requires at least two time grid points. Set opts.NtRel >= 2.');
end

tGrid = tGrid(:);
if any(diff(tGrid) <= 0)
    error('Reliability time grid must be strictly increasing.');
end

Rstate = zeros(Nt, Ns);
Rstate(1,:) = 1;

for it = 2:Nt
    t = tGrid(it);

    % Aimp * R(t) = rhs. Each row corresponds to one current starting state.
    Aimp = eye(Ns);
    rhs = zeros(Ns,1);

    for is = 1:Ns
        s = states(is);
        tau = tauVec(is);

        % If the horizon is before the scheduled reliability-triggered
        % inspection, the no-shock/no-failure path is still active.
        if t < tau
            rhs(is) = exp(-par.lambda * t) * ...
                gammaCdfRate(par.L - s, par.alpha * t, par.beta);
            uMax = t;
        else
            uMax = tau;
        end

        % First-shock contribution over u in (0,min(t,tau)).
        if uMax > 0 && par.lambda > 0
            Nu = max(1, round(opts.NuRel));
            uGrid = ((1:Nu)' - 0.5) / Nu * uMax;
            wu = uMax / Nu;

            for iu = 1:Nu
                u = uGrid(iu);
                Prow = shockSurvivalRowAtU( ...
                    s, u, z, wz, v, wv, pmSurvVec, pmPdfMat, M, par, opts.Nconv);

                tq = t - u;
                [knownPart, currentWeight] = reliabilityAtIntermediateTime( ...
                    tq, tGrid, Rstate, it);

                rhs(is) = rhs(is) + wu * (Prow * knownPart(:));

                if currentWeight > 0
                    Aimp(is,:) = Aimp(is,:) - wu * currentWeight * Prow;
                end
            end
        end

        % Scheduled reliability-triggered branch, active when t >= tau.
        if t >= tau
            tq = t - tau;
            [knownPart, currentWeight] = reliabilityAtIntermediateTime( ...
                tq, tGrid, Rstate, it);

            rhs(is) = rhs(is) + PrelSurv(is,:) * knownPart(:);

            if currentWeight > 0
                Aimp(is,:) = Aimp(is,:) - currentWeight * PrelSurv(is,:);
            end
        end
    end

    if rcond(Aimp) < 1e-12
        warning('Reliability recursion matrix is poorly conditioned at t = %.6g.', t);
    end

    Rnow = Aimp \ rhs;
    Rnow = min(max(real(Rnow), 0), 1);
    Rstate(it,:) = Rnow';
end
end

function [knownPart, currentWeight] = reliabilityAtIntermediateTime(tq, tGrid, Rstate, it)
% Decompose R(tq) into
%   knownPart + currentWeight * R(t_it)
% using linear interpolation.
%
% For tq <= t_{it-1}, the value is completely known from previous time
% points. For tq in (t_{it-1},t_it], the unresolved current-time value is
% retained implicitly instead of replacing it by R(t_{it-1}).

Ns = size(Rstate,2);
tNow = tGrid(it);
tPrev = tGrid(it-1);

if tq <= 0
    knownPart = ones(1,Ns);
    currentWeight = 0;
    return;
end

tol = 10 * eps(max(1, abs(tNow)));

if tq <= tPrev + tol
    if it == 2
        % The only previously known point is t=0.
        knownPart = Rstate(1,:);
    else
        tqUse = min(max(tq, tGrid(1)), tPrev);
        knownPart = interp1(tGrid(1:it-1), Rstate(1:it-1,:), tqUse, 'linear');
    end
    currentWeight = 0;
    return;
end

if tq <= tNow + tol
    dt = tNow - tPrev;
    currentWeight = min(max((tq - tPrev) / dt, 0), 1);
    knownPart = (1 - currentWeight) * Rstate(it-1,:);
    return;
end

error('Reliability recursion requested a future value beyond the current time grid point.');
end

function Prow = shockSurvivalRowAtU(s, u, z, wz, v, wv, pmSurvVec, pmPdfMat, M, par, Nconv)
% Density-in-time survival transition row at a first shock occurring at u.
% Corrective-replacement paths are excluded. The conditional state
% transition components are corrected in exactly the same way as in the
% one-step kernel used for stationary cost and MTTF calculations.

comp = shockTransitionComponentsAtU( ...
    s, u, z, wz, v, wv, pmSurvVec, pmPdfMat, M, par, Nconv);

coeff = par.lambda * exp(-par.lambda * u);

Prow = zeros(1,numel(z)+1);
Prow(1) = coeff * comp.Ppr;
Prow(2:end) = coeff * ...
    ((comp.pmSuccessDensity + comp.doNothingDensity) .* wz)';
end

function comp = shockTransitionComponentsAtU(s, u, z, wz, v, wv, pmSurvVec, pmPdfMat, M, par, Nconv)
% Conditional transition components given that the first shock occurs at u.
% This helper is shared by the one-step kernel and the reliability recursion
% so both calculations use exactly the same PM/do-nothing mass correction.

% The Beta post-maintenance density matrix is precomputed per policy.

% Partition of the pre-maintenance state at the shock inspection:
%   below M, in [M,L), or failed at/above L.
survPair = sumGammaSurvRate([par.L - s; M - s], u, par, Nconv);
survL = survPair(1);
survM = survPair(2);

Pcr = min(max(survL,0),1);
Ppm = max(survM - survL, 0);
Pdo = max(1 - survM, 0);

% Small numerical inconsistencies in the two survival evaluations are
% removed by normalizing the three mutually exclusive conditional events.
eventMass = Pcr + Ppm + Pdo;
if eventMass > 0
    Pcr = Pcr / eventMass;
    Ppm = Ppm / eventMass;
    Pdo = Pdo / eventMass;
end

% Evaluate both Gamma-sum PDFs with one shared mixture expansion.
Nv = numel(v);
guBoth = sumGammaPdfRate([v - s; z - s], u, par);
gu_v = guBoth(1:Nv);
gu_z = guBoth(Nv+1:end);
gu_v(v <= s) = 0;
gu_z(z <= s) = 0;

% Raw PR/successful-PM split on the v-grid.
PprRaw = sum(wv .* pmSurvVec .* gu_v);
pmSuccessProbRaw = sum(wv .* (1 - pmSurvVec) .* gu_v);
pmRaw = PprRaw + pmSuccessProbRaw;

if Ppm <= 0
    Ppr = 0;
    pmSuccessTarget = 0;
elseif pmRaw > 0
    pmScale = Ppm / pmRaw;
    Ppr = PprRaw * pmScale;
    pmSuccessTarget = pmSuccessProbRaw * pmScale;
else
    error('Shock-triggered PM target is positive but the v-grid gives zero PM mass.');
end

pmSuccessDensityRaw = pmPdfMat * (wv .* gu_v);
[pmSuccessDensity, pmSuccessDensityRawMass] = ...
    scaleDensityToTarget(pmSuccessDensityRaw, wz, pmSuccessTarget);

doNothingDensityRaw = gu_z;
[doNothingDensity, doNothingDensityRawMass] = ...
    scaleDensityToTarget(doNothingDensityRaw, wz, Pdo, find(z >= s, 1, 'first'));

comp = struct();
comp.Pcr = Pcr;
comp.Ppm = Ppm;
comp.Pdo = Pdo;
comp.Ppr = Ppr;
comp.pmSuccessTarget = pmSuccessTarget;
comp.pmSuccessDensity = pmSuccessDensity;
comp.doNothingDensity = doNothingDensity;

% Raw quantities for diagnostics in the integrated one-step kernel.
comp.PprRaw = PprRaw;
comp.pmSuccessProbRaw = pmSuccessProbRaw;
comp.pmSuccessDensityRawMass = pmSuccessDensityRawMass;
comp.doNothingDensityRawMass = doNothingDensityRawMass;
end

function [density, rawMass] = scaleDensityToTarget(rawDensity, weights, targetMass, fallbackIndex)
% Scale one transition-density component to its own theoretical mass.
% PM-success and do-nothing densities are corrected separately.
%
% For a do-nothing transition from the highest state-grid midpoint, the
% support (s,M) may lie entirely inside the same terminal state cell, so a
% midpoint density evaluation can have zero raw mass. In that special case
% fallbackIndex assigns the target probability to the cell containing s.

if nargin < 4
    fallbackIndex = [];
end

rawDensity = max(real(rawDensity),0);
rawMass = sum(rawDensity(:) .* weights(:));

if targetMass <= 0
    density = zeros(size(rawDensity));
    return;
end

if rawMass > 0
    density = rawDensity * (targetMass / rawMass);
    return;
end

if ~isempty(fallbackIndex)
    density = zeros(size(rawDensity));
    density(fallbackIndex) = targetMass / weights(fallbackIndex);
    return;
end

error('Positive transition target mass encountered with zero raw density mass.');
end
function tab = gridSearchPolicy(RIGrid, MGrid, par, opts, label)
% Keep numerical diagnostics, but suppress repetitive per-policy warnings
% during screening; the final optimum is fully evaluated with warnings on.
optsGrid = opts;
optsGrid.warnOnMass = false;
[RIList, MList] = ndgrid(RIGrid(:), MGrid(:));
RIList = RIList(:);
MList  = MList(:);
n = numel(RIList);

CostRate = nan(n,1);
Availability = nan(n,1);
Ecycle = nan(n,1);
MTTFnew = nan(n,1);
MTTFstat = nan(n,1);
Ppm = nan(n,1);
Ppr = nan(n,1);
Pcr = nan(n,1);
Down = nan(n,1);
Err = strings(n,1);

usePar = opts.useParallel && license('test','Distrib_Computing_Toolbox');
if usePar
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('Starting parallel pool...\n');
        pool = parpool;
    end
    fprintf('Parallel pool active with %d workers.\n', pool.NumWorkers);
    dq = parallel.pool.DataQueue;
    progress = makeProgressPrinter(n, label, opts.progressEverySeconds);
    afterEach(dq, @(~) progress());

    parfor k = 1:n
        try
            rep = evaluatePolicy(RIList(k), MList(k), par, optsGrid, false, false);
            CostRate(k) = rep.costRate;
            Availability(k) = rep.averageAvailability;
            Ecycle(k) = rep.expectedCycleLength;
            MTTFnew(k) = rep.MTTF_fromNew;
            MTTFstat(k) = rep.MTTF_stationary;
            Ppm(k) = rep.expectedCountsPerCycle.PM;
            Ppr(k) = rep.expectedCountsPerCycle.preventiveReplacement;
            Pcr(k) = rep.expectedCountsPerCycle.correctiveReplacement;
            Down(k) = rep.expectedDowntimePerCycle;
        catch ME
            Err(k) = string(ME.message);
        end
        send(dq, k);
    end
else
    fprintf('Running serial grid search.\n');
    t0 = tic;
    lastPrint = 0;
    for k = 1:n
        try
            rep = evaluatePolicy(RIList(k), MList(k), par, optsGrid, false, false);
            CostRate(k) = rep.costRate;
            Availability(k) = rep.averageAvailability;
            Ecycle(k) = rep.expectedCycleLength;
            MTTFnew(k) = rep.MTTF_fromNew;
            MTTFstat(k) = rep.MTTF_stationary;
            Ppm(k) = rep.expectedCountsPerCycle.PM;
            Ppr(k) = rep.expectedCountsPerCycle.preventiveReplacement;
            Pcr(k) = rep.expectedCountsPerCycle.correctiveReplacement;
            Down(k) = rep.expectedDowntimePerCycle;
        catch ME
            Err(k) = string(ME.message);
        end
        if toc(t0)-lastPrint > opts.progressEverySeconds || k == n
            lastPrint = toc(t0);
            fprintf('[%s] %4d/%4d policies evaluated (%.1f%%), elapsed %.1fs\n', ...
                label, k, n, 100*k/n, toc(t0));
        end
    end
end

tab = table(RIList, MList, CostRate, Availability, Ecycle, MTTFnew, MTTFstat, ...
    Ppm, Ppr, Pcr, Down, Err, ...
    'VariableNames', {'RI','M','CostRate','Availability','ExpectedCycleLength', ...
    'MTTF_fromNew','MTTF_stationary','ExpectedPM','ExpectedPR','ExpectedCR', ...
    'ExpectedDowntime','Error'});
end

function progress = makeProgressPrinter(n, label, printEverySeconds)
count = 0;
t0 = tic;
lastPrintTime = 0;
progress = @update;
    function update()
        count = count + 1;
        elapsed = toc(t0);
        if elapsed - lastPrintTime >= printEverySeconds || count == n
            lastPrintTime = elapsed;
            rate = count / max(elapsed, eps);
            remain = (n - count) / max(rate, eps);
            fprintf('[%s] %4d/%4d policies evaluated (%.1f%%), elapsed %.1fs, ETA %.1fs\n', ...
                label, count, n, 100*count/n, elapsed, remain);
        end
    end
end

%% ========================================================================
%  Numerical building blocks
% ========================================================================
function tau = tauFromReliabilityThreshold(s, RI, par)
% Solve P(X_tau < L-s)=RI, where X_t~Gamma(alpha*t,beta).
margin = par.L - s;
if margin <= 0
    tau = 0;
    return;
end
fun = @(t) gammaCdfRate(margin, par.alpha * t, par.beta) - RI;
lo = 0;
hi = max(1, margin * par.beta / max(par.alpha, eps));
% Increase hi until reliability falls below RI.
while fun(hi) > 0
    hi = hi * 2;
    if hi > 1e10
        error('Cannot bracket tau. Check RI and degradation parameters.');
    end
end
% fzero uses safeguarded interpolation with the existing sign bracket;
% avoids always evaluating gammainc 80 times for each initial state.
if fun(hi) == 0
    tau = hi;
else
    tau = fzero(fun, [lo, hi]);
end
end

function d = downtimeIntegral(s, tau, par, Nt)
if tau <= 0
    d = 0;
    return;
end
t = ((1:Nt)' - 0.5) / Nt * tau;
wt = tau / Nt;
shape = par.alpha .* t;
survFail = gammaSurvRate(par.L - s, shape, par.beta);
d = sum(wt .* exp(-par.lambda .* t) .* survFail);
end

function b = mttfSourceIntegral(s, tau, par, Nt)
if tau <= 0
    b = 0;
    return;
end
t = ((1:Nt)' - 0.5) / Nt * tau;
wt = tau / Nt;
shape = par.alpha .* t;
R = gammaCdfRate(par.L - s, shape, par.beta); % P(X_t < L-s)
b = sum(wt .* exp(-par.lambda .* t) .* R);
end

function pdf = gammaPdfRate(x, a, beta)
% Gamma PDF with rate beta. Supports scalar-vector broadcasting.
[xb, ab] = commonSize(x, a);
pdf = zeros(size(xb));
idx = xb > 0 & ab > 0;
if any(idx(:))
    aa = ab(idx);
    xx = xb(idx);
    logp = aa .* log(beta) - gammaln(aa) + (aa - 1).*log(xx) - beta.*xx;
    pdf(idx) = exp(logp);
end
% Degenerate case a=0 is not represented by a density and is intentionally 0.
pdf(~isfinite(pdf)) = 0;
end

function cdf = gammaCdfRate(x, a, beta)
% Gamma CDF P(X <= x) with rate beta. Supports scalar-vector broadcasting.
[xb, ab] = commonSize(x, a);
cdf = zeros(size(xb));
idxDeg = xb > 0 & (ab <= 0);
cdf(idxDeg) = 1; % X_0 = 0, so P(X_0 < x)=1 for x>0
idx = xb > 0 & (ab > 0);
if any(idx(:))
    cdf(idx) = gammainc(beta .* xb(idx), ab(idx), 'lower');
end
cdf(isinf(xb) & xb > 0) = 1;
cdf = min(max(cdf,0),1);
end

function surv = gammaSurvRate(x, a, beta)
% Gamma survival P(X >= x) with rate beta. Supports scalar-vector broadcasting.
[xb, ab] = commonSize(x, a);
surv = zeros(size(xb));
surv(xb <= 0) = 1;
idxDeg = xb > 0 & (ab <= 0);
surv(idxDeg) = 0; % X_0=0, so P(X_0>=x)=0 for x>0
idx = xb > 0 & (ab > 0);
if any(idx(:))
    surv(idx) = gammainc(beta .* xb(idx), ab(idx), 'upper');
end
surv = min(max(surv,0),1);
end

function [xb, ab] = commonSize(x, a)
% Minimal scalar expansion helper for MATLAB releases without full implicit indexing support.
if isscalar(x) && isscalar(a)
    xb = x;
    ab = a;
elseif isscalar(x)
    ab = a;
    xb = ones(size(a)) .* x;
elseif isscalar(a)
    xb = x;
    ab = ones(size(x)) .* a;
else
    if ~isequal(size(x), size(a))
        error('Inputs x and a must be scalar or have the same size.');
    end
    xb = x;
    ab = a;
end
end

function pdf = sumGammaPdfRate(u, y, par)
% Exact nonnegative Gamma-mixture expansion for the PDF of X_y + W.
% X_y~Gamma(alpha*y, rate beta), W~Gamma(alphaS, rate betaS).
% Uses a negative-binomial mixture of Gamma(shapeBase+k, rateHigh),
% with a rigorous 1e-12 bound on the omitted mixture probability.
pdf = zeros(size(u));
idx = u > 0;
if ~any(idx(:))
    return;
end
uPos = u(idx);
[rateHigh, shapeBase, shapeLow, prob] = gammaSumMixtureParams(y, par);
shape = shapeBase;
weights = exp(shapeLow * log(prob));
weightTotal = weights;
component = gammaPdfRate(uPos, shape, rateHigh);
result = weights .* component;
for k = 1:2000
    if weightTotal >= 1-1e-12
        break;
    end
    weights = weights .* ((shapeLow + k-1) / k) .* (1-prob);
    component = component .* (rateHigh .* uPos ./ (shapeBase + k-1));
    result = result + weights .* component;
    weightTotal = weightTotal + weights;
end
if weightTotal < 1-1e-9
    error('Gamma PDF mixture did not converge: check shock/degradation parameters.');
end
pdf(idx) = max(result, 0);
end

function surv = sumGammaSurvRate(u, y, par, Nconv) %#ok<INUSD>
% Exact nonnegative Gamma-mixture expansion for P(X_y+W >= u).
% Accepts scalar or vector u so both threshold tails can be found together.
surv = ones(size(u));
idx = u > 0;
if ~any(idx(:))
    return;
end
uPos = u(idx);
[rateHigh, shapeBase, shapeLow, prob] = gammaSumMixtureParams(y, par);
weights = exp(shapeLow * log(prob));
weightTotal = weights;
componentSurv = gammaSurvRate(uPos, shapeBase, rateHigh);
componentPDF = gammaPdfRate(uPos, shapeBase, rateHigh);
result = weights .* componentSurv;
for k = 1:2000
    if weightTotal >= 1-1e-12
        break;
    end
    weights = weights .* ((shapeLow + k-1) / k) .* (1-prob);
    componentSurv = min(1, componentSurv + (uPos ./ (shapeBase+k-1)) .* componentPDF);
    componentPDF = componentPDF .* (rateHigh .* uPos ./ (shapeBase+k-1));
    result = result + weights .* componentSurv;
    weightTotal = weightTotal + weights;
end
if weightTotal < 1-1e-9
    error('Gamma survival mixture did not converge: check parameters.');
end
surv(idx) = min(max(result,0),1);
end

function [rateHigh, shapeBase, shapeLow, prob] = gammaSumMixtureParams(y, par)
% If rates are unequal, rewrite the lower-rate Gamma random variable
% as a negative-binomial mixture of Gammas at the larger rate.
a = par.alpha*y;
c = par.alphaS;
b = par.beta;
d = par.betaS;
shapeBase = a+c;
rateHigh = max(b,d);
prob = min(b,d)/rateHigh;
if b <= d
    shapeLow = a;
else
    shapeLow = c;
end
end

function Bsurv = betaSurvRate(x, p, q)
% P(Beta(p,q) >= x).
Bsurv = zeros(size(x));
Bsurv(x <= 0) = 1;
Bsurv(x >= 1) = 0;
idx = x > 0 & x < 1;
if any(idx(:))
    Bsurv(idx) = 1 - betainc(x(idx), p, q);
end
Bsurv = min(max(Bsurv,0),1);
end

function pdf = betaPdfUnit(x, p, q)
pdf = zeros(size(x));
idx = x > 0 & x < 1;
if any(idx(:))
    xx = x(idx);
    pdf(idx) = exp((p-1).*log(xx) + (q-1).*log1p(-xx) - betaln(p,q));
end
pdf(~isfinite(pdf)) = 0;
end

function mat = pmPdfMatrix(z, v, p, q)
% mat(i,j) = f_{S*}(z_i | v_j), where S*=v_j*B, B~Beta(p,q).
z = z(:);
v = v(:)';
r = z ./ v;
mat = betaPdfUnit(r, p, q) ./ v;
mat(r <= 0 | r >= 1) = 0;
mat(~isfinite(mat)) = 0;
end

%% ========================================================================
%  Reporting and plotting
% ========================================================================

function validatePolicy(RI, M, par)
if ~(RI > 0 && RI < 1)
    error('RI must satisfy 0 < RI < 1.');
end
if ~(M > 0 && M < par.L)
    error('M must satisfy 0 < M < L.');
end
end
