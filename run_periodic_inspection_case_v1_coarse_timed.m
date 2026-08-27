%RUN_PERIODIC_INSPECTION_CASE_V1_COARSE_TIMED
% Script version for an equal-interval inspection-CBM policy.
%
% Policy definition:
%   - Inspections are performed every fixed interval T.
%   - Shock arrivals do NOT trigger inspections.
%   - At an inspection:
%       * if degradation < M: do nothing;
%       * if M <= degradation < L: perform imperfect PM;
%       * if PM result remains >= M: preventive replacement;
%       * if degradation >= L: corrective replacement.
%
% The degradation increment over any interval of length t is
%
%   Delta S_t = X_t + sum_{j=1}^{N(t)} W_j.
%
% Decision variables:
%   T  - equal inspection interval
%   M  - preventive-maintenance threshold
%
% The physical and cost parameters below are kept identical to the
% currently supplied hybrid-policy code. Only the policy is changed.
%
% How to run:
%   >> run_periodic_inspection_case_v1_coarse_timed
%
% Main output variables:
%   results       - complete result structure
%   coarseTab     - coarse-grid search table
%   repBest       - evaluation at the coarse-grid optimum
%   timingLog     - runtime table
%
% Reliability:
%   The final evaluation below is set to false by default, matching the
%   current optimization-only workflow. Change the final argument of
%   evaluatePeriodicPolicy from false to true to generate reliability.

clearvars;
clc;
close all;

scriptTimer = tic;
timingNames = {};
timingSeconds = [];

fprintf('\n============================================================\n');
fprintf(' Equal-interval inspection policy -- V1 coarse timed\n');
fprintf('============================================================\n');

%% --------------------- Model parameters ---------------------
par.alpha  = 1.02;       % Gamma-process shape rate alpha
par.beta   = 1.20;       % Gamma-process rate beta
par.lambda = 0.25;        % HPP shock arrival intensity

% Shock magnitude W ~ Gamma(alphaS,betaS), betaS is RATE.
par.alphaS = 2.0;
par.betaS  = 1.5;

% Imperfect PM ratio B ~ Beta(p,q), post-PM degradation S*=B*v.
par.p = 6;
par.q = 5;

% Failure threshold.
par.L = 12.5;

%% --------------------- Cost parameters ----------------------
par.Ci = 10;     % inspection cost
par.Cm = 100;    % imperfect PM cost
par.Cp = 280;    % preventive replacement cost
par.Cc = 350;    % corrective replacement cost
par.Cd = 50;     % downtime cost rate

%% --------------------- Numerical options --------------------
opts.Nz       = 60;     % grid points on post-maintenance interval (0,M)
opts.Nv       = 80;     % quadrature points for v in [M,L]
opts.NtInt    = 80;     % time quadrature points for downtime/MTTF source

% Reliability recursion. It is only used when doReliabilityCurve=true.
opts.NtRel          = 301;
opts.maxNtRel       = 4000;
opts.dtRelMaximum   = 0.20;
opts.tRelMaxFactor  = 6.0;
opts.reliabilityTailWarning = 1e-5;

% Compound-degradation table.
opts.NtCompound   = 1201;
opts.NuCompound   = 1201;
opts.NfftCompound = 8192;

opts.useParallel = true;
opts.showProgress = true;
opts.progressEverySeconds = 3;

% Warning threshold for RAW, pre-correction quadrature errors.
opts.minMassWarning = 1e-3;

%% --------------------- Coarse grid search -------------------
% The M-grid is kept the same as in the supplied hybrid-policy code.
coarse.MGrid = linspace(0.75,0.99,25) * par.L;

% T is the new decision variable. The range below is deliberately broad
% for the current parameter set. If the optimum is on either boundary,
% enlarge the range before reporting the final result.
coarse.TGrid = linspace(0.25,8.00,25);

% The compound table only needs to cover the largest candidate interval.
opts.tCompoundMax = max(coarse.TGrid);

%% --------------------- Precompute compound distribution -----
fprintf('\n[0] Precomputing compound-degradation PDF and CDF table...\n');

functionTimer = tic;
compound = buildCompoundDistributionTable(par,opts);
timingNames{end+1,1} = 'buildCompoundDistributionTable';
timingSeconds(end+1,1) = toc(functionTimer);

fprintf('Compound table completed. t in [0, %.4g], u in [0, %.4g].\n', ...
    compound.tGrid(end),compound.uGrid(end));
fprintf('Minimum retained probability at FFT upper bound = %.8f\n', ...
    min(compound.totalRetainedMass));
fprintf('Maximum probability in the upper 2%% of the FFT grid = %.3e\n', ...
    max(compound.upperTailMass));

%% --------------------- Coarse-grid optimization -------------
fprintf('\n[1] Periodic-policy coarse grid search: %d x %d = %d policies\n', ...
    numel(coarse.TGrid),numel(coarse.MGrid), ...
    numel(coarse.TGrid)*numel(coarse.MGrid));

functionTimer = tic;
coarseTab = gridSearchPeriodicPolicy( ...
    coarse.TGrid,coarse.MGrid,par,opts,compound,'periodic coarse');
timingNames{end+1,1} = 'gridSearchPeriodicPolicy';
timingSeconds(end+1,1) = toc(functionTimer);

validCost = isfinite(coarseTab.CostRate);

if ~any(validCost)
    error('No valid periodic-policy candidate was obtained. Inspect coarseTab.Error.');
end

validIndices = find(validCost);
[~,localBest] = min(coarseTab.CostRate(validCost));
idxBest = validIndices(localBest);

bestPolicy.T = coarseTab.T(idxBest);
bestPolicy.M = coarseTab.M(idxBest);
bestPolicy.CostRate = coarseTab.CostRate(idxBest);

fprintf('\n[2] Re-evaluating the periodic-policy optimum without reliability curve...\n');

functionTimer = tic;
repBest = evaluatePeriodicPolicy( ...
    bestPolicy.T,bestPolicy.M,par,opts,compound,false);
% Change the final false to true when the reliability curve is needed.
timingNames{end+1,1} = 'evaluatePeriodicPolicy at coarse optimum';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
printPeriodicPolicyReport( ...
    repBest,'Final equal-interval policy from coarse-grid optimization');
timingNames{end+1,1} = 'printPeriodicPolicyReport';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
checkBoundaryPeriodicPolicy( ...
    bestPolicy.T,coarse.TGrid,bestPolicy.M,coarse.MGrid,'periodic coarse');
timingNames{end+1,1} = 'checkBoundaryPeriodicPolicy';
timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Plots --------------------------------
functionTimer = tic;
try
    plotPeriodicCostContour(coarseTab,par,repBest);
catch ME
    warning('Periodic-policy cost-contour plotting failed: %s',ME.message);
end
timingNames{end+1,1} = 'plotPeriodicCostContour';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
try
    plotPeriodicReliabilityCurve(repBest);
catch ME
    warning('Periodic-policy reliability plotting failed: %s',ME.message);
end
timingNames{end+1,1} = 'plotPeriodicReliabilityCurve';
timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Results and timing --------------------
timingNames{end+1,1} = 'total script time before save';
timingSeconds(end+1,1) = toc(scriptTimer);

timingLog = table( ...
    timingNames,timingSeconds, ...
    'VariableNames',{'FunctionOrStep','Seconds'});

fprintf('\n---------------- Runtime log for main-script function calls ----------------\n');
disp(timingLog);

results = struct();
results.policyName = 'Equal-interval inspection';
results.parameters = par;
results.options = opts;
results.coarseGrid = coarse;
results.coarseGridTable = coarseTab;
results.bestPolicy = bestPolicy;
results.finalBest = repBest;
results.timingLog = timingLog;

results.compoundDiagnostics = struct( ...
    'tMaximum',compound.tGrid(end), ...
    'uMaximumForPolicyCalculations',compound.uGrid(end), ...
    'fftUpperBound',compound.fftUpperBound, ...
    'minimumRetainedMass',min(compound.totalRetainedMass), ...
    'maximumUpperTailMass',max(compound.upperTailMass));

save('periodic_inspection_results_v1_coarse_timed.mat','results');
fprintf('\nSaved results to periodic_inspection_results_v1_coarse_timed.mat\n');
fprintf('Done.\n');

%% ========================================================================
%  Main policy evaluator
% ========================================================================
function rep = evaluatePeriodicPolicy( ...
    T,M,par,opts,compound,doReliabilityCurve)

if nargin < 6
    doReliabilityCurve = false;
end

validatePeriodicPolicy(T,M,par,compound);

% Midpoint grid on (0,M), representing the continuous post-maintenance state.
Nz = opts.Nz;
z = ((1:Nz)' - 0.5) / Nz * M;
wz = ones(Nz,1) * (M/Nz);
states = [0;z];
Ns = numel(states);

% Midpoint grid on [M,L] for maintenance-state integrals.
Nv = opts.Nv;
v = M + ((1:Nv)' - 0.5) / Nv * (par.L-M);
wv = ones(Nv,1) * ((par.L-M)/Nv);

pmSurvVec = betaSurvRate(M./v,par.p,par.q);
pmPdfMat = pmPdfMatrix(z,v,par.p,par.q);

% Full transition kernel and first-failure survival kernel.
Kdens = zeros(Ns,Nz);
Q1full = zeros(Ns,1);
Psurv = zeros(Ns,Ns);

% RAW and corrected probability-mass diagnostics.
PrawMassError = zeros(Ns,1);
PcorrectedMassError = zeros(Ns,1);
PdoNothingMassError = zeros(Ns,1);
PpmSuccessMassError = zeros(Ns,1);
PpmEventMassError = zeros(Ns,1);

cond = initConditionalStructPeriodic(Ns);

for i = 1:Ns
    row = buildPeriodicOneStepRow( ...
        states(i),T,M,par,opts,compound, ...
        z,wz,v,wv,pmSurvVec,pmPdfMat);

    Kdens(i,:) = row.fullInteriorDensity;
    Q1full(i) = row.fullMass0;
    Psurv(i,:) = row.PsurvRow;

    PrawMassError(i) = row.rawMassError;
    PcorrectedMassError(i) = row.correctedMassError;
    PdoNothingMassError(i) = row.doNothingMassError;
    PpmSuccessMassError(i) = row.pmSuccessMassError;
    PpmEventMassError(i) = row.pmEventMassError;

    cond.mu(i) = row.mu;
    cond.downtime(i) = row.downtime;
    cond.Ppm(i) = row.Ppm;
    cond.Ppr(i) = row.Ppr;
    cond.Pcr(i) = row.Pcr;
    cond.mttfSource(i) = row.mttfSource;
end

% Stationary distribution:
%   pi(ds)=a0*delta_0(ds)+a0*H(s)ds.
A = zeros(Nz,Nz);

for iz = 1:Nz
    for js = 1:Nz
        A(iz,js) = Kdens(js+1,iz)*wz(js);
    end
end

b = Kdens(1,:)';
H = (eye(Nz)-A)\b;

H(H < 0 & H > -1e-10) = 0;

if any(H < -1e-8)
    warning('Stationary density H contains noticeable negative values. Increase grid sizes.');
end

H = max(H,0);

a0 = 1/(1+sum(H.*wz));
piVec = [a0; a0*H.*wz];
piVec = piVec/sum(piVec);

% For an equal-interval policy, the cycle length is exactly T for every s.
Ecycle = sum(piVec.*cond.mu);

Epm = sum(piVec.*cond.Ppm);
Epr = sum(piVec.*cond.Ppr);
Ecr = sum(piVec.*cond.Pcr);
Edown = sum(piVec.*cond.downtime);

costInspection = par.Ci/Ecycle;
costPM = par.Cm*Epm/Ecycle;
costPR = par.Cp*Epr/Ecycle;
costCR = par.Cc*Ecr/Ecycle;
costDowntime = par.Cd*Edown/Ecycle;

Cinf = costInspection+costPM+costPR+costCR+costDowntime;

share = [ ...
    costInspection,costPM,costPR,costCR,costDowntime]/Cinf;

availability = 1-Edown/Ecycle;

% MTTF:
%   m = b + Psurv*m, hence (I-Psurv)m=b.
I = eye(Ns);
mVec = (I-Psurv)\cond.mttfSource;

MTTF_stationary = sum(piVec.*mVec);
MTTF_fromNew = mVec(1);

% Optional first-failure reliability curve.
rel = struct('t',[],'Rstationary',[],'RfromNew',[],'RbyState',[]);

if doReliabilityCurve
    tMax = opts.tRelMaxFactor* ...
        max([MTTF_stationary,MTTF_fromNew,Ecycle]);

    if ~isfinite(tMax) || tMax <= 0
        tMax = 1e3;
    end

    tMax = min(tMax,1e6);

    % Make T an exact integer multiple of the reliability time step. This
    % removes interpolation error from the fixed-delay recursion.
    nCycles = max(1,ceil(tMax/T));
    nSub = max(1,ceil(T/max(opts.dtRelMaximum,eps)));

    if nCycles*nSub+1 > opts.maxNtRel
        nSub = max(1,floor((opts.maxNtRel-1)/nCycles));

        warning(['Periodic reliability grid was coarsened to satisfy ', ...
            'opts.maxNtRel=%d. Increase opts.maxNtRel for higher accuracy.'], ...
            opts.maxNtRel);
    end

    % Also respect the requested minimum number of points when possible.
    nSubForMinimum = ceil((max(round(opts.NtRel),2)-1)/nCycles);
    nSub = max(nSub,nSubForMinimum);

    if nCycles*nSub+1 > opts.maxNtRel
        nSub = max(1,floor((opts.maxNtRel-1)/nCycles));
    end

    dt = T/nSub;
    tGrid = (0:nCycles*nSub)'*dt;

    Rstate = reliabilityCurvePeriodic( ...
        tGrid,states,T,Psurv,par,compound);

    rel.t = tGrid;
    rel.RbyState = Rstate;
    rel.Rstationary = enforceReliabilityVector(Rstate*piVec);
    rel.RfromNew = enforceReliabilityVector(Rstate(:,1));
end

rep = struct();

rep.policy.T = T;
rep.policy.inspectionInterval = T;
rep.policy.M = M;

rep.costRate = Cinf;
rep.expectedCycleLength = Ecycle;

rep.expectedCountsPerCycle = struct( ...
    'inspection',1, ...
    'PM',Epm, ...
    'preventiveReplacement',Epr, ...
    'correctiveReplacement',Ecr);

rep.expectedDowntimePerCycle = Edown;

rep.costRateComponents = struct( ...
    'inspection',costInspection, ...
    'PM',costPM, ...
    'preventiveReplacement',costPR, ...
    'correctiveReplacement',costCR, ...
    'downtime',costDowntime);

rep.costShare = struct( ...
    'inspection',share(1), ...
    'PM',share(2), ...
    'preventiveReplacement',share(3), ...
    'correctiveReplacement',share(4), ...
    'downtime',share(5));

rep.averageAvailability = availability;
rep.MTTF_stationary = MTTF_stationary;
rep.MTTF_fromNew = MTTF_fromNew;

rep.reliability = rel;
rep.parameters = par;

rep.stationary = struct( ...
    'a0',a0, ...
    'z',z, ...
    'wz',wz, ...
    'H',H, ...
    'piVec',piVec);

rep.conditional = cond;

maximumUncorrectedComponentError = max([ ...
    max(abs(PrawMassError)), ...
    max(abs(PdoNothingMassError)), ...
    max(abs(PpmSuccessMassError)), ...
    max(abs(PpmEventMassError))]);

rep.diagnostics = struct( ...
    'maxRawTransitionMassError',max(abs(PrawMassError)), ...
    'meanRawTransitionMassError',mean(abs(PrawMassError)), ...
    'maxCorrectedTransitionMassError',max(abs(PcorrectedMassError)), ...
    'maxDoNothingMassError',max(abs(PdoNothingMassError)), ...
    'maxPMSuccessMassError',max(abs(PpmSuccessMassError)), ...
    'maxPMEventMassError',max(abs(PpmEventMassError)), ...
    'maxUncorrectedComponentError',maximumUncorrectedComponentError, ...
    'minMTTFLinearMatrixReciprocalCondition',rcond(I-Psurv), ...
    'stationaryTransitionMassCheck',sum(piVec.*(Q1full+sum(Kdens.*wz',2))));

if doReliabilityCurve && ~isempty(rel.t)
    rep.diagnostics.MTTFfromReliabilityCurve = ...
        trapz(rel.t,rel.Rstationary);
    rep.diagnostics.reliabilityTail = rel.Rstationary(end);
    rep.diagnostics.relativeMTTFcurveError = ...
        abs(rep.diagnostics.MTTFfromReliabilityCurve-MTTF_stationary) / ...
        max(MTTF_stationary,eps);
else
    rep.diagnostics.MTTFfromReliabilityCurve = NaN;
    rep.diagnostics.reliabilityTail = NaN;
    rep.diagnostics.relativeMTTFcurveError = NaN;
end

if rep.diagnostics.maxUncorrectedComponentError > opts.minMassWarning
    warning(['Periodic-policy maximum uncorrected component error is %.3g. ', ...
        'Increase opts.Nz, opts.Nv, or the compound grids.'], ...
        rep.diagnostics.maxUncorrectedComponentError);
end

if doReliabilityCurve && ...
        rep.diagnostics.reliabilityTail > opts.reliabilityTailWarning

    warning(['Periodic-policy reliability at the final time is %.3g. ', ...
        'Increase opts.tRelMaxFactor or opts.maxNtRel.'], ...
        rep.diagnostics.reliabilityTail);
end
end

%% ========================================================================
%  One-step transition kernel for the fixed inspection interval T
% ========================================================================
function row = buildPeriodicOneStepRow( ...
    s,T,M,par,opts,compound,z,wz,v,wv,pmSurvVec,pmPdfMat)

% Inspection-state degradation:
%   V = s + Delta S_T,
%   Delta S_T = X_T + sum_{j=1}^{N(T)} W_j.
g_v = compoundPdfRate(v-s,T,compound);
g_v(v <= s) = 0;

g_z = compoundPdfRate(z-s,T,compound);
g_z(z <= s) = 0;

F_M = compoundCdfRate(M-s,T,compound);
F_L = compoundCdfRate(par.L-s,T,compound);

F_M = min(max(F_M,0),1);
F_L = min(max(F_L,F_M),1);

% Inspection outcomes:
%   V<M      : do nothing
%   M<=V<L   : imperfect PM
%   V>=L     : corrective replacement
Pcr = 1-F_L;
Ppm = F_L-F_M;

% Split the PM event into:
%   successful PM -> post-PM state in (0,M)
%   preventive replacement -> state 0
PprRaw = sum(wv.*pmSurvVec.*g_v);

pmSuccessProbabilityRaw = ...
    sum(wv.*(1-pmSurvVec).*g_v);

pmEventMassRaw = PprRaw+pmSuccessProbabilityRaw;

if Ppm <= 0
    Ppr = 0;
    pmSuccessTarget = 0;
elseif pmEventMassRaw > 0
    pmEventScale = Ppm/pmEventMassRaw;
    Ppr = max(PprRaw,0)*pmEventScale;
    pmSuccessTarget = max(pmSuccessProbabilityRaw,0)*pmEventScale;
else
    Ppr = Ppm;
    pmSuccessTarget = 0;
end

Ppr = min(max(Ppr,0),Ppm);
pmSuccessTarget = max(Ppm-Ppr,0);

% Raw continuous transition densities.
pmSuccessDensityRaw = pmPdfMat*(wv.*g_v);
doNothingDensityRaw = g_z;

% Theoretical masses of the two interior components.
doNothingTarget = F_M;

doNothingMassRaw = sum(doNothingDensityRaw.*wz);
pmSuccessMassRaw = sum(pmSuccessDensityRaw.*wz);

doNothingMassError = doNothingMassRaw-doNothingTarget;
pmSuccessMassError = ...
    pmSuccessMassRaw-pmSuccessProbabilityRaw;
pmEventMassError = pmEventMassRaw-Ppm;

rawMassError = ...
    Pcr+PprRaw+doNothingMassRaw+pmSuccessMassRaw-1;

% Correct do-nothing and successful-PM components separately.
doNothingDensity = doNothingDensityRaw;

if doNothingTarget <= 0
    doNothingDensity(:) = 0;
elseif doNothingMassRaw > 0
    doNothingDensity = doNothingDensity*(doNothingTarget/doNothingMassRaw);
else
    index = find(z > s,1,'first');

    if isempty(index)
        index = numel(z);
    end

    doNothingDensity(:) = 0;
    doNothingDensity(index) = doNothingTarget/wz(index);
end

pmSuccessDensity = pmSuccessDensityRaw;

if pmSuccessTarget <= 0
    pmSuccessDensity(:) = 0;
elseif pmSuccessMassRaw > 0
    pmSuccessDensity = ...
        pmSuccessDensity*(pmSuccessTarget/pmSuccessMassRaw);
else
    fallbackDensity = pmPdfMat*wv;
    fallbackMass = sum(fallbackDensity.*wz);

    if fallbackMass > 0
        pmSuccessDensity = ...
            fallbackDensity*(pmSuccessTarget/fallbackMass);
    else
        pmSuccessDensity(:) = 0;
        pmSuccessDensity(1) = pmSuccessTarget/wz(1);
    end
end

interiorDensity = doNothingDensity+pmSuccessDensity;

mass0Full = Pcr+Ppr;

correctedMassError = ...
    mass0Full+sum(interiorDensity.*wz)-1;

% First-failure survival kernel:
% corrective replacement is excluded because first failure has occurred;
% preventive replacement is retained because it occurs before failure.
PsurvRow = zeros(1,numel(z)+1);
PsurvRow(1) = Ppr;
PsurvRow(2:end) = (interiorDensity.*wz)';

% Fixed cycle length.
mu = T;

% b(s)=int_0^T P(T_f>t|s)dt
%     =int_0^T F_Delta(L-s;t)dt.
mttfSource = compoundReliabilityIntegral( ...
    par.L-s,T,compound,opts.NtInt);

% d(s)=int_0^T P(T_f<=t|s)dt = T-b(s).
downtime = max(T-mttfSource,0);

row = struct();
row.mu = mu;
row.Ppm = Ppm;
row.Ppr = Ppr;
row.Pcr = Pcr;
row.downtime = downtime;
row.mttfSource = mttfSource;

row.fullMass0 = mass0Full;
row.fullInteriorDensity = interiorDensity(:)';
row.PsurvRow = PsurvRow;

row.rawMassError = rawMassError;
row.correctedMassError = correctedMassError;
row.doNothingMassError = doNothingMassError;
row.pmSuccessMassError = pmSuccessMassError;
row.pmEventMassError = pmEventMassError;
end

%% ========================================================================
%  Conditional structure
% ========================================================================
function cond = initConditionalStructPeriodic(Ns)
cond = struct();
cond.mu = zeros(Ns,1);
cond.downtime = zeros(Ns,1);
cond.Ppm = zeros(Ns,1);
cond.Ppr = zeros(Ns,1);
cond.Pcr = zeros(Ns,1);
cond.mttfSource = zeros(Ns,1);
end

%% ========================================================================
%  Reliability recursion for fixed interval T
% ========================================================================
function Rstate = reliabilityCurvePeriodic( ...
    tGrid,states,T,Psurv,par,compound)

Ns = numel(states);
Nt = numel(tGrid);

if Nt < 2
    error('Periodic reliability recursion requires at least two time points.');
end

dt = tGrid(2)-tGrid(1);
nDelay = round(T/dt);

if nDelay < 1 || abs(nDelay*dt-T) > 1e-8*max(1,T)
    error('Reliability time grid must make T an integer multiple of dt.');
end

Rstate = zeros(Nt,Ns);
Rstate(1,:) = 1;

for it = 2:Nt
    t = tGrid(it);

    if t < T-1e-12
        for is = 1:Ns
            Rstate(it,is) = compoundCdfRate( ...
                par.L-states(is),t,compound);
        end
    else
        pastIndex = it-nDelay;

        if pastIndex < 1
            error('Invalid fixed-delay reliability index.');
        end

        Rstate(it,:) = ...
            (Psurv*Rstate(pastIndex,:).').';
    end

    Rstate(it,:) = min(max(real(Rstate(it,:)),0),1);
end

% Remove only very small numerical upward reversals.
Rstate = cummin(Rstate,1);
Rstate(1,:) = 1;
end

function R = enforceReliabilityVector(R)
R = min(max(real(R),0),1);
R = cummin(R);
R(1) = 1;
end

%% ========================================================================
%  Coarse-grid search
% ========================================================================
function tab = gridSearchPeriodicPolicy( ...
    TGrid,MGrid,par,opts,compound,label)

[TList,MList] = ndgrid(TGrid(:),MGrid(:));
TList = TList(:);
MList = MList(:);
n = numel(TList);

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

usePar = opts.useParallel && ...
    license('test','Distrib_Computing_Toolbox');

if usePar
    pool = gcp('nocreate');

    if isempty(pool)
        fprintf('Starting parallel pool...\n');
        pool = parpool;
    end

    fprintf('Parallel pool active with %d workers.\n',pool.NumWorkers);

    dq = parallel.pool.DataQueue;
    progress = makeProgressPrinterPeriodic( ...
        n,label,opts.progressEverySeconds);

    afterEach(dq,@(~) progress());

    parfor k = 1:n
        try
            rep = evaluatePeriodicPolicy( ...
                TList(k),MList(k),par,opts,compound,false);

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

        send(dq,k);
    end
else
    fprintf('Running serial periodic-policy grid search.\n');

    t0 = tic;
    lastPrint = 0;

    for k = 1:n
        try
            rep = evaluatePeriodicPolicy( ...
                TList(k),MList(k),par,opts,compound,false);

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
                label,k,n,100*k/n,toc(t0));
        end
    end
end

tab = table( ...
    TList,MList,CostRate,Availability,Ecycle,MTTFnew,MTTFstat, ...
    Ppm,Ppr,Pcr,Down,Err, ...
    'VariableNames',{ ...
    'T','M','CostRate','Availability','ExpectedCycleLength', ...
    'MTTF_fromNew','MTTF_stationary','ExpectedPM','ExpectedPR', ...
    'ExpectedCR','ExpectedDowntime','Error'});
end

function progress = makeProgressPrinterPeriodic( ...
    n,label,printEverySeconds)

count = 0;
t0 = tic;
lastPrintTime = 0;
progress = @update;

    function update()
        count = count+1;
        elapsed = toc(t0);

        if elapsed-lastPrintTime >= printEverySeconds || count == n
            lastPrintTime = elapsed;
            rate = count/max(elapsed,eps);
            remain = (n-count)/max(rate,eps);

            fprintf(['[%s] %4d/%4d policies evaluated (%.1f%%), ', ...
                'elapsed %.1fs, ETA %.1fs\n'], ...
                label,count,n,100*count/n,elapsed,remain);
        end
    end
end

%% ========================================================================
%  Compound-degradation distribution table
% ========================================================================
function compound = buildCompoundDistributionTable(par,opts)

Nt = opts.NtCompound;
Nu = opts.NuCompound;
Nfft = 2^nextpow2(opts.NfftCompound);
tMaximum = opts.tCompoundMax;

if Nt < 2 || Nu < 2
    error('NtCompound and NuCompound must be at least 2.');
end

meanRate = ...
    par.alpha/par.beta + ...
    par.lambda*par.alphaS/par.betaS;

varianceRate = ...
    par.alpha/par.beta^2 + ...
    par.lambda*( ...
    par.alphaS/par.betaS^2 + ...
    (par.alphaS/par.betaS)^2);

meanMaximum = meanRate*tMaximum;
sdMaximum = sqrt(max(varianceRate*tMaximum,0));

fftUpperBound = max([ ...
    10*par.L, ...
    meanMaximum+12*sdMaximum+2*par.L, ...
    50]);

du = fftUpperBound/Nfft;
fftEdges = (0:Nfft)*du;
fftCenters = ((1:Nfft)-0.5)*du; %#ok<NASGU>

shockMass = diff(gammaCdfRate( ...
    fftEdges,par.alphaS,par.betaS));

shockMass = max(real(shockMass),0);
shockTransform = fft(shockMass,Nfft);

tGrid = linspace(0,tMaximum,Nt);
uGrid = linspace(0,par.L,Nu);

uBinIndex = floor(uGrid/du)+1;
uBinIndex = min(max(uBinIndex,1),Nfft);

cdfRaw = zeros(Nt,Nu);
pdfRaw = zeros(Nt,Nu);
retainedMass = zeros(Nt,1);
upperTailMass = zeros(Nt,1);

nTail = max(10,round(0.02*Nfft));

for it = 1:Nt
    t = tGrid(it);

    if t == 0
        cdfRaw(it,:) = 1;
        cdfRaw(it,1) = 0;
        pdfRaw(it,:) = 0;
        retainedMass(it) = 1;
        upperTailMass(it) = 0;
        continue;
    end

    continuousMass = diff(gammaCdfRate( ...
        fftEdges,par.alpha*t,par.beta));

    continuousMass = max(real(continuousMass),0);

    compoundShockTransform = ...
        exp(par.lambda*t*(shockTransform-1));

    totalMass = real(ifft( ...
        fft(continuousMass,Nfft).*compoundShockTransform));

    totalMass(abs(totalMass) < 1e-15) = 0;
    totalMass(totalMass < 0) = 0;

    retainedMass(it) = sum(totalMass);
    upperTailMass(it) = sum(totalMass(end-nTail+1:end));

    cdfAtEdges = [0,cumsum(totalMass)];
    pdfAtCenters = totalMass/du;

    cdfRaw(it,:) = interp1( ...
        fftEdges,cdfAtEdges,uGrid,'linear',0);

    pdfRaw(it,:) = pdfAtCenters(uBinIndex);

    cdfRaw(it,:) = min(max(cdfRaw(it,:),0),1);
    pdfRaw(it,:) = max(pdfRaw(it,:),0);
end

compound.tGrid = tGrid;
compound.uGrid = uGrid;
compound.cdfRaw = cdfRaw;
compound.pdfRaw = pdfRaw;

compound.fftUpperBound = fftUpperBound;
compound.fftStep = du;
compound.totalRetainedMass = retainedMass;
compound.upperTailMass = upperTailMass;
end

function F = compoundCdfRate(x,t,compound)
[xb,tb] = commonSize(x,t);
F = zeros(size(xb));

for k = 1:numel(xb)
    if xb(k) <= 0
        F(k) = 0;
    elseif tb(k) <= 0
        F(k) = 1;
    else
        if tb(k) > compound.tGrid(end)+1e-12
            error('Requested compound CDF time %.6g exceeds table maximum %.6g.', ...
                tb(k),compound.tGrid(end));
        end

        tt = min(max(tb(k),compound.tGrid(1)),compound.tGrid(end));
        xx = min(xb(k),compound.uGrid(end));

        row = interp1( ...
            compound.tGrid,compound.cdfRaw,tt,'linear');

        F(k) = interp1( ...
            compound.uGrid,row,xx,'linear');

        if xb(k) > compound.uGrid(end)
            F(k) = row(end);
        end
    end
end

F = min(max(real(F),0),1);
end

function pdf = compoundPdfRate(x,t,compound)
pdf = zeros(size(x));

if t <= 0
    return;
end

if t > compound.tGrid(end)+1e-12
    error('Requested compound PDF time %.6g exceeds table maximum %.6g.', ...
        t,compound.tGrid(end));
end

tt = min(max(t,compound.tGrid(1)),compound.tGrid(end));
row = interp1(compound.tGrid,compound.pdfRaw,tt,'linear');

valid = x > 0 & x <= compound.uGrid(end);

if any(valid(:))
    pdf(valid) = interp1( ...
        compound.uGrid,row,x(valid),'linear',0);
end

pdf = max(real(pdf),0);
end

function value = compoundReliabilityIntegral( ...
    margin,T,compound,Nt)

if T <= 0
    value = 0;
    return;
end

t = ((1:Nt)' - 0.5)/Nt*T;
wt = T/Nt;

F = compoundCdfRate( ...
    ones(size(t))*margin,t,compound);

value = sum(wt.*F);
value = min(max(value,0),T);
end

%% ========================================================================
%  Gamma and Beta numerical functions
% ========================================================================
function pdf = gammaPdfRate(x,a,beta)
[xb,ab] = commonSize(x,a);
pdf = zeros(size(xb));

idx = xb > 0 & ab > 0;

if any(idx(:))
    aa = ab(idx);
    xx = xb(idx);

    logp = aa.*log(beta)-gammaln(aa)+ ...
        (aa-1).*log(xx)-beta.*xx;

    pdf(idx) = exp(logp);
end

pdf(~isfinite(pdf)) = 0;
end

function cdf = gammaCdfRate(x,a,beta)
[xb,ab] = commonSize(x,a);
cdf = zeros(size(xb));

idxDeg = xb > 0 & ab <= 0;
cdf(idxDeg) = 1;

idx = xb > 0 & ab > 0;

if any(idx(:))
    cdf(idx) = gammainc( ...
        beta.*xb(idx),ab(idx),'lower');
end

cdf(isinf(xb) & xb > 0) = 1;
cdf = min(max(cdf,0),1);
end

function [xb,ab] = commonSize(x,a)
if isscalar(x) && isscalar(a)
    xb = x;
    ab = a;
elseif isscalar(x)
    ab = a;
    xb = ones(size(a)).*x;
elseif isscalar(a)
    xb = x;
    ab = ones(size(x)).*a;
else
    if ~isequal(size(x),size(a))
        error('Inputs must be scalar or have identical sizes.');
    end

    xb = x;
    ab = a;
end
end

function Bsurv = betaSurvRate(x,p,q)
Bsurv = zeros(size(x));
Bsurv(x <= 0) = 1;
Bsurv(x >= 1) = 0;

idx = x > 0 & x < 1;

if any(idx(:))
    Bsurv(idx) = 1-betainc(x(idx),p,q);
end

Bsurv = min(max(Bsurv,0),1);
end

function pdf = betaPdfUnit(x,p,q)
pdf = zeros(size(x));
idx = x > 0 & x < 1;

if any(idx(:))
    xx = x(idx);

    pdf(idx) = exp( ...
        (p-1).*log(xx)+ ...
        (q-1).*log1p(-xx)- ...
        betaln(p,q));
end

pdf(~isfinite(pdf)) = 0;
end

function mat = pmPdfMatrix(z,v,p,q)
z = z(:);
v = v(:)';

ratio = z./v;
mat = betaPdfUnit(ratio,p,q)./v;

mat(ratio <= 0 | ratio >= 1) = 0;
mat(~isfinite(mat)) = 0;
end

%% ========================================================================
%  Reporting and plotting
% ========================================================================
function printPeriodicPolicyReport(rep,titleStr)

fprintf('\n------------------------------------------------------------\n');
fprintf('%s\n',titleStr);
fprintf('------------------------------------------------------------\n');

fprintf('T = %.8g, M = %.8g, M/L = %.4f\n', ...
    rep.policy.T,rep.policy.M,rep.policy.M/rep.parameters.L);

fprintf('Cost rate C_inf              = %.10g\n',rep.costRate);
fprintf('Expected cycle length        = %.10g\n',rep.expectedCycleLength);
fprintf('Expected downtime/cycle      = %.10g\n',rep.expectedDowntimePerCycle);
fprintf('Average availability         = %.10f\n',rep.averageAvailability);
fprintf('MTTF, stationary initial pi  = %.10g\n',rep.MTTF_stationary);
fprintf('MTTF, new state s=0          = %.10g\n',rep.MTTF_fromNew);
fprintf('Stationary mass at state 0   = %.10f\n',rep.stationary.a0);

fprintf('\nExpected event counts per cycle:\n');
fprintf('  inspection                 = %.10g\n', ...
    rep.expectedCountsPerCycle.inspection);
fprintf('  PM                         = %.10g\n', ...
    rep.expectedCountsPerCycle.PM);
fprintf('  preventive replacement     = %.10g\n', ...
    rep.expectedCountsPerCycle.preventiveReplacement);
fprintf('  corrective replacement     = %.10g\n', ...
    rep.expectedCountsPerCycle.correctiveReplacement);

fprintf('\nCost-rate components and shares:\n');
fprintf('  inspection                 = %-12.6g  share = %6.2f%%\n', ...
    rep.costRateComponents.inspection,100*rep.costShare.inspection);
fprintf('  PM                         = %-12.6g  share = %6.2f%%\n', ...
    rep.costRateComponents.PM,100*rep.costShare.PM);
fprintf('  preventive replacement     = %-12.6g  share = %6.2f%%\n', ...
    rep.costRateComponents.preventiveReplacement, ...
    100*rep.costShare.preventiveReplacement);
fprintf('  corrective replacement     = %-12.6g  share = %6.2f%%\n', ...
    rep.costRateComponents.correctiveReplacement, ...
    100*rep.costShare.correctiveReplacement);
fprintf('  downtime                   = %-12.6g  share = %6.2f%%\n', ...
    rep.costRateComponents.downtime,100*rep.costShare.downtime);

fprintf('\nDiagnostics:\n');
fprintf('  max RAW transition error   = %.3e\n', ...
    rep.diagnostics.maxRawTransitionMassError);
fprintf('  max corrected error        = %.3e\n', ...
    rep.diagnostics.maxCorrectedTransitionMassError);
fprintf('  max do-nothing mass error  = %.3e\n', ...
    rep.diagnostics.maxDoNothingMassError);
fprintf('  max PM-success mass error  = %.3e\n', ...
    rep.diagnostics.maxPMSuccessMassError);
fprintf('  max PM-event mass error    = %.3e\n', ...
    rep.diagnostics.maxPMEventMassError);
fprintf('  rcond(I-P_survival)        = %.3e\n', ...
    rep.diagnostics.minMTTFLinearMatrixReciprocalCondition);

if ~isempty(rep.reliability.t)
    fprintf('\nReliability curve generated with %d time points.\n', ...
        numel(rep.reliability.t));
    fprintf('  R_new(t_end)               = %.6g\n', ...
        rep.reliability.RfromNew(end));
    fprintf('  R_stationary(t_end)        = %.6g\n', ...
        rep.reliability.Rstationary(end));
    fprintf('  integral R_stationary      = %.10g\n', ...
        rep.diagnostics.MTTFfromReliabilityCurve);
    fprintf('  relative MTTF curve error  = %.3e\n', ...
        rep.diagnostics.relativeMTTFcurveError);
end
end

function checkBoundaryPeriodicPolicy( ...
    bestT,TGrid,bestM,MGrid,label)

if abs(bestT-min(TGrid)) < 1e-12 || ...
   abs(bestT-max(TGrid)) < 1e-12 || ...
   abs(bestM-min(MGrid)) < 1e-12 || ...
   abs(bestM-max(MGrid)) < 1e-12

    fprintf(['\nWarning: %s-grid optimum is on the search boundary. ', ...
        'Consider enlarging the grid range.\n'],label);
end
end

function plotPeriodicCostContour(tab,par,repBest)

valid = isfinite(tab.CostRate);

if ~any(valid)
    warning('No valid periodic-policy results are available for contour plotting.');
    return;
end

tab = tab(valid,:);

Tvals = unique(tab.T);
Mvals = unique(tab.M);
xVals = Mvals/par.L;

Z = nan(numel(Tvals),numel(Mvals));

for i = 1:numel(Tvals)
    for j = 1:numel(Mvals)
        idx = abs(tab.T-Tvals(i)) < 1e-12 & ...
              abs(tab.M-Mvals(j)) < 1e-12;

        if any(idx)
            Z(i,j) = tab.CostRate(find(idx,1));
        end
    end
end

finiteZ = Z(isfinite(Z));

if isempty(finiteZ)
    return;
end

zMin = min(finiteZ);
zMax = max(finiteZ);

if abs(zMax-zMin) < eps
    return;
end

nLevels = 26;
u = linspace(0,1,nLevels).^2.2;
levels = unique(zMin+(zMax-zMin)*u);

fig = figure( ...
    'Name','Periodic-policy coarse-grid cost contour', ...
    'Color','w', ...
    'Position',[100,100,780,560]);

ax = axes('Parent',fig);
hold(ax,'on');

[C,h] = contour( ...
    ax,xVals,Tvals,Z,levels,'LineWidth',0.90);

colormap(ax,jet(256));

cb = colorbar(ax);
ylabel(cb,'Long-run cost rate');

caxis(ax,[zMin,zMax]);

labelLevels = levels(3:4:end);

if ~isempty(labelLevels)
    clabel(C,h,labelLevels, ...
        'FontSize',8, ...
        'Color','k', ...
        'LabelSpacing',520);
end

xBest = repBest.policy.M/par.L;
yBest = repBest.policy.T;

plot(ax,xBest,yBest,'kp', ...
    'MarkerFaceColor','y', ...
    'MarkerSize',12, ...
    'LineWidth',1.2);

text(ax,xBest,yBest, ...
    sprintf('  optimum %.4g',repBest.costRate), ...
    'FontSize',9, ...
    'FontWeight','bold', ...
    'VerticalAlignment','bottom');

xlabel(ax,'Preventive maintenance threshold M/L');
ylabel(ax,'Inspection interval T');
title(ax,'Equal-interval policy long-run cost-rate contour');

box(ax,'on');
set(ax,'Layer','top');

savefig(fig,'periodic_policy_cost_contour_v1.fig');

try
    exportgraphics( ...
        fig,'periodic_policy_cost_contour_v1.png','Resolution',300);
catch
    print(fig,'periodic_policy_cost_contour_v1.png','-dpng','-r300');
end

fprintf(['Saved periodic-policy contour plot to ', ...
    'periodic_policy_cost_contour_v1.fig and ', ...
    'periodic_policy_cost_contour_v1.png\n']);
end

function plotPeriodicReliabilityCurve(repBest)

if ~isfield(repBest,'reliability') || ...
        isempty(repBest.reliability.t)

    fprintf('Periodic-policy reliability curve is empty; plot skipped.\n');
    return;
end

fig = figure( ...
    'Name','Periodic-policy reliability function', ...
    'Color','w');

ax = axes('Parent',fig);

plot(ax, ...
    repBest.reliability.t, ...
    repBest.reliability.Rstationary, ...
    'LineWidth',1.5);

hold(ax,'on');

plot(ax, ...
    repBest.reliability.t, ...
    repBest.reliability.RfromNew, ...
    '--','LineWidth',1.5);

grid(ax,'on');

xlabel(ax,'Time');
ylabel(ax,'Reliability');

legend(ax, ...
    'stationary initial distribution', ...
    'new state s=0', ...
    'Location','best');

title(ax,'Reliability under optimal equal-interval policy');

savefig(fig,'periodic_policy_reliability_v1.fig');

try
    exportgraphics( ...
        fig,'periodic_policy_reliability_v1.png','Resolution',300);
catch
    print(fig,'periodic_policy_reliability_v1.png','-dpng','-r300');
end

fprintf(['Saved periodic-policy reliability plot to ', ...
    'periodic_policy_reliability_v1.fig and ', ...
    'periodic_policy_reliability_v1.png\n']);
end

%% ========================================================================
%  Validation
% ========================================================================
function validatePeriodicPolicy(T,M,par,compound)

if ~(T > 0)
    error('Inspection interval T must satisfy T>0.');
end

if T > compound.tGrid(end)+1e-12
    error('T=%.6g exceeds compound-table maximum %.6g.', ...
        T,compound.tGrid(end));
end

if ~(M > 0 && M < par.L)
    error('M must satisfy 0<M<L.');
end
end
