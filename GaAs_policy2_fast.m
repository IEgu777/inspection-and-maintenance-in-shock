% GaAs LASER CASE (Example 2) -- POLICY2
% Adapted directly from run_policy2_case_v3_fast.m; local numerical functions are retained.
%
% Time: 1 model unit = 250 operating hours.
% Degradation: Gamma process alpha=7.426 and rate beta=1/0.0674.
% Shock size: Gamma(shape, RATE) parameters re-estimated from NEW HPP
%   shock arrival observations; read via GaAs_case_parameters.m.
% Injected shocks follow Bernoulli grid-event sampling (p=0.25); the
% policy model uses estimated mean event frequency as HPP approximation.
% Beta imperfect PM: p=6, q=5 (assumed, not measured for GaAs).
% Cost inputs are centralized in GaAs_case_parameters.m. The default
%   relative cost structure is anchored at Cp=50 (Zhang et al., 2014)
%   and uses Section 6.1 ratios; no real cost invoices are available.
%
% Run mode: 'quick' for a pilot; 'full' for a 0.01 resolution paper grid.
% Accuracy: the grid optimum is approximate, not proof of global optimality.
% After optimization, examine mass-error diagnostics and boundary warnings.
% Requires MATLAB; a parallel toolbox is optional.

% FAST VERSION: preserves Policy 2 equations, with GaAs case inputs.
% Optimizations: skip MTTF linear solve during cost search; final optimum retains it.
% Other accelerations are policy-specific and documented below.
% Verify results by comparing quick/full grids and quadrature refinements.
% GAAS_POLICY2_FAST
% Policy 2 benchmark for the GaAs laser case in Section 6.2.
%
% Policy 2:
%   Inspections are conducted only at shock arrival times. This is the
%   boundary policy corresponding to RI=0, and RI is not optimized.
%
% At each inspection, the same CBM rules as the proposed policy are used.
%
% Main features:
%   1) Only M is optimized by a coarse-grid search.
%   2) The final policy is evaluated directly at the coarse-grid optimum.
%   3) This file is a script, not a main function.
%   4) Runtime is recorded for every local function called by the script.
%   5) The reliability curve starts from the stationary post-maintenance
%      distribution of Policy 2.
%
% V2 corrections:
%   - do-nothing and successful-PM transition densities are corrected
%     separately to their theoretical probability masses;
%   - RAW transition errors are retained as diagnostics;
%   - the reliability function is solved on a uniform Volterra grid using
%     a precomputed first-shock survival kernel;
%   - the reliability horizon is based on the stationary MTTF.
%
% How to run:
%   >> GaAs_policy2_fast
%
% Main output variables:
%   results       - complete Policy 2 result structure
%   coarseTab     - coarse-grid search table
%   repBest       - final evaluation at the coarse-grid optimum
%   timingLog     - runtime table

clearvars;
clc;
close all;

% Switch to 'full' before reporting final numerical results in the paper.
runMode = getenv('GAAS_POLICY_RUN_MODE');
if isempty(runMode)
    runMode = 'full'; % Matches the full-grid paper experiment.
end
assert(any(strcmpi(runMode, {'quick','full'})), ...
    'runMode must be quick or full.');

scriptTimer = tic;
timingNames = {};
timingSeconds = [];

fprintf('\n============================================================\n');
fprintf(' GaAs laser: Policy 2 (shock-triggered only)\n');
fprintf('============================================================\n');

%% --------------------- Model and cost parameters -----------------------
% Shared by all four policy implementations; change GaAs_case_parameters.m.
par = GaAs_case_parameters();

%% --------------------- Numerical options --------------------
opts.Nz       = 60;     % grid points on post-maintenance interval (0,M)
opts.Nv       = 80;     % quadrature points for v in [M,L]
opts.Nlag     = 48;     % Gauss-Laguerre points for Exp(lambda) integration

% Reliability recursion on a uniform time grid.
opts.NtRel          = 501;   % minimum reliability time-grid size
opts.maxNtRel       = 2200;  % maximum number of time points
opts.dtRelMaximum   = 0.25;  % desired largest time step
opts.tRelMaxFactor  = 6.0;   % tMax >= factor*max(MTTFs,Ecycle)
opts.kernelTailProbability = 1e-8;
opts.reliabilityTailWarning = 1e-5;

% Distribution of X_Y+W at a shock-triggered inspection, Y~Exp(lambda).
opts.NuInspection   = 1201;
opts.NfftInspection = 8192;
opts.inspectionUpperFactor = 12;

opts.useParallel = true;
opts.showProgress = true;
opts.progressEverySeconds = 3;

% This threshold is applied to RAW, pre-correction quadrature errors.
opts.minMassWarning = 1e-3;

%% --------------------- Precompute inspection distribution ---
fprintf('\n[0] Precomputing the distribution of X_Y+W, Y~Exp(lambda)...\n');

functionTimer = tic;
inspectionDist = buildPolicy2InspectionDistribution(par, opts);
timingNames{end+1,1} = 'buildPolicy2InspectionDistribution';
timingSeconds(end+1,1) = toc(functionTimer);

fprintf('Inspection distribution completed.\n');
fprintf('Retained probability at FFT upper bound = %.10f\n', ...
    inspectionDist.totalRetainedMass);
fprintf('Probability in the upper 2%% of the FFT grid = %.3e\n', ...
    inspectionDist.upperTailMass);

%% --------------------- Coarse grid search only --------------
if strcmpi(runMode, 'quick')
    coarse.MGrid = (0.02:0.02:0.98) * par.L;
else
    coarse.MGrid = (0.01:0.01:0.99) * par.L;
end

fprintf('\n[1] Policy 2 coarse grid search: %d policies\n', ...
    numel(coarse.MGrid));

functionTimer = tic;
coarseTab = gridSearchPolicy2( ...
    coarse.MGrid, par, opts, inspectionDist, 'policy2 coarse');
timingNames{end+1,1} = 'gridSearchPolicy2';
timingSeconds(end+1,1) = toc(functionTimer);

validCost = isfinite(coarseTab.CostRate);
if ~any(validCost)
    error('No valid Policy 2 candidate was obtained. Inspect coarseTab.Error.');
end

validIndices = find(validCost);
[~, localBest] = min(coarseTab.CostRate(validCost));
idxBest = validIndices(localBest);

bestPolicy.RI = NaN;
bestPolicy.M = coarseTab.M(idxBest);
bestPolicy.CostRate = coarseTab.CostRate(idxBest);

fprintf('\n[2] Re-evaluating the Policy 2 optimum with reliability curve...\n');

functionTimer = tic;
repBest = evaluatePolicy2( ...
    bestPolicy.M, par, opts, inspectionDist, false, true);
timingNames{end+1,1} = 'evaluatePolicy2 at coarse optimum';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
printPolicy2Report(repBest, 'Final Policy 2 from coarse-grid optimization');
timingNames{end+1,1} = 'printPolicy2Report';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
checkBoundaryPolicy2( ...
    bestPolicy.M, coarse.MGrid, 'Policy 2 coarse');
timingNames{end+1,1} = 'checkBoundaryPolicy2';
timingSeconds(end+1,1) = toc(functionTimer);

% To compute a final reliability curve, enable the optional flag in the
% final evaluation call (a computationally costly calculation).
%% --------------------- Plots --------------------------------
functionTimer = tic;
try
    plotPolicy2CostCurve(coarseTab, par, repBest);
catch ME
    warning('Policy 2 cost-curve plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotPolicy2CostCurve';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
try
    plotPolicy2ReliabilityCurve(repBest);
catch ME
    warning('Policy 2 reliability plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotPolicy2ReliabilityCurve';
timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Results and timing --------------------
timingNames{end+1,1} = 'total script time before save';
timingSeconds(end+1,1) = toc(scriptTimer);

timingLog = table( ...
    timingNames, timingSeconds, ...
    'VariableNames', {'FunctionOrStep','Seconds'});

fprintf('\n---------------- Runtime log for main-script function calls ----------------\n');
disp(timingLog);

results = struct();
results.caseName = 'GaAs laser, Example 2';
results.runMode = runMode;
results.timeUnitHours = 250;
results.costUnit = 'dimensionless relative cost';
results.shockArrivalNote = 'lambda is estimated grid-event frequency; HPP frequency-matching approximation.';
results.policyName      = 'Policy 2';
results.parameters      = par;
results.options         = opts;
results.coarseGrid      = coarse;
results.coarseGridTable = coarseTab;
results.bestPolicy      = bestPolicy;
results.finalBest       = repBest;
results.timingLog       = timingLog;
results.inspectionDistributionDiagnostics = struct( ...
    'fftUpperBound', inspectionDist.fftUpperBound, ...
    'fftStep', inspectionDist.fftStep, ...
    'totalRetainedMass', inspectionDist.totalRetainedMass, ...
    'upperTailMass', inspectionDist.upperTailMass);

save('GaAs_policy2_results.mat', 'results');
fprintf('\nSaved Policy 2 results to GaAs_policy2_results.mat\n');
fprintf('Done.\n');

%% ========================================================================
%  Main Policy 2 evaluator
% ========================================================================
function rep = evaluatePolicy2( ...
    M, par, opts, inspectionDist, doReliabilityCurve, doMTTF)

if nargin < 5
    doReliabilityCurve = false;
end
if nargin < 6, doMTTF = true; end
if doReliabilityCurve && ~doMTTF, error('Reliability curve requires MTTF.'); end

validatePolicy2(M, par);

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

pmSurvVec = betaSurvRate(M ./ v, par.p, par.q);
pmPdfMat  = pmPdfMatrix(z, v, par.p, par.q);

% Transition kernels and conditional cycle quantities.
Kdens = zeros(Ns,Nz);
Q1full = zeros(Ns,1);
Psurv = zeros(Ns,Ns);

% Diagnostics are recorded BEFORE and AFTER component-wise correction.
PrawMassError = zeros(Ns,1);
PcorrectedMassError = zeros(Ns,1);
PdoNothingMassError = zeros(Ns,1);
PpmSuccessMassError = zeros(Ns,1);
PpmEventMassError = zeros(Ns,1);

cond = initConditionalStructPolicy2(Ns);

for i = 1:Ns
    row = buildPolicy2OneStepRow( ...
        states(i), M, par, opts, inspectionDist, ...
        z, wz, v, wv, pmSurvVec, pmPdfMat);

    Kdens(i,:) = row.fullInteriorDensity;
    Q1full(i)  = row.fullMass0;
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
% pi(ds) = a0*delta0(ds) + a0*H(s)ds.
A = zeros(Nz,Nz);

for iz = 1:Nz
    for js = 1:Nz
        A(iz,js) = Kdens(js+1,iz) * wz(js);
    end
end

b = Kdens(1,:)';
H = (eye(Nz)-A) \ b;

H(H < 0 & H > -1e-10) = 0;
if any(H < -1e-8)
    warning('Policy 2 stationary density H contains noticeable negative values.');
end
H = max(H,0);

a0 = 1 / (1 + sum(H .* wz));
piVec = [a0; a0 * H .* wz];
piVec = piVec / sum(piVec);

% Expected cycle quantities.
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
share = [costInspection,costPM,costPR,costCR,costDowntime] ./ Cinf;

availability = 1-Edown/Ecycle;

% MTTF:
% m = source + P_survival*m.
I = eye(Ns);
MTTF_stationary = NaN;
MTTF_fromNew = NaN;
mtRcond = NaN;
if doMTTF
    mVec = (I-Psurv) \ cond.mttfSource;
    MTTF_stationary = sum(piVec .* mVec);
    MTTF_fromNew = mVec(1);
    mtRcond = rcond(I-Psurv);
end

% Reliability curve from the stationary distribution and from state 0.
rel = struct('t', [], 'Rstationary', [], 'RfromNew', [], 'RbyState', []);

if doReliabilityCurve
    tMax = opts.tRelMaxFactor * ...
        max([MTTF_stationary,MTTF_fromNew,Ecycle]);

    if ~isfinite(tMax) || tMax <= 0
        tMax = 1e3;
    end

    tMax = min(tMax,1e6);

    NtRel = max(2,round(opts.NtRel));
    NtRequired = ceil(tMax/max(opts.dtRelMaximum,eps))+1;
    NtRel = max(NtRel,NtRequired);

    if NtRel > opts.maxNtRel
        warning(['Policy 2 reliability grid was capped at %d points. ', ...
            'Increase opts.maxNtRel for a finer Volterra recursion.'], ...
            opts.maxNtRel);
        NtRel = opts.maxNtRel;
    end

    tGrid = linspace(0,tMax,NtRel);

    [Rstate,kernelDiagnostics] = reliabilityCurvePolicy2( ...
        tGrid, states, z, wz, v, wv, ...
        pmSurvVec, pmPdfMat, M, par, opts, inspectionDist);

    rel.t = tGrid(:);
    rel.RbyState = Rstate;
    rel.Rstationary = enforceReliabilityVector(Rstate*piVec);
    rel.RfromNew = enforceReliabilityVector(Rstate(:,1));
    rel.kernelDiagnostics = kernelDiagnostics;
end

rep = struct();
rep.policy.RI = NaN;
rep.policy.M  = M;

rep.costRate = Cinf;
rep.expectedCycleLength = Ecycle;

rep.expectedCountsPerCycle = struct( ...
    'inspection', 1, ...
    'PM', Epm, ...
    'preventiveReplacement', Epr, ...
    'correctiveReplacement', Ecr);

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

rep.stationary = struct( ...
    'a0', a0, ...
    'z', z, ...
    'wz', wz, ...
    'H', H, ...
    'piVec', piVec);

rep.conditional = cond;

maximumUncorrectedComponentError = max([ ...
    max(abs(PrawMassError)), ...
    max(abs(PdoNothingMassError)), ...
    max(abs(PpmSuccessMassError)), ...
    max(abs(PpmEventMassError))]);

rep.diagnostics = struct( ...
    'maxRawTransitionMassError', max(abs(PrawMassError)), ...
    'meanRawTransitionMassError', mean(abs(PrawMassError)), ...
    'maxCorrectedTransitionMassError', max(abs(PcorrectedMassError)), ...
    'maxDoNothingMassError', max(abs(PdoNothingMassError)), ...
    'maxPMSuccessMassError', max(abs(PpmSuccessMassError)), ...
    'maxPMEventMassError', max(abs(PpmEventMassError)), ...
    'maxUncorrectedComponentError', maximumUncorrectedComponentError, ...
    'minMTTFLinearMatrixReciprocalCondition', mtRcond);

if doReliabilityCurve && ~isempty(rel.t)
    rep.diagnostics.MTTFfromReliabilityCurve = ...
        trapz(rel.t,rel.Rstationary);
    rep.diagnostics.reliabilityTail = rel.Rstationary(end);
    rep.diagnostics.relativeMTTFcurveError = ...
        abs(rep.diagnostics.MTTFfromReliabilityCurve-MTTF_stationary) / ...
        max(MTTF_stationary,eps);
    rep.diagnostics.reliabilityKernel = rel.kernelDiagnostics;
else
    rep.diagnostics.MTTFfromReliabilityCurve = NaN;
    rep.diagnostics.reliabilityTail = NaN;
    rep.diagnostics.relativeMTTFcurveError = NaN;
    rep.diagnostics.reliabilityKernel = struct();
end

if rep.diagnostics.maxUncorrectedComponentError > opts.minMassWarning
    warning(['Policy 2 maximum uncorrected component error is %.3g. ', ...
        'Increase opts.Nz, opts.Nv, or the inspection grids.'], ...
        rep.diagnostics.maxUncorrectedComponentError);
end

if doReliabilityCurve && ...
        rep.diagnostics.reliabilityTail > opts.reliabilityTailWarning
    warning(['Policy 2 reliability at the final time is %.3g. ', ...
        'Increase opts.tRelMaxFactor or opts.maxNtRel.'], ...
        rep.diagnostics.reliabilityTail);
end
end

%% ========================================================================
%  Policy 2 one-step kernel
% ========================================================================
function row = buildPolicy2OneStepRow( ...
    s, M, par, opts, inspectionDist, ...
    z, wz, v, wv, pmSurvVec, pmPdfMat)

components = buildPolicy2TransitionComponents( ...
    s,M,par,inspectionDist, ...
    z,wz,v,wv,pmSurvVec,pmPdfMat);

% The time to the next inspection is Exp(lambda).
mu = 1/par.lambda;

% Source term:
% integral_0^infinity exp(-lambda*t)F_X(L-s;t)dt.
mttfSource = policy2OperatingTime( ...
    par.L-s, par, ...
    inspectionDist.laguerreX, ...
    inspectionDist.laguerreW);

downtime = max(mu-mttfSource,0);

row = components;
row.mu = mu;
row.downtime = downtime;
row.mttfSource = mttfSource;
end


function row = buildPolicy2TransitionComponents( ...
    s,M,par,incrementDist, ...
    z,wz,v,wv,pmSurvVec,pmPdfMat)

L = par.L;

% Distribution of the degradation increment at the next shock.
g_v = inspectionPdfRate(v-s,incrementDist);
g_v(v <= s) = 0;

g_z = inspectionPdfRate(z-s,incrementDist);
g_z(z <= s) = 0;

F_M = inspectionCdfRate(M-s,incrementDist);
F_L = inspectionCdfRate(L-s,incrementDist);

F_M = min(max(F_M,0),1);
F_L = min(max(F_L,F_M),1);

Pcr = 1-F_L;
Ppm = F_L-F_M;

PprRaw = sum(wv .* pmSurvVec .* g_v);

pmSuccessProbabilityRaw = ...
    sum(wv .* (1-pmSurvVec) .* g_v);

pmEventMassRaw = ...
    PprRaw+pmSuccessProbabilityRaw;

if Ppm <= 0
    Ppr = 0;
    pmSuccessTarget = 0;
elseif pmEventMassRaw > 0
    pmEventScale = Ppm/pmEventMassRaw;
    Ppr = max(PprRaw,0)*pmEventScale;
    pmSuccessTarget = ...
        max(pmSuccessProbabilityRaw,0)*pmEventScale;
else
    Ppr = Ppm;
    pmSuccessTarget = 0;
end

Ppr = min(max(Ppr,0),Ppm);
pmSuccessTarget = max(Ppm-Ppr,0);

pmSuccessDensityRaw = pmPdfMat * (wv .* g_v);
doNothingDensityRaw = g_z;

doNothingTarget = F_M;

doNothingMassRaw = sum(doNothingDensityRaw .* wz);
pmSuccessMassRaw = sum(pmSuccessDensityRaw .* wz);

doNothingMassError = ...
    doNothingMassRaw-doNothingTarget;

pmSuccessMassError = ...
    pmSuccessMassRaw-pmSuccessProbabilityRaw;

pmEventMassError = ...
    pmEventMassRaw-Ppm;

mass0Full = Pcr+Ppr;

rawMassError = ...
    Pcr+PprRaw+doNothingMassRaw+pmSuccessMassRaw-1;

doNothingDensity = doNothingDensityRaw;

if doNothingTarget <= 0
    doNothingDensity(:) = 0;
elseif doNothingMassRaw > 0
    doNothingDensity = doNothingDensity * ...
        (doNothingTarget/doNothingMassRaw);
else
    index = find(z > s,1,'first');

    if isempty(index)
        index = numel(z);
    end

    doNothingDensity(:) = 0;
    doNothingDensity(index) = ...
        doNothingTarget/wz(index);
end

pmSuccessDensity = pmSuccessDensityRaw;

if pmSuccessTarget <= 0
    pmSuccessDensity(:) = 0;
elseif pmSuccessMassRaw > 0
    pmSuccessDensity = pmSuccessDensity * ...
        (pmSuccessTarget/pmSuccessMassRaw);
else
    fallbackDensity = pmPdfMat*wv;
    fallbackMass = sum(fallbackDensity .* wz);

    if fallbackMass > 0
        pmSuccessDensity = fallbackDensity * ...
            (pmSuccessTarget/fallbackMass);
    else
        pmSuccessDensity(:) = 0;
        pmSuccessDensity(1) = ...
            pmSuccessTarget/wz(1);
    end
end

interiorDensity = ...
    doNothingDensity+pmSuccessDensity;

correctedMassError = ...
    mass0Full+sum(interiorDensity .* wz)-1;

PsurvRow = zeros(1,numel(z)+1);
PsurvRow(1) = Ppr;
PsurvRow(2:end) = (interiorDensity .* wz)';

row = struct();
row.Ppm = Ppm;
row.Ppr = Ppr;
row.Pcr = Pcr;

row.fullMass0 = mass0Full;
row.fullInteriorDensity = interiorDensity(:)';
row.PsurvRow = PsurvRow;

row.rawMassError = rawMassError;
row.correctedMassError = correctedMassError;
row.doNothingMassError = doNothingMassError;
row.pmSuccessMassError = pmSuccessMassError;
row.pmEventMassError = pmEventMassError;
end


function cond = initConditionalStructPolicy2(Ns)
cond = struct();
cond.mu = zeros(Ns,1);
cond.downtime = zeros(Ns,1);
cond.Ppm = zeros(Ns,1);
cond.Ppr = zeros(Ns,1);
cond.Pcr = zeros(Ns,1);
cond.mttfSource = zeros(Ns,1);
end

function value = policy2OperatingTime( ...
    margin, par, laguerreX, laguerreW)

t = laguerreX/par.lambda;

survivalBeforeShock = gammaCdfRate( ...
    ones(size(t))*margin, ...
    par.alpha.*t, ...
    par.beta);

value = sum(laguerreW .* survivalBeforeShock) / par.lambda;
value = min(max(value,0),1/par.lambda);
end

%% ========================================================================
%  Policy 2 reliability recursion
% ========================================================================
function [Rstate,diagnostics] = reliabilityCurvePolicy2( ...
    tGrid, states, z, wz, v, wv, ...
    pmSurvVec, pmPdfMat, M, par, opts, inspectionDist)

Ns = numel(states);
Nt = numel(tGrid);

if Nt < 2
    error('Policy 2 reliability recursion requires at least two time points.');
end

dt = tGrid(2)-tGrid(1);

if any(abs(diff(tGrid)-dt) > 1e-10*max(1,dt))
    error('Policy 2 reliability recursion requires a uniform time grid.');
end

kernelMaximum = min( ...
    tGrid(end), ...
    -log(opts.kernelTailProbability)/par.lambda);

nKernel = max(1,floor(kernelMaximum/dt));
uKernel = (0:nKernel)'*dt;
kernelMaximum = uKernel(end);

% K(:,:,k) is the density of the non-failure transition kernel at u_k:
%
%   lambda*exp(-lambda*u_k)*P_survival(u_k).
K = zeros(Ns,Ns,nKernel+1);

kernelRawError = zeros(nKernel+1,1);
kernelCorrectedError = zeros(nKernel+1,1);
kernelDoNothingError = zeros(nKernel+1,1);
kernelPMSuccessError = zeros(nKernel+1,1);
kernelPMEventError = zeros(nKernel+1,1);
kernelRetainedMass = zeros(nKernel+1,1);
kernelUpperTailMass = zeros(nKernel+1,1);

fprintf('  Building Policy 2 reliability kernel at %d time points...\n', ...
    nKernel+1);

progressStep = max(1,round((nKernel+1)/20));

for ku = 1:nKernel+1
    u = uKernel(ku);

    conditionalDist = ...
        conditionalInspectionDistributionAtTime( ...
        u,par,inspectionDist);

    kernelRetainedMass(ku) = ...
        conditionalDist.totalRetainedMass;

    kernelUpperTailMass(ku) = ...
        conditionalDist.upperTailMass;

    coefficient = ...
        par.lambda*exp(-par.lambda*u);

    maximumRaw = 0;
    maximumCorrected = 0;
    maximumDoNothing = 0;
    maximumPMSuccess = 0;
    maximumPMEvent = 0;

    for is = 1:Ns
        components = buildPolicy2TransitionComponents( ...
            states(is),M,par,conditionalDist, ...
            z,wz,v,wv,pmSurvVec,pmPdfMat);

        K(is,:,ku) = ...
            coefficient*components.PsurvRow;

        maximumRaw = max( ...
            maximumRaw,abs(components.rawMassError));

        maximumCorrected = max( ...
            maximumCorrected, ...
            abs(components.correctedMassError));

        maximumDoNothing = max( ...
            maximumDoNothing, ...
            abs(components.doNothingMassError));

        maximumPMSuccess = max( ...
            maximumPMSuccess, ...
            abs(components.pmSuccessMassError));

        maximumPMEvent = max( ...
            maximumPMEvent, ...
            abs(components.pmEventMassError));
    end

    kernelRawError(ku) = maximumRaw;
    kernelCorrectedError(ku) = maximumCorrected;
    kernelDoNothingError(ku) = maximumDoNothing;
    kernelPMSuccessError(ku) = maximumPMSuccess;
    kernelPMEventError(ku) = maximumPMEvent;

    if mod(ku,progressStep) == 0 || ku == nKernel+1
        fprintf('    kernel progress: %5.1f%%\n', ...
            100*ku/(nKernel+1));
    end
end

Rstate = zeros(Nt,Ns);
Rstate(1,:) = 1;

% Trapezoidal discretization of the Volterra equation contains
% 0.5*dt*K(0)*R(t), so each time step uses an implicit linear solve.
implicitMatrix = ...
    eye(Ns)-0.5*dt*K(:,:,1);

if rcond(implicitMatrix) < 1e-12
    warning(['Policy 2 reliability implicit matrix is poorly conditioned. ', ...
        'Reduce opts.dtRelMaximum.']);
end

for it = 2:Nt
    step = it-1;
    t = tGrid(it);

    directSurvival = zeros(Ns,1);

    for is = 1:Ns
        directSurvival(is) = ...
            exp(-par.lambda*t) * ...
            gammaCdfRate( ...
            par.L-states(is), ...
            par.alpha*t, ...
            par.beta);
    end

    rhs = directSurvival;

    if step <= nKernel
        % Interior trapezoidal points 0<u<t.
        for k = 1:step-1
            rhs = rhs + ...
                dt*K(:,:,k+1)*Rstate(it-k,:)';
        end

        % Endpoint u=t, where R(0)=1.
        rhs = rhs + ...
            0.5*dt*K(:,:,step+1)*ones(Ns,1);

    else
        % The first-shock kernel is truncated at kernelMaximum. The omitted
        % probability is bounded by exp(-lambda*kernelMaximum).
        for k = 1:nKernel-1
            rhs = rhs + ...
                dt*K(:,:,k+1)*Rstate(it-k,:)';
        end

        rhs = rhs + ...
            0.5*dt*K(:,:,nKernel+1) ...
            * Rstate(it-nKernel,:)';
    end

    Rstate(it,:) = ...
        (implicitMatrix\rhs)';

    Rstate(it,:) = ...
        min(max(real(Rstate(it,:)),0),1);
end

% Remove very small numerical increases.
Rstate = cummin(Rstate,1);
Rstate(1,:) = 1;

diagnostics = struct();
diagnostics.dt = dt;
diagnostics.kernelMaximum = kernelMaximum;
diagnostics.omittedFirstShockTailBound = ...
    exp(-par.lambda*kernelMaximum);
diagnostics.rcondImplicitMatrix = ...
    rcond(implicitMatrix);
diagnostics.maxRawTransitionMassError = ...
    max(kernelRawError);
diagnostics.maxCorrectedTransitionMassError = ...
    max(kernelCorrectedError);
diagnostics.maxDoNothingMassError = ...
    max(kernelDoNothingError);
diagnostics.maxPMSuccessMassError = ...
    max(kernelPMSuccessError);
diagnostics.maxPMEventMassError = ...
    max(kernelPMEventError);
diagnostics.minimumConditionalRetainedMass = ...
    min(kernelRetainedMass);
diagnostics.maximumConditionalUpperTailMass = ...
    max(kernelUpperTailMass);
end


function dist = conditionalInspectionDistributionAtTime( ...
    t,par,inspectionDist)

Nfft = inspectionDist.Nfft;
Nlinear = inspectionDist.Nlinear;
edges = inspectionDist.fftEdges;
centers = inspectionDist.fftCenters;
du = inspectionDist.fftStep;

if t <= 0
    continuousMass = zeros(1,Nfft);
    continuousMass(1) = 1;
else
    continuousMass = diff(gammaCdfRate( ...
        edges,par.alpha*t,par.beta));

    continuousMass = max(real(continuousMass),0);
end

convolutionMass = real(ifft( ...
    fft(continuousMass,Nlinear) ...
    .* inspectionDist.shockTransform));

convolutionMass = convolutionMass(1:Nfft);
convolutionMass(abs(convolutionMass) < 1e-15) = 0;
convolutionMass(convolutionMass < 0) = 0;

cdfAtEdges = [0,cumsum(convolutionMass)];
pdfAtCenters = convolutionMass/du;

dist.uGrid = inspectionDist.uGrid;

dist.cdf = interp1( ...
    edges,cdfAtEdges,dist.uGrid,'linear',0);

% Use the same finite-volume bins as the inspection distribution.
dist.pdf = pdfAtCenters(inspectionDist.uBinIndex);

dist.cdf = min(max(dist.cdf,0),1);
dist.pdf = max(dist.pdf,0);

dist.totalRetainedMass = sum(convolutionMass);

nTail = max(10,round(0.02*Nfft));
dist.upperTailMass = ...
    sum(convolutionMass(end-nTail+1:end));
end


function R = enforceReliabilityVector(R)
R = min(max(real(R(:)),0),1);
R = cummin(R);
R(1) = 1;
end

%% ========================================================================
%  Grid search with parallel computing and progress display
% ========================================================================
function tab = gridSearchPolicy2( ...
    MGrid, par, opts, inspectionDist, label)

MList = MGrid(:);
n = numel(MList);

RIList = nan(n,1);
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
    progress = makeProgressPrinterPolicy2( ...
        n,label,opts.progressEverySeconds);

    afterEach(dq,@(~) progress());

    parfor k = 1:n
        try
            rep = evaluatePolicy2( ...
                MList(k),par,opts,inspectionDist,false,false);

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
    fprintf('Running serial Policy 2 grid search.\n');

    t0 = tic;
    lastPrint = 0;

    for k = 1:n
        try
            rep = evaluatePolicy2( ...
                MList(k),par,opts,inspectionDist,false,false);

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
    RIList,MList,CostRate,Availability,Ecycle,MTTFnew,MTTFstat, ...
    Ppm,Ppr,Pcr,Down,Err, ...
    'VariableNames',{ ...
    'RI','M','CostRate','Availability','ExpectedCycleLength', ...
    'MTTF_fromNew','MTTF_stationary','ExpectedPM','ExpectedPR', ...
    'ExpectedCR','ExpectedDowntime','Error'});
end

function progress = makeProgressPrinterPolicy2( ...
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
%  Distribution of X_Y+W, where Y~Exp(lambda)
% ========================================================================
function dist = buildPolicy2InspectionDistribution(par,opts)
Nfft = 2^nextpow2(opts.NfftInspection);
Nlinear = 2*Nfft;

fftUpperBound = max( ...
    opts.inspectionUpperFactor*par.L, ...
    50);

du = fftUpperBound/Nfft;
edges = (0:Nfft)*du;
centers = ((1:Nfft)-0.5)*du;

shockMass = diff(gammaCdfRate( ...
    edges,par.alphaS,par.betaS));

shockMass = max(real(shockMass),0);
shockTransform = fft(shockMass,Nlinear);

[laguerreX,laguerreW] = gaussLaguerre(opts.Nlag);

mixtureMass = zeros(1,Nfft);

for k = 1:numel(laguerreX)
    y = laguerreX(k)/par.lambda;

    if y <= 0
        continuousMass = zeros(1,Nfft);
        continuousMass(1) = 1;
    else
        continuousMass = diff(gammaCdfRate( ...
            edges,par.alpha*y,par.beta));

        continuousMass = max(real(continuousMass),0);
    end

    convolutionMass = real(ifft( ...
        fft(continuousMass,Nlinear).*shockTransform));

    convolutionMass = convolutionMass(1:Nfft);
    convolutionMass(abs(convolutionMass) < 1e-15) = 0;
    convolutionMass(convolutionMass < 0) = 0;

    mixtureMass = mixtureMass + ...
        laguerreW(k)*convolutionMass;
end

uGrid = linspace(0,par.L,opts.NuInspection);

uBinIndex = floor(uGrid/du)+1;
uBinIndex = min(max(uBinIndex,1),Nfft);

cdfAtEdges = [0,cumsum(mixtureMass)];
pdfAtCenters = mixtureMass/du;

dist.uGrid = uGrid;
dist.uBinIndex = uBinIndex;

dist.cdf = interp1( ...
    edges,cdfAtEdges,uGrid,'linear',0);

% Piecewise-constant finite-volume density.
dist.pdf = pdfAtCenters(uBinIndex);

dist.cdf = min(max(dist.cdf,0),1);
dist.pdf = max(dist.pdf,0);

dist.laguerreX = laguerreX;
dist.laguerreW = laguerreW;

dist.Nfft = Nfft;
dist.Nlinear = Nlinear;
dist.fftEdges = edges;
dist.fftCenters = centers;
dist.shockTransform = shockTransform;

dist.fftUpperBound = fftUpperBound;
dist.fftStep = du;
dist.totalRetainedMass = sum(mixtureMass);

nTail = max(10,round(0.02*Nfft));
dist.upperTailMass = ...
    sum(mixtureMass(end-nTail+1:end));
end

function F = inspectionCdfRate(x,dist)
F = zeros(size(x));

F(x <= 0) = 0;

valid = x > 0 & x <= dist.uGrid(end);

if any(valid(:))
    F(valid) = interp1( ...
        dist.uGrid,dist.cdf,x(valid),'linear');
end

F(x > dist.uGrid(end)) = dist.cdf(end);
F = min(max(real(F),0),1);
end

function pdf = inspectionPdfRate(x,dist)
pdf = zeros(size(x));

valid = x > 0 & x <= dist.uGrid(end);

if any(valid(:))
    pdf(valid) = interp1( ...
        dist.uGrid,dist.pdf,x(valid),'linear',0);
end

pdf = max(real(pdf),0);
end

%% ========================================================================
%  Numerical building blocks
% ========================================================================
function pdf = gammaPdfRate(x,a,beta)
[xb,ab] = commonSize(x,a);
pdf = zeros(size(xb));

idx = xb > 0 & ab > 0;

if any(idx(:))
    aa = ab(idx);
    xx = xb(idx);

    logp = aa.*log(beta) - gammaln(aa) + ...
        (aa-1).*log(xx) - beta.*xx;

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
        (p-1).*log(xx) + ...
        (q-1).*log1p(-xx) - ...
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

function [x,w] = gaussLaguerre(n)
diagonal = (2*(1:n)-1).';
offDiagonal = (1:n-1).';

matrix = ...
    diag(diagonal) + ...
    diag(offDiagonal,1) + ...
    diag(offDiagonal,-1);

[V,D] = eig(matrix);

x = diag(D);
[x,index] = sort(x);
V = V(:,index);

w = (V(1,:).^2).';
end

%% ========================================================================
%  Reporting and plotting
% ========================================================================
function printPolicy2Report(rep,titleStr)
fprintf('\n------------------------------------------------------------\n');
fprintf('%s\n',titleStr);
fprintf('------------------------------------------------------------\n');

fprintf('RI = --, M = %.8g, M/L = %.4f\n', ...
    rep.policy.M,rep.policy.M/rep.parameters.L);

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
    fprintf('  omitted kernel tail bound  = %.3e\n', ...
        rep.diagnostics.reliabilityKernel.omittedFirstShockTailBound);
    fprintf('  max conditional FFT tail   = %.3e\n', ...
        rep.diagnostics.reliabilityKernel.maximumConditionalUpperTailMass);
end
end

function checkBoundaryPolicy2(bestM,MGrid,label)
if abs(bestM-min(MGrid)) < 1e-12 || ...
        abs(bestM-max(MGrid)) < 1e-12

    fprintf(['\nWarning: %s-grid optimum is on the search boundary. ', ...
        'Consider enlarging the grid range.\n'],label);
end
end

function plotPolicy2CostCurve(tab,par,repBest)
valid = isfinite(tab.CostRate);

if ~any(valid)
    warning('No valid Policy 2 results are available for plotting.');
    return;
end

tab = tab(valid,:);
[~,order] = sort(tab.M);
tab = tab(order,:);

fig = figure( ...
    'Name','Policy 2 coarse-grid cost curve', ...
    'Color','w', ...
    'Position',[100,100,760,520]);

ax = axes('Parent',fig);

plot(ax, ...
    tab.M/par.L, ...
    tab.CostRate, ...
    'LineWidth',1.5);

hold(ax,'on');

plot(ax, ...
    repBest.policy.M/par.L, ...
    repBest.costRate, ...
    'kp', ...
    'MarkerFaceColor','y', ...
    'MarkerSize',12, ...
    'LineWidth',1.2);

xlabel(ax,'Preventive maintenance threshold M/L');
ylabel(ax,'Long-run cost rate');
title(ax,'Policy 2 long-run cost rate');

grid(ax,'on');
box(ax,'on');

savefig(fig,'GaAs_policy2_cost_curve.fig');

try
    exportgraphics( ...
        fig,'GaAs_policy2_cost_curve.png','Resolution',300);
catch
    print(fig,'GaAs_policy2_cost_curve.png','-dpng','-r300');
end

fprintf(['Saved Policy 2 cost curve to ', ...
    'GaAs_policy2_cost_curve.fig and GaAs_policy2_cost_curve.png\n']);
end

function plotPolicy2ReliabilityCurve(repBest)
if ~isfield(repBest,'reliability') || ...
        isempty(repBest.reliability.t)

    fprintf('Policy 2 reliability curve is empty; plot skipped.\n');
    return;
end

fig = figure( ...
    'Name','Policy 2 reliability function', ...
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

title(ax,'Reliability function under optimal Policy 2');

savefig(fig,'GaAs_policy2_reliability.fig');

try
    exportgraphics( ...
        fig,'GaAs_policy2_reliability.png','Resolution',300);
catch
    print(fig,'GaAs_policy2_reliability.png','-dpng','-r300');
end

fprintf(['Saved Policy 2 reliability plot to ', ...
    'GaAs_policy2_reliability.fig and GaAs_policy2_reliability.png\n']);
end

function validatePolicy2(M,par)
if ~(M > 0 && M < par.L)
    error('M must satisfy 0 < M < L.');
end

if ~(par.lambda > 0)
    error('Policy 2 requires lambda > 0.');
end
end
