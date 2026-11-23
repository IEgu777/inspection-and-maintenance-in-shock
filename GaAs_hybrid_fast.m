% GaAs LASER CASE (Example 2) -- HYBRID
% Adapted directly from run_shock_cbm_case_v6_fast.m; local numerical functions are retained.
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

% GAAS_HYBRID_FAST (adapted from run_shock_cbm_case_v6_fast.m)
% Faster script version for the hybrid inspection-CBM policy in a shock environment.
% Fast version changes relative to V5:
%   1) Reuse the Beta post-maintenance density matrix within each policy.
%   2) Replace costly hypergeom/convolution quadrature by a convergent,
%      nonnegative Gamma-mixture expansion for sums of unequal-rate Gammas.
%   3) Use fzero instead of 80 fixed bisection steps for inspection time.
%   4) Skip MTTF during cost-only policy screening; compute it for the final policy.
%   5) Search an admissible interior grid (quick or full mode).
% The numerical model and quadrature sizes are retained; parameters and costs
% are adapted to the GaAs case. Validate the outputs in MATLAB.
%
% Main changes relative to the uploaded V4 version:
%   1) Only coarse-grid optimization is performed. The fine-grid stage is removed.
%   2) The final policy evaluation is performed directly at the coarse-grid optimum.
%   3) This file is a script, not a main function. Local functions are kept below.
%   4) Runtime is recorded for every local function called directly by the main script.
%   5) The cost contour is drawn from the coarse grid only. Contour levels are denser
%      near the low-cost region, but only selected contour lines are labelled.
%   6) Unused plotting functions from the previous version are removed.
%   7) Corrected component-wise transition-mass normalization and reliability recursion.
%      Model parameters, costs and search ranges are adapted for the GaAs case.
%
% How to run:
%   >> GaAs_hybrid_fast
%
% Main output variables in the MATLAB workspace:
%   results       - struct containing parameters, options, grid table, best policy, timing log
%   coarseTab     - coarse-grid search results
%   repBest       - full evaluation of the coarse-grid optimum
%   timingLog     - runtime table for main-script function calls

clearvars;
clc;
close all;

% Switch to 'full' before reporting final numerical results in the paper.
runMode = getenv('GAAS_POLICY_RUN_MODE');
if isempty(runMode)
    runMode = 'full'; % Default to the paper grid when run individually.
end
assert(any(strcmpi(runMode, {'quick','full'})), ...
    'runMode must be quick or full.');

scriptTimer = tic;
timingNames = {};
timingSeconds = [];

fprintf('\n============================================================\n');
fprintf(' GaAs laser: hybrid inspection-CBM (250-hour unit)\n');
fprintf('============================================================\n');

%% --------------------- Model and cost parameters -----------------------
% Shared by all four policy implementations; change GaAs_case_parameters.m.
par = GaAs_case_parameters();

%% --------------------- Numerical options --------------------
opts.Nz       = 60;     % grid points on post-maintenance state interval (0,M)
opts.Nv       = 80;     % quadrature points for v in [M,L]
opts.Ny       = 18;     % quadrature points for shock-arrival time y in [0,tau]
opts.NtInt    = 80;     % quadrature points for time integrals in d(s) and MTTF source term
opts.Nconv    = 60;     % quadrature points for survival of X+W
opts.NuRel    = 18;     % quadrature points in reliability recursion
opts.NtRel    = 60;     % reliability time-grid size
opts.tRelMaxFactor = 2.5;
opts.computeReliabilityCurve = false;  % true: generate curve only at final optimum
opts.useParallel = true;
opts.showProgress = true;
opts.progressEverySeconds = 3;
opts.minMassWarning = 5e-2;
opts.warnOnMass = true;  % warnings enabled for the final full-policy evaluation

%% --------------------- Coarse grid search only --------------
% Exclude invalid boundary values, where the model requires 0<RI<1, 0<M<L.
if strcmpi(runMode, 'quick')
    coarse.RIGrid = [0.05:0.05:0.95, 0.96:0.01:0.99];
    coarse.MGrid  = (0.05:0.05:0.95) * par.L;
else
    coarse.RIGrid = [0.01:0.01:0.98, 0.981:0.001:0.999];
    coarse.MGrid  = (0.01:0.01:0.99) * par.L;
end

if strcmpi(runMode, 'full')
    assert(numel(coarse.RIGrid)==117 && numel(coarse.MGrid)==99, ...
        'Incorrect full search-grid size: expected 117 RI x 99 M.');
    assert(abs(coarse.RIGrid(end)-0.999)<1e-12 && ...
        any(abs(coarse.RIGrid-0.995)<1e-12), ...
        'Dense RI grid must contain 0.995 and end at 0.999.');
    fprintf('RI DENSE GRID ACTIVE: FULL | 117 RI x 99 M = 11583 candidates\n');
end

fprintf('\n[1] Coarse grid search: %d x %d = %d policies\n', ...
    numel(coarse.RIGrid), numel(coarse.MGrid), numel(coarse.RIGrid)*numel(coarse.MGrid));

functionTimer = tic;
coarseTab = gridSearchPolicy(coarse.RIGrid, coarse.MGrid, par, opts, 'coarse');
timingNames{end+1,1} = 'gridSearchPolicy';
timingSeconds(end+1,1) = toc(functionTimer);

validPolicyIdx = find(isfinite(coarseTab.CostRate));
if isempty(validPolicyIdx)
    error('Every grid policy evaluation failed. Inspect coarseTab.Error.');
end
[~, idxInValid] = min(coarseTab.CostRate(validPolicyIdx));
idxBest = validPolicyIdx(idxInValid);
bestPolicy.RI = coarseTab.RI(idxBest);
bestPolicy.M  = coarseTab.M(idxBest);
bestPolicy.CostRate = coarseTab.CostRate(idxBest);
%%要改true
fprintf('\n[2] Re-evaluating coarse-grid optimum with full MTTF metrics...\n');
functionTimer = tic;
% Full evaluation (including MTTF) is performed only for the best policy.
repBest = evaluatePolicy(bestPolicy.RI, bestPolicy.M, par, opts, opts.computeReliabilityCurve, true);
timingNames{end+1,1} = 'evaluatePolicy at coarse optimum';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
printPolicyReport(repBest, 'Final policy from coarse-grid optimization');
if ~isfinite(repBest.MTTF_fromNew) || ~isfinite(repBest.MTTF_stationary) || ...
        ~isfinite(repBest.averageAvailability)
    warning(['Hybrid final metrics contain NaN/Inf. Do not report them. ', ...
        'Check transition-mass diagnostics and grid convergence.']);
end
timingNames{end+1,1} = 'printPolicyReport';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
checkBoundary(bestPolicy.RI, coarse.RIGrid, bestPolicy.M, coarse.MGrid, 'coarse');
timingNames{end+1,1} = 'checkBoundary';
timingSeconds(end+1,1) = toc(functionTimer);

% To compute a final reliability curve, enable the optional flag in the
% final evaluation call (a computationally costly calculation).
%% --------------------- Plots --------------------------------
functionTimer = tic;
try
    plotCostContourCoarseStyle(coarseTab, par, repBest);
catch ME
    warning('Cost contour plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotCostContourCoarseStyle';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
try
    plotFinalReliabilityCurve(repBest);
catch ME
    warning('Reliability plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotFinalReliabilityCurve';
timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Results and timing --------------------
timingNames{end+1,1} = 'total script time before save';
timingSeconds(end+1,1) = toc(scriptTimer);
timingLog = table(timingNames, timingSeconds, 'VariableNames', {'FunctionOrStep','Seconds'});

fprintf('\n---------------- Runtime log for main-script function calls ----------------\n');
disp(timingLog);

results = struct();
results.caseName = 'GaAs laser, Example 2';
results.runMode = runMode;
results.timeUnitHours = 250;
results.costUnit = 'dimensionless relative cost';
results.shockArrivalNote = 'lambda is estimated grid-event frequency; HPP frequency-matching approximation.';
results.parameters      = par;
results.options         = opts;
results.coarseGrid      = coarse;
results.coarseGridTable = coarseTab;
results.bestPolicy      = bestPolicy;
results.finalBest       = repBest;
results.timingLog       = timingLog;

save('GaAs_hybrid_results.mat', 'results');
fprintf('\nSaved results to GaAs_hybrid_results.mat\n');
fprintf('Done.\n');

%% ========================================================================
%  Main policy evaluator
% ========================================================================
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
function printPolicyReport(rep, titleStr)
fprintf('\n------------------------------------------------------------\n');
fprintf('%s\n', titleStr);
fprintf('------------------------------------------------------------\n');
fprintf('RI = %.8f, M = %.8g, M/L = %.4f\n', rep.policy.RI, rep.policy.M, rep.policy.M / rep.parameters.L);
fprintf('Cost rate C_inf              = %.10g\n', rep.costRate);
fprintf('Expected cycle length        = %.10g\n', rep.expectedCycleLength);
fprintf('Expected downtime/cycle      = %.10g\n', rep.expectedDowntimePerCycle);
fprintf('Average availability         = %.10f\n', rep.averageAvailability);
fprintf('MTTF, stationary initial pi  = %.10g\n', rep.MTTF_stationary);
fprintf('MTTF, new state s=0          = %.10g\n', rep.MTTF_fromNew);
fprintf('Stationary mass at state 0   = %.10f\n', rep.stationary.a0);

fprintf('\nExpected event counts per cycle:\n');
fprintf('  inspection                 = %.10g\n', rep.expectedCountsPerCycle.inspection);
fprintf('  PM                         = %.10g\n', rep.expectedCountsPerCycle.PM);
fprintf('  preventive replacement     = %.10g\n', rep.expectedCountsPerCycle.preventiveReplacement);
fprintf('  corrective replacement     = %.10g\n', rep.expectedCountsPerCycle.correctiveReplacement);

fprintf('\nCost-rate components and shares:\n');
fprintf('  inspection                 = %-12.6g  share = %6.2f%%\n', rep.costRateComponents.inspection, 100*rep.costShare.inspection);
fprintf('  PM                         = %-12.6g  share = %6.2f%%\n', rep.costRateComponents.PM, 100*rep.costShare.PM);
fprintf('  preventive replacement     = %-12.6g  share = %6.2f%%\n', rep.costRateComponents.preventiveReplacement, 100*rep.costShare.preventiveReplacement);
fprintf('  corrective replacement     = %-12.6g  share = %6.2f%%\n', rep.costRateComponents.correctiveReplacement, 100*rep.costShare.correctiveReplacement);
fprintf('  downtime                   = %-12.6g  share = %6.2f%%\n', rep.costRateComponents.downtime, 100*rep.costShare.downtime);

fprintf('\nDiagnostics:\n');
fprintf('  max transition mass error  = %.3e\n', rep.diagnostics.maxFullTransitionMassError);
fprintf('  rcond(I-P_survival)        = %.3e\n', rep.diagnostics.minMTTFLinearMatrixReciprocalCondition);

if ~isempty(rep.reliability.t)
    fprintf('\nReliability curve generated with %d time points.\n', numel(rep.reliability.t));
    fprintf('  R_new(t_end)               = %.6g\n', rep.reliability.RfromNew(end));
    fprintf('  R_stationary(t_end)        = %.6g\n', rep.reliability.Rstationary(end));
end
end

function checkBoundary(bestRI, RIGrid, bestM, MGrid, label)
if abs(bestRI - min(RIGrid)) < 1e-12 || abs(bestRI - max(RIGrid)) < 1e-12 || ...
   abs(bestM  - min(MGrid))  < 1e-12 || abs(bestM  - max(MGrid))  < 1e-12
    fprintf('\nWarning: %s-grid optimum is on the search boundary. Consider enlarging the grid range.\n', label);
end
end


function plotCostContourCoarseStyle(tab, par, repBest)
% Plot long-run cost-rate contours from the coarse grid only.
% Low-cost levels are intentionally denser. Only selected contour lines are
% labelled to avoid label crowding.

valid = isfinite(tab.CostRate);
if ~any(valid)
    warning('No valid grid-search results are available for contour plotting.');
    return;
end

tab = tab(valid,:);
RIvals = unique(tab.RI);
Mvals  = unique(tab.M);
xVals  = Mvals ./ par.L;

Z = nan(numel(RIvals), numel(Mvals));
for i = 1:numel(RIvals)
    for j = 1:numel(Mvals)
        idx = abs(tab.RI - RIvals(i)) < 1e-12 & abs(tab.M - Mvals(j)) < 1e-12;
        if any(idx)
            Z(i,j) = tab.CostRate(find(idx,1));
        end
    end
end

finiteZ = Z(isfinite(Z));
if isempty(finiteZ)
    warning('Cost values are not suitable for contour plotting.');
    return;
end
zMin = min(finiteZ);
zMax = max(finiteZ);
if ~(isfinite(zMin) && isfinite(zMax)) || abs(zMax-zMin) < eps
    warning('Cost values are not suitable for contour plotting.');
    return;
end

% Denser levels near the low-cost region: u.^2.2 clusters levels around zMin.
nLevels = 26;
u = linspace(0, 1, nLevels).^2.2;
levels = unique(zMin + (zMax - zMin) * u);

fig = figure('Name','Coarse-grid cost contour', ...
             'Color','w', ...
             'Position',[100, 100, 780, 560]);
ax = axes('Parent', fig);
hold(ax, 'on');

[C,h] = contour(ax, xVals, RIvals, Z, levels, 'LineWidth', 0.90);
colormap(ax, jet(256));
cb = colorbar(ax);
ylabel(cb, 'Long-run cost rate');
caxis(ax, [zMin, zMax]);

% Label only a subset of contour levels. This avoids dense overlapping labels.
labelLevels = levels(3:4:end);
if ~isempty(labelLevels)
    clabel(C, h, labelLevels, 'FontSize', 8, 'Color', 'k', 'LabelSpacing', 520);
end

% Mark the optimum explicitly, since the minimum contour can degenerate to a point.
xBest = repBest.policy.M / par.L;
yBest = repBest.policy.RI;
zBest = repBest.costRate;
plot(ax, xBest, yBest, 'kp', ...
    'MarkerFaceColor', 'y', ...
    'MarkerSize', 12, ...
    'LineWidth', 1.2);
text(ax, xBest, yBest, sprintf('  optimum %.4g', zBest), ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'Color', 'k', ...
    'VerticalAlignment', 'bottom');

xlabel(ax, 'Preventive maintenance threshold M/L');
ylabel(ax, 'Reliability threshold R_I');
title(ax, 'Long-run cost rate contour');
box(ax, 'on');
set(ax, 'Layer', 'top');

savefig(fig, 'GaAs_hybrid_cost_contour.fig');
try
    exportgraphics(fig, 'GaAs_hybrid_cost_contour.png', 'Resolution', 300);
catch
    print(fig, 'GaAs_hybrid_cost_contour.png', '-dpng', '-r300');
end
fprintf('Saved coarse contour plot to GaAs_hybrid_cost_contour.fig and GaAs_hybrid_cost_contour.png\n');
end

function plotFinalReliabilityCurve(repBest)
% Plot the two reliability curves only for the final best policy.
if ~isfield(repBest, 'reliability') || isempty(repBest.reliability.t)
    fprintf('Final reliability curve is empty; reliability plot skipped.\n');
    return;
end

fig = figure('Name','Reliability function under final policy', 'Color', 'w');
ax = axes('Parent', fig);
plot(ax, repBest.reliability.t, repBest.reliability.Rstationary, 'LineWidth', 1.5); hold(ax, 'on');
plot(ax, repBest.reliability.t, repBest.reliability.RfromNew, '--', 'LineWidth', 1.5);
grid(ax, 'on');
xlabel(ax, 'Time');
ylabel(ax, 'Reliability');
legend(ax, 'stationary initial distribution', 'new state s=0', 'Location', 'best');
title(ax, 'Reliability function under final best policy');

savefig(fig, 'GaAs_hybrid_reliability.fig');
try
    exportgraphics(fig, 'GaAs_hybrid_reliability.png', 'Resolution', 300);
catch
    print(fig, 'GaAs_hybrid_reliability.png', '-dpng', '-r300');
end
fprintf('Saved final reliability plot to GaAs_hybrid_reliability.fig and GaAs_hybrid_reliability.png\n');
end

function validatePolicy(RI, M, par)
if ~(RI > 0 && RI < 1)
    error('RI must satisfy 0 < RI < 1.');
end
if ~(M > 0 && M < par.L)
    error('M must satisfy 0 < M < L.');
end
end
