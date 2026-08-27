%RUN_POLICY1_CASE_V2_COARSE_TIMED
% Script version for Policy 1 in Section 6.4.
%
% Policy 1:
%   Inspections are conducted only when the system reliability reaches the
%   predefined threshold RI. Shock arrivals do not trigger inspections.
%
% The reliability-triggered interval is based on the compound degradation
%
%   Delta S_t = X_t + sum_{j=1}^{N(t)} W_j,
%
% as defined in Appendix B of the manuscript.
%
% Main features:
%   1) Only coarse-grid optimization is performed.
%   2) The final policy is evaluated directly at the coarse-grid optimum.
%   3) This file is a script, not a main function.
%   4) Runtime is recorded for every local function called by the script.
%   5) The reliability curve starts from the stationary post-maintenance
%      distribution of Policy 1.
%
% V2 corrections:
%   - do-nothing and successful-PM transition densities are corrected
%     separately to their theoretical probability masses;
%   - RAW transition errors are retained as diagnostics;
%   - the monotone CDF copy is used only to solve the inspection time,
%     while transition probabilities use the CDF/PDF from the same FFT mass;
%   - the reliability horizon is based on the stationary MTTF.
%
% How to run:
%   >> run_policy1_case_v2_coarse_timed
%
% Main output variables:
%   results       - complete Policy 1 result structure
%   coarseTab     - coarse-grid search table
%   repBest       - final evaluation at the coarse-grid optimum
%   timingLog     - runtime table

clearvars;
clc;
close all;

scriptTimer = tic;
timingNames = {};
timingSeconds = [];

fprintf('\n============================================================\n');
fprintf(' Policy 1: reliability-triggered inspections only -- V2 corrected\n');
fprintf('============================================================\n');

%% --------------------- Model parameters ---------------------
% These values are identical to those used by the proposed policy code.
par.alpha  = 1.02;       % Gamma-process shape rate alpha
par.beta   = 1.20;       % Gamma-process rate beta
par.lambda = 0.25;       % HPP shock arrival intensity

% Shock magnitude W ~ Gamma(alphaS,betaS), betaS is RATE.
par.alphaS = 2.5;
par.betaS  = 1.0;

% Imperfect PM ratio follows Beta(p,q).
par.p = 6;
par.q = 5;

% Failure threshold.
par.L = 12.5;

%% --------------------- Cost parameters ----------------------
par.Ci = 5;     % inspection cost
par.Cm = 100;    % imperfect PM cost
par.Cp = 280;    % preventive replacement cost
par.Cc = 350;    % corrective replacement cost
par.Cd = 50;     % downtime cost rate

%% --------------------- Numerical options --------------------
opts.Nz       = 60;     % grid points on post-maintenance interval (0,M)
opts.Nv       = 80;     % quadrature points for v in [M,L]
opts.NtInt    = 100;    % quadrature points for downtime and MTTF source

% Reliability recursion. The time horizon is based on the stationary MTTF.
opts.NtRel          = 301;   % minimum reliability time-grid size
opts.maxNtRel       = 4000;  % safeguard for automatic time-grid refinement
opts.dtRelMaximum   = 0.20;  % desired largest time step
opts.tRelMaxFactor  = 6.0;   % tMax >= factor*max(MTTFs,Ecycle)
opts.reliabilityTailWarning = 1e-5;

% Compound-degradation table used for Policy 1.
opts.NtCompound   = 1201;
opts.NuCompound   = 1201;
opts.NfftCompound = 8192;
opts.tCompoundMax = 60;

opts.useParallel = true;
opts.showProgress = true;
opts.progressEverySeconds = 3;

% This threshold is applied to RAW, pre-correction quadrature errors.
opts.minMassWarning = 1e-3;

%% --------------------- Precompute compound distribution -----
fprintf('\n[0] Precomputing compound-degradation PDF and CDF table...\n');

functionTimer = tic;
compound = buildCompoundDistributionTable(par, opts);
timingNames{end+1,1} = 'buildCompoundDistributionTable';
timingSeconds(end+1,1) = toc(functionTimer);

fprintf('Compound table completed. t in [0, %.4g], u in [0, %.4g].\n', ...
    compound.tGrid(end), compound.uGrid(end));
fprintf('Minimum retained probability at FFT upper bound = %.8f\n', ...
    min(compound.totalRetainedMass));
fprintf('Maximum probability in the upper 2%% of the FFT grid = %.3e\n', ...
    max(compound.upperTailMass));

%% --------------------- Coarse grid search only --------------
coarse.RIGrid = linspace(0.6, 0.99, 40);
coarse.MGrid  = linspace(0.6, 0.99, 40) * par.L;

fprintf('\n[1] Policy 1 coarse grid search: %d x %d = %d policies\n', ...
    numel(coarse.RIGrid), numel(coarse.MGrid), ...
    numel(coarse.RIGrid)*numel(coarse.MGrid));

functionTimer = tic;
coarseTab = gridSearchPolicy1( ...
    coarse.RIGrid, coarse.MGrid, par, opts, compound, 'policy1 coarse');
timingNames{end+1,1} = 'gridSearchPolicy1';
timingSeconds(end+1,1) = toc(functionTimer);

validCost = isfinite(coarseTab.CostRate);
if ~any(validCost)
    error('No valid Policy 1 candidate was obtained. Inspect coarseTab.Error.');
end

validIndices = find(validCost);
[~, localBest] = min(coarseTab.CostRate(validCost));
idxBest = validIndices(localBest);

bestPolicy.RI = coarseTab.RI(idxBest);
bestPolicy.M  = coarseTab.M(idxBest);
bestPolicy.CostRate = coarseTab.CostRate(idxBest);

%%改true
fprintf('\n[2] Re-evaluating the Policy 1 optimum with reliability curve...\n');

functionTimer = tic;
repBest = evaluatePolicy1( ...
    bestPolicy.RI, bestPolicy.M, par, opts, compound, false);
timingNames{end+1,1} = 'evaluatePolicy1 at coarse optimum';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
printPolicy1Report(repBest, 'Final Policy 1 from coarse-grid optimization');
timingNames{end+1,1} = 'printPolicy1Report';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
checkBoundaryPolicy1( ...
    bestPolicy.RI, coarse.RIGrid, ...
    bestPolicy.M, coarse.MGrid, 'Policy 1 coarse');
timingNames{end+1,1} = 'checkBoundaryPolicy1';
timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Plots --------------------------------
functionTimer = tic;
try
    plotPolicy1CostContour(coarseTab, par, repBest);
catch ME
    warning('Policy 1 cost-contour plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotPolicy1CostContour';
timingSeconds(end+1,1) = toc(functionTimer);

functionTimer = tic;
try
    plotPolicy1ReliabilityCurve(repBest);
catch ME
    warning('Policy 1 reliability plotting failed: %s', ME.message);
end
timingNames{end+1,1} = 'plotPolicy1ReliabilityCurve';
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
results.policyName      = 'Policy 1';
results.parameters      = par;
results.options         = opts;
results.coarseGrid      = coarse;
results.coarseGridTable = coarseTab;
results.bestPolicy      = bestPolicy;
results.finalBest       = repBest;
results.timingLog       = timingLog;
results.compoundDiagnostics = struct( ...
    'tMaximum', compound.tGrid(end), ...
    'uMaximumForPolicyCalculations', compound.uGrid(end), ...
    'fftUpperBound', compound.fftUpperBound, ...
    'minimumRetainedMass', min(compound.totalRetainedMass), ...
    'maximumUpperTailMass', max(compound.upperTailMass));

save('policy1_results_v2_coarse_timed.mat', 'results');
fprintf('\nSaved Policy 1 results to policy1_results_v2_coarse_timed.mat\n');
fprintf('Done.\n');

%% ========================================================================
%  Main Policy 1 evaluator
% ========================================================================
function rep = evaluatePolicy1(RI, M, par, opts, compound, doReliabilityCurve)
if nargin < 6
    doReliabilityCurve = false;
end

validatePolicy1(RI, M, par);

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
Kdens = zeros(Ns, Nz);
Q1full = zeros(Ns,1);
Psurv = zeros(Ns,Ns);

% Diagnostics are recorded BEFORE and AFTER component-wise correction.
PrawMassError = zeros(Ns,1);
PcorrectedMassError = zeros(Ns,1);
PdoNothingMassError = zeros(Ns,1);
PpmSuccessMassError = zeros(Ns,1);
PpmEventMassError = zeros(Ns,1);

cond = initConditionalStructPolicy1(Ns);

for i = 1:Ns
    row = buildPolicy1OneStepRow( ...
        states(i), RI, M, par, opts, compound, ...
        z, wz, v, wv, pmSurvVec, pmPdfMat);

    Kdens(i,:) = row.fullInteriorDensity;
    Q1full(i)  = row.fullMass0;
    Psurv(i,:) = row.PsurvRow;

    PrawMassError(i) = row.rawMassError;
    PcorrectedMassError(i) = row.correctedMassError;
    PdoNothingMassError(i) = row.doNothingMassError;
    PpmSuccessMassError(i) = row.pmSuccessMassError;
    PpmEventMassError(i) = row.pmEventMassError;

    cond.tau(i) = row.tau;
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
    warning('Policy 1 stationary density H contains noticeable negative values.');
end
H = max(H,0);

a0 = 1 / (1 + sum(H .* wz));
piVec = [a0; a0 * H .* wz];
piVec = piVec / sum(piVec);

% Stationary expected cycle quantities.
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

% MTTF:
% m = source + P_survival*m.
I = eye(Ns);
mVec = (I-Psurv) \ cond.mttfSource;

MTTF_stationary = sum(piVec .* mVec);
MTTF_fromNew    = mVec(1);

% Reliability curve starting from each state and from the stationary law.
rel = struct('t', [], 'Rstationary', [], 'RfromNew', [], 'RbyState', []);

if doReliabilityCurve
    tMax = opts.tRelMaxFactor * ...
        max([MTTF_stationary, MTTF_fromNew, Ecycle]);

    if ~isfinite(tMax) || tMax <= 0
        tMax = 1e3;
    end

    tMax = min(tMax,1e6);

    positiveTau = cond.tau(cond.tau > 0);
    NtRel = max(2,round(opts.NtRel));

    % Resolve both the delay tau(s) and the requested absolute time step.
    NtRequiredByStep = ceil(tMax/max(opts.dtRelMaximum,eps))+1;
    NtRel = max(NtRel,NtRequiredByStep);

    if ~isempty(positiveTau)
        desiredStep = 0.5*min(positiveTau);
        NtRequiredByDelay = ceil(tMax/max(desiredStep,eps))+1;
        NtRel = max(NtRel,NtRequiredByDelay);
    end

    if NtRel > opts.maxNtRel
        warning(['Policy 1 reliability grid was capped at %d points. ', ...
            'Increase opts.maxNtRel for a finer delay recursion.'], ...
            opts.maxNtRel);
        NtRel = opts.maxNtRel;
    end

    tGrid = linspace(0,tMax,NtRel);

    Rstate = reliabilityCurvePolicy1( ...
        tGrid, states, cond.tau, Psurv, par, compound);

    rel.t = tGrid(:);
    rel.RbyState = Rstate;
    rel.Rstationary = enforceReliabilityVector(Rstate*piVec);
    rel.RfromNew = enforceReliabilityVector(Rstate(:,1));
end

rep = struct();
rep.policy.RI = RI;
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
    'minMTTFLinearMatrixReciprocalCondition', rcond(I-Psurv));

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
    warning(['Policy 1 maximum uncorrected component error is %.3g. ', ...
        'Increase opts.Nz, opts.Nv, or the compound grids.'], ...
        rep.diagnostics.maxUncorrectedComponentError);
end

if doReliabilityCurve && ...
        rep.diagnostics.reliabilityTail > opts.reliabilityTailWarning
    warning(['Policy 1 reliability at the final time is %.3g. ', ...
        'Increase opts.tRelMaxFactor or opts.maxNtRel.'], ...
        rep.diagnostics.reliabilityTail);
end
end

%% ========================================================================
%  Policy 1 one-step kernel, Appendix B
% ========================================================================
function row = buildPolicy1OneStepRow( ...
    s, RI, M, par, opts, compound, z, wz, v, wv, ...
    pmSurvVec, pmPdfMat)

L = par.L;

tau = tauFromCompoundReliability(s, RI, par, compound);

% Compound degradation density at the inspection time.
g_v = compoundPdfRate(v-s, tau, compound);
g_v(v <= s) = 0;

g_z = compoundPdfRate(z-s, tau, compound);
g_z(z <= s) = 0;

% The RAW CDF and PDF are generated from the same FFT probability masses.
F_M = compoundCdfRate(M-s, tau, compound);
F_L = compoundCdfRate(L-s, tau, compound);

F_M = min(max(F_M,0),1);
F_L = min(max(F_L,F_M),1);

Pcr = 1-F_L;
Ppm = F_L-F_M;

% The PM event is split into successful PM and preventive replacement.
% Their ratio is first computed from the same v-quadrature and then the
% complete PM event is scaled to the CDF-based probability Ppm.
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
    % Conservative fallback for an unresolved, extremely small PM event.
    Ppr = Ppm;
    pmSuccessTarget = 0;
end

% Remove only roundoff and enforce exact PM-event balance.
Ppr = min(max(Ppr,0),Ppm);
pmSuccessTarget = max(Ppm-Ppr,0);

pmSuccessDensityRaw = pmPdfMat * (wv .* g_v);
doNothingDensityRaw = g_z;

% The two continuous transition components have different theoretical
% masses and are corrected separately:
%
%   integral Q3(s,z)dz = F1(M-s;tau),
%   integral Q2(s,z)dz = Ppm(s)-Ppr(s).
doNothingTarget = F_M;

doNothingMassRaw = sum(doNothingDensityRaw .* wz);
pmSuccessMassRaw = sum(pmSuccessDensityRaw .* wz);

doNothingMassError = doNothingMassRaw-doNothingTarget;

% This measures only the z-grid integration error in Q2. The separate
% pmEventMassError measures the v-grid error in the total PM event.
pmSuccessMassError = ...
    pmSuccessMassRaw-pmSuccessProbabilityRaw;

pmEventMassError = pmEventMassRaw-Ppm;

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
    % This fallback is used only when a very narrow interval is missed by
    % the midpoint state grid.
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
    % Shape-only fallback for an exceptionally small PM probability.
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

% Survival kernel excludes corrective replacement and retains preventive
% replacement at state 0.
PsurvRow = zeros(1,numel(z)+1);
PsurvRow(1) = Ppr;
PsurvRow(2:end) = (interiorDensity .* wz)';

% Under Policy 1, the cycle length is deterministic conditional on s.
mu = tau;

% b(s)=integral_0^tau F1(L-s;t)dt and
% d(s)=integral_0^tau G1(L-s;t)dt=tau-b(s).
mttfSource = compoundReliabilityIntegral( ...
    L-s, tau, compound, opts.NtInt);

downtime = max(tau-mttfSource,0);

row = struct();
row.tau = tau;
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


function cond = initConditionalStructPolicy1(Ns)
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
%  Policy 1 reliability recursion
% ========================================================================
function Rstate = reliabilityCurvePolicy1( ...
    tGrid, states, tauVec, Psurv, par, compound)

Ns = numel(states);
Nt = numel(tGrid);

if Nt < 2
    error('Policy 1 reliability recursion requires at least two time points.');
end

Rstate = zeros(Nt,Ns);
Rstate(1,:) = 1;

for it = 2:Nt
    t = tGrid(it);
    Rnow = zeros(1,Ns);

    for is = 1:Ns
        s = states(is);
        tau = tauVec(is);

        if t < tau
            Rnow(is) = compoundCdfRate( ...
                par.L-s, t, compound);
        else
            tRemain = t-tau;
            Rpast = interpolateReliability( ...
                tRemain, tGrid(1:it-1), Rstate(1:it-1,:));

            Rnow(is) = Psurv(is,:) * Rpast(:);
        end
    end

    Rstate(it,:) = min(max(real(Rnow),0),1);
end

% Reliability is theoretically nonincreasing. Remove very small numerical
% reversals caused by time interpolation.
Rstate = cummin(Rstate,1);
Rstate(1,:) = 1;
end

function Rv = interpolateReliability(tq, tPast, RPast)
Ns = size(RPast,2);

if tq <= 0
    Rv = ones(1,Ns);
    return;
end

if isempty(tPast) || isempty(RPast)
    Rv = ones(1,Ns);
    return;
end

tPast = tPast(:);

if numel(tPast) < 2
    Rv = RPast(end,:);
elseif tq <= tPast(1)
    Rv = RPast(1,:);
elseif tq >= tPast(end)
    Rv = RPast(end,:);
else
    Rv = zeros(1,Ns);
    for j = 1:Ns
        Rv(j) = interp1(tPast,RPast(:,j),tq,'linear');
    end
end

Rv = min(max(real(Rv),0),1);
end

function R = enforceReliabilityVector(R)
R = min(max(real(R(:)),0),1);
R = cummin(R);
R(1) = 1;
end

%% ========================================================================
%  Grid search with parallel computing and progress display
% ========================================================================
function tab = gridSearchPolicy1( ...
    RIGrid, MGrid, par, opts, compound, label)

[RIList,MList] = ndgrid(RIGrid(:),MGrid(:));
RIList = RIList(:);
MList = MList(:);
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
    progress = makeProgressPrinterPolicy1( ...
        n,label,opts.progressEverySeconds);

    afterEach(dq,@(~) progress());

    parfor k = 1:n
        try
            rep = evaluatePolicy1( ...
                RIList(k),MList(k),par,opts,compound,false);

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
    fprintf('Running serial Policy 1 grid search.\n');

    t0 = tic;
    lastPrint = 0;

    for k = 1:n
        try
            rep = evaluatePolicy1( ...
                RIList(k),MList(k),par,opts,compound,false);

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

function progress = makeProgressPrinterPolicy1( ...
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
    meanMaximum + 12*sdMaximum + 2*par.L, ...
    50]);

du = fftUpperBound/Nfft;
fftEdges = (0:Nfft)*du;
fftCenters = ((1:Nfft)-0.5)*du;

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

    % Piecewise-constant finite-volume density. This avoids losing the
    % probability in the first FFT bin when the Gamma shape is below one.
    pdfRaw(it,:) = pdfAtCenters(uBinIndex);

    cdfRaw(it,:) = min(max(cdfRaw(it,:),0),1);
    pdfRaw(it,:) = max(pdfRaw(it,:),0);
end

% Only the inspection-time root finder uses the monotone-in-time copy.
% Transition probabilities continue to use cdfRaw, which is generated
% from exactly the same FFT masses as pdfRaw.
cdfForTau = cdfRaw;

for iu = 1:Nu
    cdfForTau(:,iu) = cummin(cdfForTau(:,iu));
end

compound.tGrid = tGrid;
compound.uGrid = uGrid;
compound.cdfRaw = cdfRaw;
compound.cdfForTau = cdfForTau;
compound.pdfRaw = pdfRaw;

% Backward-compatible aliases used nowhere outside this script.
compound.cdf = cdfRaw;
compound.pdf = pdfRaw;

compound.fftUpperBound = fftUpperBound;
compound.fftStep = du;
compound.totalRetainedMass = retainedMass;
compound.upperTailMass = upperTailMass;
end

function tau = tauFromCompoundReliability( ...
    s,RI,par,compound)

margin = par.L-s;

if margin <= 0
    tau = 0;
    return;
end

Ftime = interp1( ...
    compound.uGrid,compound.cdfForTau.', ...
    margin,'linear').';

Ftime(1) = 1;
Ftime = cummin(Ftime);

index = find(Ftime <= RI,1,'first');

if isempty(index)
    error(['Cannot bracket the Policy 1 inspection time. ', ...
        'Increase opts.tCompoundMax.']);
end

if index == 1
    tau = compound.tGrid(1);
    return;
end

Fleft = Ftime(index-1);
Fright = Ftime(index);
tLeft = compound.tGrid(index-1);
tRight = compound.tGrid(index);

if abs(Fright-Fleft) < 1e-14
    tau = tRight;
else
    tau = tLeft + ...
        (RI-Fleft)*(tRight-tLeft)/(Fright-Fleft);
end

tau = min(max(tau,tLeft),tRight);
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
    margin,tau,compound,Nt)

if tau <= 0
    value = 0;
    return;
end

t = ((1:Nt)' - 0.5) / Nt * tau;
wt = tau/Nt;

F = compoundCdfRate( ...
    ones(size(t))*margin,t,compound);

value = sum(wt.*F);
value = min(max(value,0),tau);
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

%% ========================================================================
%  Reporting and plotting
% ========================================================================
function printPolicy1Report(rep,titleStr)
fprintf('\n------------------------------------------------------------\n');
fprintf('%s\n',titleStr);
fprintf('------------------------------------------------------------\n');

fprintf('RI = %.8f, M = %.8g, M/L = %.4f\n', ...
    rep.policy.RI,rep.policy.M,rep.policy.M/rep.parameters.L);

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

function checkBoundaryPolicy1( ...
    bestRI,RIGrid,bestM,MGrid,label)

if abs(bestRI-min(RIGrid)) < 1e-12 || ...
   abs(bestRI-max(RIGrid)) < 1e-12 || ...
   abs(bestM-min(MGrid)) < 1e-12 || ...
   abs(bestM-max(MGrid)) < 1e-12

    fprintf(['\nWarning: %s-grid optimum is on the search boundary. ', ...
        'Consider enlarging the grid range.\n'],label);
end
end

function plotPolicy1CostContour(tab,par,repBest)
valid = isfinite(tab.CostRate);

if ~any(valid)
    warning('No valid Policy 1 results are available for contour plotting.');
    return;
end

tab = tab(valid,:);

RIvals = unique(tab.RI);
Mvals = unique(tab.M);
xVals = Mvals/par.L;

Z = nan(numel(RIvals),numel(Mvals));

for i = 1:numel(RIvals)
    for j = 1:numel(Mvals)
        idx = abs(tab.RI-RIvals(i)) < 1e-12 & ...
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
    'Name','Policy 1 coarse-grid cost contour', ...
    'Color','w', ...
    'Position',[100,100,780,560]);

ax = axes('Parent',fig);
hold(ax,'on');

[C,h] = contour( ...
    ax,xVals,RIvals,Z,levels,'LineWidth',0.90);

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
yBest = repBest.policy.RI;

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
ylabel(ax,'Reliability threshold R_I');
title(ax,'Policy 1 long-run cost-rate contour');

box(ax,'on');
set(ax,'Layer','top');

savefig(fig,'policy1_cost_contour_v2.fig');

try
    exportgraphics( ...
        fig,'policy1_cost_contour_v2.png','Resolution',300);
catch
    print(fig,'policy1_cost_contour_v2.png','-dpng','-r300');
end

fprintf(['Saved Policy 1 contour plot to ', ...
    'policy1_cost_contour_v2.fig and policy1_cost_contour_v2.png\n']);
end

function plotPolicy1ReliabilityCurve(repBest)
if ~isfield(repBest,'reliability') || ...
        isempty(repBest.reliability.t)

    fprintf('Policy 1 reliability curve is empty; plot skipped.\n');
    return;
end

fig = figure( ...
    'Name','Policy 1 reliability function', ...
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

title(ax,'Reliability function under optimal Policy 1');

savefig(fig,'policy1_reliability_curve_v2.fig');

try
    exportgraphics( ...
        fig,'policy1_reliability_curve_v2.png','Resolution',300);
catch
    print(fig,'policy1_reliability_curve_v2.png','-dpng','-r300');
end

fprintf(['Saved Policy 1 reliability plot to ', ...
    'policy1_reliability_curve_v2.fig and policy1_reliability_curve_v2.png\n']);
end

function validatePolicy1(RI,M,par)
if ~(RI > 0 && RI < 1)
    error('RI must satisfy 0 < RI < 1.');
end

if ~(M > 0 && M < par.L)
    error('M must satisfy 0 < M < L.');
end
end
