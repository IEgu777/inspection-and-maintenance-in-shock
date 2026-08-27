%SENSITIVITY_CD_50PCT
% One-at-a-time sensitivity analysis for downtime cost rate.
%
% Only two perturbed values are optimized:
%       50% of baseline  (-50%)
%      150% of baseline  (+50%)
%
% Initial optimization grid:
%   RI  = 0.65:0.01:0.99
%   M/L = 0.50:0.01:0.80
%
% If the optimum hits a movable boundary, the corresponding range is
% expanded automatically and the optimization is repeated.
%
% beta and betaS are RATE parameters, consistent with the updated baseline.

clearvars;
clc;
close all;

scriptTimer = tic;

fprintf('\n============================================================\n');
fprintf(' Sensitivity analysis: Cd (downtime cost rate)\n');
fprintf('============================================================\n');

%% --------------------- Baseline parameters ------------------
parBase.alpha  = 1.02;
parBase.beta   = 1.2;
parBase.lambda = 0.25;

% Shock magnitude W ~ Gamma(alphaS,betaS), betaS is RATE.
parBase.alphaS = 2.5;
parBase.betaS  = 1;

% Imperfect PM ratio B ~ Beta(p,q).
parBase.p = 6;
parBase.q = 5;

parBase.L = 12.5;

parBase.Ci = 10;
parBase.Cm = 100;
parBase.Cp = 280;
parBase.Cc = 350;
parBase.Cd = 50;

targetField = 'Cd';
baselineValue = parBase.(targetField);

%% --------------------- Numerical options --------------------
opts.Nz       = 60;
opts.Nv       = 80;
opts.Ny       = 18;
opts.NtInt    = 80;
opts.Nconv    = 60;
opts.NuRel    = 18;
opts.NtRel    = 60;
opts.tRelMaxFactor = 2.5;
opts.useParallel = true;
opts.showProgress = true;
opts.progressEverySeconds = 3;
opts.minMassWarning = 5e-2;

%% --------------------- Sensitivity levels -------------------
sens.multipliers = [0.50, 1.50];

%% --------------------- Initial search ranges ----------------
sens.RILow  = 0.75;
sens.RIHigh = 0.99;
sens.RIStep = 0.01;

sens.MRatioLow  = 0.70;
sens.MRatioHigh = 0.90;
sens.MRatioStep = 0.01;

% Automatic expansion.
sens.expandRI     = 0.10;
sens.expandMRatio = 0.10;
sens.maxExpansionRounds = 3;

% Hard limits used only if the optimum reaches a boundary.
sens.RIHardLow  = 0.01;
sens.RIHardHigh = 0.99;
sens.MRatioHardLow  = 0.05;
sens.MRatioHardHigh = 0.99;

%% --------------------- Optimize -50% and +50% ---------------
nCase = numel(sens.multipliers);

Multiplier        = sens.multipliers(:);
ChangeFromBasePct = 100*(Multiplier-1);
ParameterValue    = baselineValue*Multiplier;

OptimalCostRate   = nan(nCase,1);
OptimalRI         = nan(nCase,1);
OptimalM          = nan(nCase,1);
OptimalMOverL     = nan(nCase,1);
ExpansionCount    = nan(nCase,1);
HitFinalBoundary  = false(nCase,1);

for icase = 1:nCase
    par = parBase;
    par.(targetField) = ParameterValue(icase);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('Case %d/%d: %s = %.10g (%+.0f%%%% from baseline)\n', ...
        icase, nCase, targetField, par.(targetField), ChangeFromBasePct(icase));
    fprintf('------------------------------------------------------------\n');

    label = sprintf('%s_%+.0fpercent', targetField, ChangeFromBasePct(icase));
    [best, info] = optimizeWithBoundaryExpansion(par, opts, sens, label);

    OptimalCostRate(icase) = best.CostRate;
    OptimalRI(icase)       = best.RI;
    OptimalM(icase)        = best.M;
    OptimalMOverL(icase)   = best.MRatio;
    ExpansionCount(icase)  = info.expansionCount;
    HitFinalBoundary(icase)= info.hitFinalBoundary;

    fprintf('  Optimal cost rate = %.8f\n', best.CostRate);
    fprintf('  Optimal RI        = %.4f\n', best.RI);
    fprintf('  Optimal M/L       = %.4f\n', best.MRatio);
end

Level = ["-50%"; "+50%"];

sensitivityTable = table( ...
    Level, Multiplier, ChangeFromBasePct, ParameterValue, ...
    OptimalCostRate, OptimalRI, OptimalM, OptimalMOverL, ...
    ExpansionCount, HitFinalBoundary);

fprintf('\n================ Sensitivity results: Cd ================\n');
disp(sensitivityTable);

%% --------------------- Save outputs --------------------------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

results = struct();
results.targetParameter = targetField;
results.baselineValue = baselineValue;
results.baselineParameters = parBase;
results.options = opts;
results.searchSettings = sens;
results.table = sensitivityTable;
results.totalRuntimeSeconds = toc(scriptTimer);

matFile = fullfile(thisDir,'sensitivity_Cd_50pct_results.mat');
csvFile = fullfile(thisDir,'sensitivity_Cd_50pct_results.csv');

save(matFile, 'results', 'sensitivityTable');
writetable(sensitivityTable, csvFile);

fprintf('\nSaved: %s\n', matFile);
fprintf('Saved: %s\n', csvFile);
fprintf('Total runtime: %.2f s\n', results.totalRuntimeSeconds);
fprintf('Done.\n');


%% ========================================================================
%  Adaptive grid optimization for sensitivity analysis
% ========================================================================
function [best, info] = optimizeWithBoundaryExpansion(par, opts, sens, label)

riLow  = sens.RILow;
riHigh = sens.RIHigh;
mLow   = sens.MRatioLow;
mHigh  = sens.MRatioHigh;

expansionCount = 0;
lastTab = table();

for iround = 0:sens.maxExpansionRounds

    RIGrid = makeSensitivityGrid(riLow, riHigh, sens.RIStep);
    MRatioGrid = makeSensitivityGrid(mLow, mHigh, sens.MRatioStep);
    MGrid = MRatioGrid * par.L;

    fprintf('\n  Search round %d: RI=[%.2f, %.2f], M/L=[%.2f, %.2f]\n', ...
        iround+1, RIGrid(1), RIGrid(end), MRatioGrid(1), MRatioGrid(end));

    lastTab = gridSearchPolicy(RIGrid, MGrid, par, opts, label);

    valid = isfinite(lastTab.CostRate);
    if ~any(valid)
        error('No valid policy was found in the current sensitivity grid.');
    end

    validIdx = find(valid);
    [bestCost, j] = min(lastTab.CostRate(valid));
    idxBest = validIdx(j);

    bestRI = lastTab.RI(idxBest);
    bestM  = lastTab.M(idxBest);
    bestMRatio = bestM / par.L;

    tolRI = 0.25 * sens.RIStep;
    tolM  = 0.25 * sens.MRatioStep;

    hitRILow  = abs(bestRI - RIGrid(1)) <= tolRI;
    hitRIHigh = abs(bestRI - RIGrid(end)) <= tolRI;
    hitMLow   = abs(bestMRatio - MRatioGrid(1)) <= tolM;
    hitMHigh  = abs(bestMRatio - MRatioGrid(end)) <= tolM;

    newRiLow  = riLow;
    newRiHigh = riHigh;
    newMLow   = mLow;
    newMHigh  = mHigh;

    if hitRILow && riLow > sens.RIHardLow + 1e-12
        newRiLow = max(sens.RIHardLow, riLow - sens.expandRI);
    end
    if hitRIHigh && riHigh < sens.RIHardHigh - 1e-12
        newRiHigh = min(sens.RIHardHigh, riHigh + sens.expandRI);
    end
    if hitMLow && mLow > sens.MRatioHardLow + 1e-12
        newMLow = max(sens.MRatioHardLow, mLow - sens.expandMRatio);
    end
    if hitMHigh && mHigh < sens.MRatioHardHigh - 1e-12
        newMHigh = min(sens.MRatioHardHigh, mHigh + sens.expandMRatio);
    end

    expanded = (abs(newRiLow-riLow) > 1e-12) || ...
               (abs(newRiHigh-riHigh) > 1e-12) || ...
               (abs(newMLow-mLow) > 1e-12) || ...
               (abs(newMHigh-mHigh) > 1e-12);

    if ~expanded
        break;
    end

    riLow  = newRiLow;
    riHigh = newRiHigh;
    mLow   = newMLow;
    mHigh  = newMHigh;
    expansionCount = expansionCount + 1;
end

valid = isfinite(lastTab.CostRate);
validIdx = find(valid);
[bestCost, j] = min(lastTab.CostRate(valid));
idxBest = validIdx(j);

best = struct();
best.RI       = lastTab.RI(idxBest);
best.M        = lastTab.M(idxBest);
best.MRatio   = best.M / par.L;
best.CostRate = bestCost;

finalRIGrid = unique(lastTab.RI);
finalMGrid  = unique(lastTab.M);
finalMRatioGrid = finalMGrid / par.L;

tolRI = 0.25 * sens.RIStep;
tolM  = 0.25 * sens.MRatioStep;

info = struct();
info.expansionCount = expansionCount;
info.finalRILow      = min(finalRIGrid);
info.finalRIHigh     = max(finalRIGrid);
info.finalMRatioLow  = min(finalMRatioGrid);
info.finalMRatioHigh = max(finalMRatioGrid);
info.hitFinalBoundary = ...
    abs(best.RI-info.finalRILow) <= tolRI || ...
    abs(best.RI-info.finalRIHigh) <= tolRI || ...
    abs(best.MRatio-info.finalMRatioLow) <= tolM || ...
    abs(best.MRatio-info.finalMRatioHigh) <= tolM;

if info.hitFinalBoundary
    warning(['Final optimum is still on a search boundary. ' ...
             'Inspect this case and expand the range further if needed.']);
end
end

function x = makeSensitivityGrid(lo, hi, step)
n = floor((hi-lo)/step + 1e-10);
x = lo + (0:n)*step;
if x(end) < hi - 1e-10
    x = [x, hi];
elseif abs(x(end)-hi) <= 1e-10
    x(end) = hi;
end
x = unique(round(x,12));
end

%% ========================================================================
%  Main policy evaluator
% ========================================================================
function rep = evaluatePolicy(RI, M, par, opts, doReliabilityCurve)
if nargin < 5
    doReliabilityCurve = false;
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
    row = buildOneStepRow(states(i), RI, M, par, opts, z, wz, v, wv, pmSurvVec);
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

% MTTF from Eq. (61): m = source + P_survival*m.
I = eye(Ns);
mVec = (I - Psurv) \ cond.mttfSource;
MTTF_stationary = sum(piVec .* mVec);
MTTF_fromNew    = mVec(1);

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
        RI, M, par, opts, cond.tau, PrelSurv);

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
    'minMTTFLinearMatrixReciprocalCondition', rcond(I-Psurv), ...
    'MTTFfromReliabilityCurve', MTTFcurve, ...
    'relativeMTTFcurveError', relativeMTTFcurveError, ...
    'reliabilityTail', reliabilityTail);

maxRawComponentError = max([ ...
    rep.diagnostics.maxRawTransitionMassError, ...
    rep.diagnostics.maxDoNothingMassError, ...
    rep.diagnostics.maxPMSuccessMassError, ...
    rep.diagnostics.maxPMEventMassError]);

if maxRawComponentError > opts.minMassWarning
    warning(['Max raw transition/component mass error is %.3g. ', ...
        'Consider increasing opts.Nz/Nv/Ny.'], maxRawComponentError);
end
end
function row = buildOneStepRow(s, RI, M, par, opts, z, wz, v, wv, pmSurvVec)
L = par.L;
lam = par.lambda;
Nz = numel(z);

% Reliability-triggered scheduled inspection time.
tau = tauFromReliabilityThreshold(s, RI, par);
aTau = par.alpha * tau;
eNoShock = exp(-lam * tau);

pmPdfMat = pmPdfMatrix(z, v, par.p, par.q);     % Nz x Nv, f_{S*}(z|v)

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
            s, y, z, wz, v, wv, pmSurvVec, M, par, opts.Nconv);

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
mttfSource = mttfSourceIntegral(s, tau, par, opts.NtInt);

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
function Rstate = reliabilityCurveByRecursion(tGrid, states, z, wz, v, wv, pmSurvVec, RI, M, par, opts, tauVec, PrelSurv)
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
                    s, u, z, wz, v, wv, pmSurvVec, M, par, opts.Nconv);

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

function Prow = shockSurvivalRowAtU(s, u, z, wz, v, wv, pmSurvVec, M, par, Nconv)
% Density-in-time survival transition row at a first shock occurring at u.
% Corrective-replacement paths are excluded. The conditional state
% transition components are corrected in exactly the same way as in the
% one-step kernel used for stationary cost and MTTF calculations.

comp = shockTransitionComponentsAtU( ...
    s, u, z, wz, v, wv, pmSurvVec, M, par, Nconv);

coeff = par.lambda * exp(-par.lambda * u);

Prow = zeros(1,numel(z)+1);
Prow(1) = coeff * comp.Ppr;
Prow(2:end) = coeff * ...
    ((comp.pmSuccessDensity + comp.doNothingDensity) .* wz)';
end

function comp = shockTransitionComponentsAtU(s, u, z, wz, v, wv, pmSurvVec, M, par, Nconv)
% Conditional transition components given that the first shock occurs at u.
% This helper is shared by the one-step kernel and the reliability recursion
% so both calculations use exactly the same PM/do-nothing mass correction.

pmPdfMat = pmPdfMatrix(z, v, par.p, par.q);

% Partition of the pre-maintenance state at the shock inspection:
%   below M, in [M,L), or failed at/above L.
survL = sumGammaSurvRate(par.L - s, u, par, Nconv);
survM = sumGammaSurvRate(M - s, u, par, Nconv);

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

gu_v = sumGammaPdfRate(v - s, u, par);
gu_v(v <= s) = 0;

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

gu_z = sumGammaPdfRate(z - s, u, par);
gu_z(z <= s) = 0;
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
            rep = evaluatePolicy(RIList(k), MList(k), par, opts, false);
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
            rep = evaluatePolicy(RIList(k), MList(k), par, opts, false);
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
for iter = 1:80
    mid = 0.5*(lo+hi);
    if fun(mid) > 0
        lo = mid;
    else
        hi = mid;
    end
end
tau = 0.5*(lo+hi);
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
% PDF of U = X_y + W, where X_y~Gamma(alpha*y,beta), W~Gamma(alphaS,betaS).
% All beta parameters are rates.
% The fast path uses the confluent hypergeometric expression when MATLAB's
% hypergeom is available. Otherwise, a fixed midpoint convolution fallback is used.
pdf = zeros(size(u));
idx = u > 0;
if ~any(idx(:))
    return;
end
uPos = u(idx);
a = par.alpha * y;
c = par.alphaS;
b = par.beta;
d = par.betaS;
if a <= 1e-12
    pdf(idx) = gammaPdfRate(uPos, c, d);
    return;
end

persistent useHypergeom
if isempty(useHypergeom)
    useHypergeom = (exist('hypergeom','file') == 2 || exist('hypergeom','builtin') == 5);
end

if useHypergeom
    try
        z = -(d - b) .* uPos;
        % hypergeom can return tiny imaginary parts due to numerical roundoff.
        hg = hypergeom(c, a + c, z);
        hg = real(hg);
        if any(~isfinite(hg)) || any(hg <= 0)
            error('BadHypergeomValue');
        end
        logp = a*log(b) + c*log(d) - gammaln(a+c) + (a+c-1).*log(uPos) - b.*uPos + log(hg);
        pdf(idx) = exp(logp);
        pdf(~isfinite(pdf)) = 0;
        pdf = max(pdf,0);
        return;
    catch
        % Fall through to numerical convolution.
        useHypergeom = false;
    end
end

% Fallback: fixed midpoint convolution. Slower but avoids Symbolic Math Toolbox dependency.
Npdf = 48;
out = zeros(size(uPos));
for ii = 1:numel(uPos)
    uu = uPos(ii);
    w = ((1:Npdf)' - 0.5) / Npdf * uu;
    ww = uu / Npdf;
    out(ii) = sum(ww .* gammaPdfRate(uu - w, a, b) .* gammaPdfRate(w, c, d));
end
pdf(idx) = out;
pdf(~isfinite(pdf)) = 0;
pdf = max(pdf,0);
end

function surv = sumGammaSurvRate(u, y, par, Nconv)
% Survival P(X_y + W >= u). Uses fixed midpoint quadrature:
% P(W>=u) + integral_0^u P(X_y >= u-w) f_W(w) dw.
if u <= 0
    surv = 1;
    return;
end
a = par.alpha * y;
if a <= 1e-12
    surv = gammaSurvRate(u, par.alphaS, par.betaS);
    return;
end
w = ((1:Nconv)' - 0.5) / Nconv * u;
ww = u / Nconv;
fw = gammaPdfRate(w, par.alphaS, par.betaS);
sx = gammaSurvRate(u - w, a, par.beta);
surv = gammaSurvRate(u, par.alphaS, par.betaS) + sum(ww .* sx .* fw);
surv = min(max(surv,0),1);
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

savefig(fig, 'coarse_cost_contour_v5.fig');
try
    exportgraphics(fig, 'coarse_cost_contour_v5.png', 'Resolution', 300);
catch
    print(fig, 'coarse_cost_contour_v5.png', '-dpng', '-r300');
end
fprintf('Saved coarse contour plot to coarse_cost_contour_v5.fig and coarse_cost_contour_v5.png\n');
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

savefig(fig, 'final_reliability_curve_v5.fig');
try
    exportgraphics(fig, 'final_reliability_curve_v5.png', 'Resolution', 300);
catch
    print(fig, 'final_reliability_curve_v5.png', '-dpng', '-r300');
end
fprintf('Saved final reliability plot to final_reliability_curve_v5.fig and final_reliability_curve_v5.png\n');
end

function validatePolicy(RI, M, par)
if ~(RI > 0 && RI < 1)
    error('RI must satisfy 0 < RI < 1.');
end
if ~(M > 0 && M < par.L)
    error('M must satisfy 0 < M < L.');
end
end
