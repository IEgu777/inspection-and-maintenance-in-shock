function results = run_section65_validation(quickTest)
% RUN_SECTION65_VALIDATION
% Independent, path-based Monte Carlo verification for Section 6.5.
% Companion: section65_analytical_reference.m (the original V6 formula solver).
%
% HOW TO RUN:
%   Add both .m files to the MATLAB path, then:
%       results = run_section65_validation(true);   % smoke test
%       results = run_section65_validation(false);  % full experiment
%
% The default counts below are DESIGN TARGETS for the planned experiment;
% they have NOT been executed or validated by the author of this script.
% Run the quick test first by setting cfg.quickTest = true.
%
% Key independence: the Monte Carlo code samples Gamma-process increments,
% Poisson shocks and Beta PM effects. It never samples maintenance outcomes
% from probabilities computed by the analytical model.
%
% Important distinctions:
%   * P_pm includes cases that subsequently lead to preventive replacement.
%   * Downtime begins at the actual continuous-degradation threshold crossing,
%     even when the failure is only detected at the next inspection.
%   * Reliability/MTTF STOP at the FIRST failure (replacement does not
%     restore first-failure reliability).
%   * The stationary initial distribution for MC reliability is sampled
%     from independently simulated long-run post-maintenance states.
%
% Outputs in ./section65_results:
%   Table_section65_single_cycle.csv and .tex
%   Table_section65_transition_cdf.csv
%   Table_section65_system_performance.csv and .tex
%   Figure_section65_reliability.png
%   Data_section65_reliability.csv
%   Table_section65_grid_stability.csv (if enabled)
%   section65_all_results.mat

clc;
if nargin < 1, quickTest = false; end

%% ------------------------ Experiment configuration -----------------------
cfg.randomSeed = 20261118;
cfg.quickTest = logical(quickTest); % true: smoke test only, NOT paper results
cfg.runGridStability = true;      % numerical Eq. (48) convergence study
cfg.showProgress = true;

% Policy fixed at Section 6.2 values; no grid-search optimization here.
policy.RI = 0.920;
policy.M  = 9.625;

% Parameterization matches run_shock_cbm_case_v6_fast.m exactly.
% Gamma(shape, RATE): use randg(shape)/rate for random generation.
par.alpha = 1.02;
par.beta = 1.20;
par.lambda = 0.25;
par.alphaS = 2.5;
par.betaS = 1.0;
par.p = 6;
par.q = 5;
par.L = 12.5;
par.Ci = 10;
par.Cm = 100;
par.Cp = 280;
par.Cc = 350;
par.Cd = 50;

% Use four representative starting states. Easily replace these later.
% All states must lie in [0,M), because these are POST-maintenance states.
cfg.sValues = [0, 0.25, 0.50, 0.85] * policy.M;

% Target simulation sizes for the manuscript, to be confirmed by running.
cfg.nSinglePerState = 100000;
cfg.nLongRep = 100;
cfg.longHorizon = 10000;       % time units AFTER burn-in
cfg.burnInHorizon = 1000;     % discard initial transient cycles
cfg.nFirstFailureNew = 20000;
cfg.nFirstFailureStationary = 20000;

% Gamma bridge brackets a continuous-degradation first-passage time.
% Interval uncertainty <= inspectionInterval/2^bridgeIterations.
cfg.bridgeIterations = 18;

% State-dependent scheduled inspections: precompute tau(s) defined by
% P[s + X_tau < L] = RI. Interpolation is only an acceleration technique;
% it does NOT use analytical one-step transition probabilities.
cfg.nScheduleNodes = 2401;
cfg.nStationaryReservoir = 80000;
cfg.maxCyclesToFirstFailure = 50000;

if cfg.quickTest
    cfg.nSinglePerState = 1000;
    cfg.nLongRep = 4;
    cfg.longHorizon = 150;
    cfg.burnInHorizon = 50;
    cfg.nFirstFailureNew = 400;
    cfg.nFirstFailureStationary = 400;
    cfg.nScheduleNodes = 241;
    cfg.nStationaryReservoir = 2000;
    cfg.runGridStability = false;
end

assert(policy.RI>0 && policy.RI<1 && policy.M>0 && policy.M<par.L);
assert(all(cfg.sValues>=0 & cfg.sValues<policy.M));

% Analytical discretization: same baseline settings as the supplied V6.
opts.Nz = 60;
opts.Nv = 80;
opts.Ny = 18;
opts.NtInt = 80;
opts.Nconv = 60;
opts.NuRel = 18;
opts.NtRel = 60;
opts.tRelMaxFactor = 2.5;
opts.computeReliabilityCurve = true;
opts.useParallel = false;
opts.showProgress = false;
opts.progressEverySeconds = 3;
opts.minMassWarning = 5e-2;
opts.warnOnMass = true;

outputDir = fullfile(fileparts(mfilename('fullpath')), 'section65_results');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
rng(cfg.randomSeed, 'twister');
clockAll = tic;

fprintf('\n=== Section 6.5: independent simulation verification ===\n');
fprintf('RI=%.3f, M=%.3f, L=%.3f; alpha=%.3f, beta(rate)=%.3f\n', ...
    policy.RI, policy.M, par.L, par.alpha, par.beta);
fprintf('Shock: lambda=%.3f, alphaS=%.3f, betaS(rate)=%.3f\n', ...
    par.lambda, par.alphaS, par.betaS);
if cfg.quickTest
    fprintf('*** QUICK TEST MODE: do NOT publish these results. ***\n');
end

%% ----------------- A. Existing analytical formulas (V6) -----------------
fprintf('\n[1/5] Evaluating V6 analytical formulas at the FIXED policy ...\n');
t0 = tic;
[ana, conditional] = section65_analytical_reference( ...
    policy.RI, policy.M, par, opts, cfg.sValues, true);
analyticalSeconds = toc(t0);
fprintf('  V6 analytical cost %.8f, availability %.8f, MTTF %.8f; time %.1fs\n', ...
    ana.costRate, ana.averageAvailability, ana.MTTF_stationary, analyticalSeconds);

% The published Section 6.2 values should be checked, NOT silently forced.
if abs(ana.costRate - 37.623) > 0.20
    warning(['The V6 cost (%.6f) differs from the manuscript value 37.623 ', ...
        'by more than 0.20. Check source parameters, code version and grids.'], ...
        ana.costRate);
end

%% -------- B. Independently calculate the inspection schedule tau(s) ----
fprintf('\n[2/5] Preparing the inspection-time schedule tau(s) ...\n');
schedule = makeInspectionSchedule(par, policy, cfg.nScheduleNodes);
% Quantify schedule interpolation error independently of the MC simulation.
checkS = linspace(0.003*policy.M, 0.997*policy.M, 49);
tauExact = arrayfun(@(s) solveTauExactly(s,par,policy.RI), checkS);
tauApprox = scheduleTau(checkS,schedule);
maxScheduleError = max(abs(tauExact-tauApprox));
fprintf('  Max absolute tau interpolation error (49 test states) = %.3e\n', ...
    maxScheduleError);
if maxScheduleError > 1e-3
    warning('Scheduled inspection lookup is not accurate enough; increase nScheduleNodes.');
end

%% ------------------- C. Conditional single-cycle checks -----------------
fprintf('\n[3/5] Independent MC for single-cycle quantities ...\n');
nS = numel(cfg.sValues);
nMeasures = 5;
measureLabels = {'mu(s)'; 'P_pm(s)'; 'P_pr(s)'; 'P_cr(s)'; 'd(s)'};
startingState = zeros(nS*nMeasures,1);
measureName = cell(nS*nMeasures,1);
analyticalValue = nan(nS*nMeasures,1);
mcValue = nan(nS*nMeasures,1);
absoluteError = nan(nS*nMeasures,1);
relativeErrorPct = nan(nS*nMeasures,1);
singleMC = struct([]);
cdfRows = zeros(nS*3,5); % s, cut-off, analytical CDF, MC CDF, abs error

for k = 1:nS
    s = cfg.sValues(k);
    % For a fixed s use EXACT root, eliminating schedule interpolation bias.
    tauFixed = solveTauExactly(s,par,policy.RI);
    sums = zeros(1,5);
    nextStateSamples = zeros(cfg.nSinglePerState,1);
    for j = 1:cfg.nSinglePerState
        c = sampleInspectionCycle(s, tauFixed, par, policy, cfg.bridgeIterations);
        sums = sums + [c.dt, c.pm, c.pr, c.cr, c.downtime];
        nextStateSamples(j) = c.nextState;
    end
    mc = sums / cfg.nSinglePerState;
    ref = [conditional(k).mu, conditional(k).Ppm, conditional(k).Ppr, ...
        conditional(k).Pcr, conditional(k).downtime];
    singleMC(k).s = s; %#ok<AGROW>
    singleMC(k).tau = tauFixed;
    singleMC(k).estimates = mc;
    singleMC(k).reference = ref;
    % Additional independent verification of the transition kernel P(s,dz):
    % compare empirical post-maintenance CDF with the V6 numerical CDF.
    for h=1:3
        x=policy.M * [0.2,0.5,0.8];
        cutoff=x(h);
        cdfAnalytic=conditional(k).fullMass0 + ...
            sum(conditional(k).fullInteriorDensity(:).* ...
                conditional(k).wz(:).*(conditional(k).z(:)<=cutoff));
        cdfMC=mean(nextStateSamples<=cutoff);
        cdfRows((k-1)*3+h,:)=[s,cutoff,cdfAnalytic,cdfMC, ...
            abs(cdfAnalytic-cdfMC)];
    end

    for m = 1:nMeasures
        row = (k-1)*nMeasures + m;
        startingState(row) = s;
        measureName{row} = measureLabels{m};
        analyticalValue(row) = ref(m);
        mcValue(row) = mc(m);
        absoluteError(row) = abs(ref(m)-mc(m));
        if abs(ref(m)) >= 1e-8
            relativeErrorPct(row) = 100 * absoluteError(row)/abs(ref(m));
        end
    end
    fprintf('  s=%7.4f, tau=%.5f: mu=%.5f(MC %.5f), d=%.5f(MC %.5f)\n', ...
        s,tauFixed,ref(1),mc(1),ref(5),mc(5));
end

tableSingle = table(startingState, measureName, analyticalValue, mcValue, ...
    absoluteError, relativeErrorPct, 'VariableNames', ...
    {'InitialState','Measure','Analytical','MonteCarlo','AbsoluteError','RelativeErrorPct'});
writetable(tableSingle,fullfile(outputDir,'Table_section65_single_cycle.csv'));
writeLatexSingleCycle(tableSingle, cfg.sValues, outputDir);
tableTransitionCDF=array2table(cdfRows,'VariableNames', ...
    {'InitialState','StateCutoff','AnalyticalCDF','MonteCarloCDF','AbsoluteError'});
writetable(tableTransitionCDF,fullfile(outputDir,'Table_section65_transition_cdf.csv'));

%% ------------ D. Long-run cost, availability, stationary states ---------
fprintf('\n[4/5] Independent long-run MC and FIRST-failure MC ...\n');
t0 = tic;
[longMC, empiricalStates] = simulateLongRun(par,policy,schedule,cfg);

% Reliability from a brand-new condition s=0.
firstFailureNew = simulateFirstFailureSample( ...
    zeros(cfg.nFirstFailureNew,1), par,policy,schedule,cfg);

% MC stationary initial states: sample from simulated inspection-epoch
% states. In particular, do NOT use ana.stationary.piVec to generate them.
empiricalIdx = randi(numel(empiricalStates), cfg.nFirstFailureStationary, 1);
initialStationaryStates = empiricalStates(empiricalIdx);
firstFailureStationary = simulateFirstFailureSample( ...
    initialStationaryStates, par,policy,schedule,cfg);
failureSeconds = toc(t0);

mcMTTFnew = mean(firstFailureNew);
mcMTTFstationary = mean(firstFailureStationary);

systemMeasure = {'C_infty'; 'AverageAvailability'; 'MTTF_stationary'; 'MTTF_new'; ...
    'MeanCycleLength'; 'P_pm_per_cycle'; 'P_pr_per_cycle'; ...
    'P_cr_per_cycle'; 'Downtime_per_cycle'};
refLong = [ana.costRate; ana.averageAvailability; ana.MTTF_stationary; ...
    ana.MTTF_fromNew; ana.expectedCycleLength; ...
    ana.expectedCountsPerCycle.PM; ana.expectedCountsPerCycle.preventiveReplacement; ...
    ana.expectedCountsPerCycle.correctiveReplacement; ana.expectedDowntimePerCycle];
mcLong = [longMC.costRate; longMC.availability; mcMTTFstationary; ...
    mcMTTFnew; longMC.meanCycleLength; longMC.pmPerCycle; longMC.prPerCycle; ...
    longMC.crPerCycle; longMC.downPerCycle];
absLong = abs(refLong-mcLong);
relLong = 100*absLong./max(abs(refLong),1e-8);
tableSystem = table(systemMeasure,refLong,mcLong,absLong,relLong, ...
    'VariableNames',{'Measure','Analytical','MonteCarlo','AbsoluteError','RelativeErrorPct'});
writetable(tableSystem,fullfile(outputDir,'Table_section65_system_performance.csv'));
writeLatexSystem(tableSystem,outputDir);

fprintf('  Cost: analytical %.6f versus MC %.6f\n',ana.costRate,longMC.costRate);
fprintf('  Availability: analytical %.6f versus MC %.6f\n', ...
    ana.averageAvailability,longMC.availability);
fprintf('  State-0 stationary atom: analytical %.6f versus MC %.6f\n', ...
    ana.stationary.a0,longMC.stationaryMassAtZero);
fprintf('  MTTF (new): analytical %.6f versus MC %.6f\n', ...
    ana.MTTF_fromNew,mcMTTFnew);
fprintf('  MTTF (stationary): analytical %.6f versus MC %.6f\n', ...
    ana.MTTF_stationary,mcMTTFstationary);

if isempty(ana.reliability.t)
    error('Analytical reliability curve was not computed.');
end
relTime = ana.reliability.t(:);
mcRnew = arrayfun(@(t) mean(firstFailureNew>t), relTime);
mcRstationary = arrayfun(@(t) mean(firstFailureStationary>t), relTime);
relData = table(relTime, ana.reliability.RfromNew(:), mcRnew, ...
    ana.reliability.Rstationary(:), mcRstationary, ...
    'VariableNames',{'Time','AnalyticalNew','MonteCarloNew', ...
    'AnalyticalStationary','MonteCarloStationary'});
writetable(relData,fullfile(outputDir,'Data_section65_reliability.csv'));

fig = figure('Name','Section 6.5 reliability validation','Color','w', ...
    'Position',[100,100,1040,420]);
subplot(1,2,1);
plot(relTime,ana.reliability.RfromNew,'-','LineWidth',1.7); hold on;
plot(relTime,mcRnew,'o','MarkerSize',4,'LineWidth',1);
xlabel('Time'); ylabel('Reliability');
title('New system (s = 0)'); legend('Analytical','Monte Carlo','Location','best');
grid on; ylim([0,1.03]);
subplot(1,2,2);
plot(relTime,ana.reliability.Rstationary,'-','LineWidth',1.7); hold on;
plot(relTime,mcRstationary,'o','MarkerSize',4,'LineWidth',1);
xlabel('Time'); ylabel('Reliability');
title('Stationary initial state'); legend('Analytical','Monte Carlo','Location','best');
grid on; ylim([0,1.03]);
try
    exportgraphics(fig,fullfile(outputDir,'Figure_section65_reliability.png'), ...
        'Resolution',300);
catch
    print(fig,fullfile(outputDir,'Figure_section65_reliability.png'),'-dpng','-r300');
end

%% ------------ E. Grid refinement for Eq. (48) (NOT Monte Carlo) ----------
fprintf('\n[5/5] Eq. (48) grid-refinement test ...\n');
gridTable = table();
if cfg.runGridStability
    gridTable = runReliabilityGridStudy(par,policy,opts,ana,analyticalSeconds);
    writetable(gridTable,fullfile(outputDir,'Table_section65_grid_stability.csv'));
else
    fprintf('  Skipped (cfg.runGridStability=false).\n');
end

results = struct();
results.cfg = cfg;
results.policy = policy;
results.parameters = par;
results.analyticalOptions = opts;
results.analytical = ana;
results.conditionalAnalytical = conditional;
results.singleCycleMonteCarlo = singleMC;
results.longRunMonteCarlo = longMC;
results.stationaryStatePool = empiricalStates;
results.firstFailureNew = firstFailureNew;
results.firstFailureStationary = firstFailureStationary;
results.singleCycleTable = tableSingle;
results.transitionCDFTable = tableTransitionCDF;
results.systemPerformanceTable = tableSystem;
results.reliabilityTable = relData;
results.gridStabilityTable = gridTable;
results.maxScheduleInterpolationError = maxScheduleError;
results.analyticalSeconds = analyticalSeconds;
results.simulationSeconds = failureSeconds;
results.totalSeconds = toc(clockAll);
save(fullfile(outputDir,'section65_all_results.mat'),'results','-v7.3');

fprintf('\n--- Results written to %s ---\n', outputDir);
fprintf('Total elapsed time: %.1f seconds\n',results.totalSeconds);
fprintf('For publication, verify numerical agreement BEFORE writing conclusions.\n');
end

%% ========================================================================
%  Independent path simulator. No analytic event-probability sampling here.
% ========================================================================
function c = sampleInspectionCycle(s,tau,par,policy,bridgeIterations)
% Sample one cycle from the post-maintenance state s through the NEXT
% inspection and the instantaneous maintenance action at that inspection.
% Samples are exact at inspection time. Continuous-failure time is
% localized by a conditional Gamma bridge, not a fixed Euler time step.
if par.lambda > 0
    nextShock = -log(max(realmin,rand())) / par.lambda;
else
    nextShock = inf;
end
shockTriggered = (nextShock < tau);
dt = min(nextShock,tau);

% Independent Gamma-process increment over dt, shape alpha*dt, rate beta.
if dt > 0
    wear = randg(par.alpha * dt) / par.beta;
else
    wear = 0;
end
preShockLevel = s + wear;

% First crossing due to continuous degradation may PRECEDE the inspection.
firstFailure = inf;
downtime = 0;
if preShockLevel >= par.L
    crossAt = gammaBridgeFirstPassage( ...
        par.L-s,wear,dt,par.alpha,bridgeIterations);
    firstFailure = crossAt;
    downtime = max(0,dt-crossAt);
end

% At a shock-triggered inspection the shock occurs at the inspection time.
if shockTriggered
    shockDamage = randg(par.alphaS)/par.betaS;
else
    shockDamage = 0;
end
observedState = preShockLevel + shockDamage;

% An instantaneous shock can cause first failure AT inspection, with no
% preceding hidden downtime (inspection/repair are instantaneous).
if isinf(firstFailure) && observedState >= par.L
    firstFailure = dt;
end

% PM happens before the decision to replace preventively. Thus one cycle
% may have BOTH pm=1 and pr=1; these events are NOT mutually exclusive.
pm = 0; pr = 0; cr = 0;
if observedState >= par.L
    cr = 1;
    nextState = 0;
elseif observedState >= policy.M
    pm = 1;
    a = randg(par.p);
    b = randg(par.q);
    pmRatio = a/(a+b);    % Beta(p,q), direct random generation
    postPM = observedState * pmRatio;
    if postPM >= policy.M
        pr = 1;
        nextState = 0;
    else
        nextState = postPM;
    end
else
    nextState = observedState;
end

c.dt = dt;
c.pm = pm;
c.pr = pr;
c.cr = cr;
c.downtime = downtime;
c.nextState = nextState;
c.observedState = observedState;
c.shockTriggered = shockTriggered;
c.firstFailureOffset = firstFailure;
c.cost = par.Ci + par.Cm*pm + par.Cp*pr + ...
    par.Cc*cr + par.Cd*downtime;
end

function crossing = gammaBridgeFirstPassage(target,totalWear,totalTime,alpha,nSteps)
% Given totalWear = X(totalTime), sample the conditional Gamma bridge.
% The Gamma-process bridge property is
% (X(mid)-X(left))/(X(right)-X(left)) ~ Beta(alpha*(mid-left),
%                                                  alpha*(right-mid)).
% Binary subdivisions bracket the first hitting time. This is NOT a
% maintenance-event kernel or any other analytical formula being tested.
assert(target>0 && totalWear>=target && totalTime>0);
leftT=0; rightT=totalTime;
leftWear=0; rightWear=totalWear;
for k=1:nSteps
    midT=(leftT+rightT)/2;
    aa=alpha*(midT-leftT);
    bb=alpha*(rightT-midT);
    ga=randg(aa);
    gb=randg(bb);
    share=ga/(ga+gb);
    if ~isfinite(share)
        % Rare numerical underflow for very small shapes; use their limit.
        share = double(rand() < aa/(aa+bb));
    end
    midWear=leftWear+(rightWear-leftWear)*share;
    if midWear>=target
        rightT=midT;
        rightWear=midWear;
    else
        leftT=midT;
        leftWear=midWear;
    end
end
crossing=(leftT+rightT)/2;
end

%% ========================================================================
%  Cached scheduled inspection time; root equation is the policy definition.
% ========================================================================
function schedule=makeInspectionSchedule(par,policy,nNodes)
schedule.states = linspace(0,policy.M,nNodes)';
schedule.taus = zeros(nNodes,1);
for i=1:nNodes
    schedule.taus(i) = solveTauExactly(schedule.states(i),par,policy.RI);
end
schedule.functionHandle = griddedInterpolant( ...
    schedule.states,schedule.taus,'pchip','nearest');
end

function tau=scheduleTau(s,schedule)
tau=schedule.functionHandle(s);
end

function tau=solveTauExactly(s,par,RI)
margin=par.L-s;
if margin<=0, tau=0; return; end
fun=@(t) gammainc(par.beta*margin,par.alpha*t,'lower')-RI;
lo=1e-10;
hi=max(1,margin*par.beta/max(par.alpha,eps));
while fun(hi)>0
    hi=2*hi;
    if hi>1e10, error('Unable to bracket scheduled time tau.'); end
end
if abs(fun(hi))<1e-14
    tau=hi;
else
    tau=fzero(fun,[lo,hi]);
end
end

%% ========================================================================
%  Long-run MC: cycles continue across PM/PR/CR; no independent resets.
% ========================================================================
function [out,stationaryStates]=simulateLongRun(par,policy,schedule,cfg)
accTime=0; accDowntime=0; accCost=0; accCycles=0;
accPM=0; accPR=0; accCR=0;
stateSeen=0;
stationaryStates=zeros(cfg.nStationaryReservoir,1);
stateFilled=0;
repCosts=zeros(cfg.nLongRep,1);
repAvail=zeros(cfg.nLongRep,1);
for r=1:cfg.nLongRep
    s=0;
    t=0;
    repTime=0; repDown=0; repCost=0;
    % Discard cycles which start before burnInHorizon. The next complete
    % cycles form the long-run measurement window (no partial cycle bias).
    while t < cfg.burnInHorizon+cfg.longHorizon
        tau=scheduleTau(s,schedule);
        c=sampleInspectionCycle(s,tau,par,policy,cfg.bridgeIterations);
        if t>=cfg.burnInHorizon
            accTime=accTime+c.dt;
            accDowntime=accDowntime+c.downtime;
            accCost=accCost+c.cost;
            accCycles=accCycles+1;
            accPM=accPM+c.pm;
            accPR=accPR+c.pr;
            accCR=accCR+c.cr;
            repTime=repTime+c.dt;
            repDown=repDown+c.downtime;
            repCost=repCost+c.cost;

            % Reservoir sampling of the EMBEDDED inspection-epoch chain.
            % This approximates the stationary pi(ds) independently of V6.
            stateSeen=stateSeen+1;
            if stateFilled<cfg.nStationaryReservoir
                stateFilled=stateFilled+1;
                stationaryStates(stateFilled)=s;
            else
                if rand()<cfg.nStationaryReservoir/stateSeen
                    idx=randi(cfg.nStationaryReservoir);
                    stationaryStates(idx)=s;
                end
            end
        end
        t=t+c.dt;
        s=c.nextState;
    end
    repCosts(r)=repCost/max(repTime,eps);
    repAvail(r)=1-repDown/max(repTime,eps);
    if cfg.showProgress && (mod(r,max(1,round(cfg.nLongRep/10)))==0 || r==1)
        fprintf('  Long-run %d/%d: cost %.5f, avail %.5f\n', ...
            r,cfg.nLongRep,repCosts(r),repAvail(r));
    end
end
stationaryStates=stationaryStates(1:stateFilled);
assert(~isempty(stationaryStates),'No post-burn-in states sampled.');
out.costRate=accCost/accTime;
out.availability=1-accDowntime/accTime;
out.meanCycleLength=accTime/accCycles;
out.pmPerCycle=accPM/accCycles;
out.prPerCycle=accPR/accCycles;
out.crPerCycle=accCR/accCycles;
out.downPerCycle=accDowntime/accCycles;
out.totalTime=accTime;
out.totalDowntime=accDowntime;
out.totalCost=accCost;
out.totalCycles=accCycles;
out.replicationCosts=repCosts;
out.replicationAvailabilities=repAvail;
out.reservoirStateCount=stateFilled;
out.stationaryMassAtZero=mean(stationaryStates==0);
end

%% ========================================================================
%  First-failure MC: stop at actual threshold crossing, even during a cycle.
% ========================================================================
function tf=simulateFirstFailureSample(initialStates,par,policy,schedule,cfg)
n=numel(initialStates);
tf=zeros(n,1);
for r=1:n
    s=initialStates(r);
    t=0;
    for k=1:cfg.maxCyclesToFirstFailure
        tau=scheduleTau(s,schedule);
        c=sampleInspectionCycle(s,tau,par,policy,cfg.bridgeIterations);
        if isfinite(c.firstFailureOffset)
            tf(r)=t+c.firstFailureOffset;
            break;
        end
        t=t+c.dt;
        s=c.nextState;
    end
    if tf(r)==0
        error('No first failure within %d cycles; increase the limit.', ...
            cfg.maxCyclesToFirstFailure);
    end
end
end

%% ========================================================================
%  Numerical stability: vary BOTH state and time quadrature resolutions.
% ========================================================================
function outTable=runReliabilityGridStudy(par,policy,baselineOpts,baselineRep,baselineRuntime)
% Coarse, baseline, fine - cost/MTTF and Eq.(48) reliability on all levels.
% Fine setting is deliberately expensive; disable cfg.runGridStability
% during MC development or adjust resolutions according to CPU resources.
coarse=baselineOpts;
coarse.Nz=30; coarse.Nv=40; coarse.Ny=10;
coarse.NtInt=40; coarse.Nconv=30;
coarse.NuRel=10; coarse.NtRel=40;
coarse.warnOnMass=false;
fine=baselineOpts;
fine.Nz=90; fine.Nv=120; fine.Ny=28;
fine.NtInt=120; fine.Nconv=90;
fine.NuRel=28; fine.NtRel=90;
fine.warnOnMass=false;

levels={'Coarse';'Baseline';'Fine'};
reps=cell(3,1);
runtime=zeros(3,1);
for k=1:3
    t0=tic;
    if k==1
        [reps{k},~]=section65_analytical_reference( ...
            policy.RI,policy.M,par,coarse,[],true);
    elseif k==2
        reps{k}=baselineRep;
    else
        [reps{k},~]=section65_analytical_reference( ...
            policy.RI,policy.M,par,fine,[],true);
    end
    runtime(k)=toc(t0);
    if k==2, runtime(k)=baselineRuntime; end
    fprintf('  %s grid: %d states, %d time points; runtime %.1fs\n', ...
        levels{k}, numel(reps{k}.stationary.piVec), ...
        numel(reps{k}.reliability.t),runtime(k));
end
% Reinterpolate all curves to the SAME time interval before comparing.
tmax=min(cellfun(@(r) r.reliability.t(end),reps));
commonTime=linspace(0,tmax,120)';
fineNew=interp1(reps{3}.reliability.t,reps{3}.reliability.RfromNew, ...
    commonTime,'linear');
fineStat=interp1(reps{3}.reliability.t,reps{3}.reliability.Rstationary, ...
    commonTime,'linear');
cost=zeros(3,1); avail=zeros(3,1); mttf=zeros(3,1);
maxRnewDiff=zeros(3,1); maxRstatDiff=zeros(3,1);
rawMassError=zeros(3,1); mttfRcond=zeros(3,1);
for k=1:3
    r=reps{k};
    cost(k)=r.costRate;
    avail(k)=r.averageAvailability;
    mttf(k)=r.MTTF_stationary;
    rawMassError(k)=r.diagnostics.maxRawTransitionMassError;
    mttfRcond(k)=r.diagnostics.minMTTFLinearMatrixReciprocalCondition;
    newCurve=interp1(r.reliability.t,r.reliability.RfromNew,commonTime,'linear');
    statCurve=interp1(r.reliability.t,r.reliability.Rstationary,commonTime,'linear');
    maxRnewDiff(k)=max(abs(newCurve-fineNew));
    maxRstatDiff(k)=max(abs(statCurve-fineStat));
end
outTable=table(levels,cost,avail,mttf,maxRnewDiff,maxRstatDiff, ...
    rawMassError,mttfRcond,runtime, ...
    'VariableNames',{'Grid','CostRate','Availability','MTTFstationary', ...
    'MaxDiffRnewVsFine','MaxDiffRstatVsFine','MaxRawTransitionMassError', ...
    'MTTFLinearSystemRcond','Seconds'});
end


%% ========================================================================
%  Ready-to-paste LaTeX tables, matching Section 6.5 text.
% ========================================================================
function writeLatexSingleCycle(T,sValues,outDir)
fid=fopen(fullfile(outDir,'Table_section65_single_cycle.tex'),'w');
if fid<0, error('Cannot write single-cycle LaTeX table.'); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n');
fprintf(fid,'\\caption{Analytical and Monte Carlo comparison of single-cycle measures under different initial states}\n');
fprintf(fid,'\\label{tab:single_cycle_validation}\n');
fprintf(fid,'\\begin{tabular}{ccccc}\n\\toprule\n');
fprintf(fid,'$s$ & Measure & Analytical & Monte Carlo & Error (\\%%) \\\\\n\\midrule\n');
labels={'$\mu(s)$','$P_{\mathrm{pm}}(s)$', ...
    '$P_{\mathrm{pr}}(s)$','$P_{\mathrm{cr}}(s)$','$d(s)$'};
for i=1:numel(sValues)
    for j=1:5
        k=(i-1)*5+j;
        if j==1
            fprintf(fid,'\\multirow{5}{*}{$%.4f$}',sValues(i));
        else
            fprintf(fid,' ');
        end
        fprintf(fid,' & %s & %.6f & %.6f & %.3f \\\\\n', ...
            labels{j},T.Analytical(k),T.MonteCarlo(k),T.RelativeErrorPct(k));
    end
    if i<numel(sValues), fprintf(fid,'\\midrule\n'); end
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function writeLatexSystem(T,outDir)
fid=fopen(fullfile(outDir,'Table_section65_system_performance.tex'),'w');
if fid<0, error('Cannot write system-performance LaTeX table.'); end
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'\\begin{table}[htbp]\n\\centering\n');
fprintf(fid,'\\caption{Comparison of analytical and Monte Carlo system performance}\n');
fprintf(fid,'\\label{tab:system_validation}\n');
fprintf(fid,'\\begin{tabular}{lccc}\n\\toprule\n');
fprintf(fid,'Measure & Analytical & Monte Carlo & Error (\\%%) \\\\\n\\midrule\n');
labels={'$C_{\infty}$','$\bar A$','$\mathrm{MTTF}_{\pi}$', ...
    '$\mathrm{MTTF}_{0}$','$E_{\pi}[\Delta T]$', ...
    '$E_{\pi}[N_{\mathrm{pm}}]$','$E_{\pi}[N_{\mathrm{pr}}]$', ...
    '$E_{\pi}[N_{\mathrm{cr}}]$','$E_{\pi}[D]$'};
for k=1:height(T)
    fprintf(fid,'%s & %.6f & %.6f & %.3f \\\\\n', ...
        labels{k},T.Analytical(k),T.MonteCarlo(k),T.RelativeErrorPct(k));
end
fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end
