function Run_GaAs_Grid_EM()
% RUN_GAAS_GRID_EM
% Sparse shocks at randomly selected existing 250-hour observation times + semi-synthetic
% GaAs laser mixed-degradation observations + SHAPE/RATE Gamma shock EM.
%
% DATA: GaAsLaser.csv, 15 real laser degradation histories (250 h grid).
% ARRIVALS: independent Bernoulli selection at each fixed 250-hour observation endpoint.
% SHOCK MAGNITUDES: Gamma(shape=5.0, RATE=7.5).
%   The preset grid shocks are fixed in GaAs_FixedInjection.csv for reproducibility.
% EXPERIMENTAL DEGRADATION AT SHOCK TIMES:
%   Shock times are grid points, so original experimental degradation readings
%   are used DIRECTLY. No Gamma-bridge interpolation is necessary.
% OBSERVATIONS: EM uses only increments of the 250-h monitoring interval
%   ENDING at a selected shock point (U,Y), INCLUDING post-threshold points.
%   U_i = [X(t_i)-X(t_i-250h)] + W_i, with Y_i=1.
%   Inter-shock times are recorded separately for audit and NOT used as Y.
%   All 15*16=240 intervals remain eligible; L=10 is NOT an EM censoring rule.
%   This is an uncensored semi-synthetic sensitivity experiment; observations
%   beyond L must not be described as operating data after physical failure.
% EM: uses ONLY (U_i,Y_i) and known X Gamma parameters, never hidden W_i.
% TIME: 1 unit=250 h; Bernoulli p=0.25 at each grid point; average 1 shock/1000 h.
% WARNING: A Bernoulli grid process is NOT an exact continuous-time HPP.
%          BOTH policy-optimization groups use ESTIMATED lambda=N_shocks/exposure.
%          shockProb=0.25 controls only semi-synthetic data generation.
% SHOCK SCHEDULE: provided by a reproducible fixed injection CSV; do not rerandomize when estimating.
%
% Outputs generated next to this script:
%   GaAs_Grid_ShockSchedule_All.csv         (synthetic truth/audit)
%   GaAs_Grid_EM_Input_ONLY.csv             (ONLY (U,Y) + timestamps)
%   GaAs_Grid_Truth_DO_NOT_FEED_EM.csv      (validation only)
%   GaAs_Grid_GridMeasurements.csv          (fixed-grid comparison for plots)
%   GaAs_Grid_EM_Iterations.csv
%   GaAs_Grid_EM_Truth_vs_Estimate.csv
%   GaAs_Grid_EM_Summary.csv                (read by policy optimization)
%   GaAs_Grid_Unit101_Comparison.png / .fig
%
% MATLAB without Statistics Toolbox: uses randg, integral, fzero, psi.

rootDir = fileparts(mfilename('fullpath'));
sourceFile = fullfile(rootDir,'GaAsLaser.csv');
if ~isfile(sourceFile)
    error('Missing GaAsLaser.csv in %s',rootDir);
end

%% Parameters shared with policy optimization
H = 250;               % hours per model time unit
shockProb = 0.25;      % injection design probability, NOT assumed true HPP lambda
alphaX = 7.426;        % Gamma-process shape RATE per model time unit
betaX = 1/0.0674;     % Gamma RATE (not scale)
shapeS_True = 5.0;   % fixed semi-synthetic example, SHAPE
scaleS_True = 1/7.5; % Gamma SCALE used for generation; RATE = 7.5
rateS_True = 1/scaleS_True;
L = 10;
rngSeed = 20261122;
%% Injection record is fixed in GaAs_FixedInjection.csv (NumPy PCG64, seed 20261122).
% rngSeed is recorded for reproducibility; MATLAB's twister is not used.


fprintf('\n========================================================\n');
fprintf(' GaAs FIXED random-grid shocks: p=%.2f, shock shape=%.2f\n',...
    shockProb,shapeS_True);
fprintf(' Time unit: 250 h; mean grid-event interval: %.0f h\n',H/shockProb);
fprintf(' NO truncation of shock injection or EM data at L=%.1f\n',L);
fprintf('========================================================\n');

%% Read ORIGINAL experimental trajectories and fixed pre-generated random shocks
Raw = readtable(sourceFile,'VariableNamingRule','preserve');
required = {'Percent Increase','Unit Number','Hours'};
assert(all(ismember(required,Raw.Properties.VariableNames)), ...
    'GaAsLaser.csv must have Percent Increase, Unit Number, Hours.');
idRaw = double(Raw.('Unit Number'));
tRaw = double(Raw.Hours);
xRaw = double(Raw.('Percent Increase'));
[~,sortOrder] = sortrows([idRaw,tRaw],[1,2]);
idRaw = idRaw(sortOrder); tRaw = tRaw(sortOrder); xRaw = xRaw(sortOrder);
unitIDs = unique(idRaw,'stable');
assert(numel(unitIDs)==15,'Expect 15 units in original GaAs data.');
% Fixed schedule: all 240 endpoints are specified explicitly.
% This CSV is used only to construct semi-synthetic U. The EM step never
% reads TrueShockDamage, only the generated (U,Y) table below.
nIntervals = 16; % 0 to 4000 hours, sampled every 250 hours
presetFile=fullfile(rootDir,'GaAs_FixedInjection.csv');
assert(isfile(presetFile),'Missing GaAs_FixedInjection.csv');
Preset=readtable(presetFile);
assert(height(Preset)==numel(unitIDs)*nIntervals, ...
    'Injection table should contain 240 monitoring intervals.');
assert(all(abs(Preset.AlphaS_True-shapeS_True)<1e-12) && ...
       all(abs(Preset.BetaS_Rate_True-rateS_True)<1e-12), ...
    'Injection CSV parameters differ from shapeS_True / rateS_True.');
shockMask=false(numel(unitIDs),nIntervals);
shockMagnitude=zeros(numel(unitIDs),nIntervals);
for iu=1:numel(unitIDs)
    for jj=1:nIntervals
        rr=(Preset.UnitID==unitIDs(iu) & ...
            Preset.ShockTimeHours==jj*H);
        assert(nnz(rr)==1,'Missing or duplicate preset shock slot.');
        shockMask(iu,jj)=logical(Preset.ShockInjected(rr));
        shockMagnitude(iu,jj)=Preset.TrueShockDamage(rr);
        assert(shockMagnitude(iu,jj)>=0 && ...
            (shockMask(iu,jj) || shockMagnitude(iu,jj)==0), ...
            'Invalid shock value in preset file.');
    end
end
fprintf('Preset random shocks loaded: %d of %d intervals.\n', ...
        nnz(shockMask),numel(shockMask));

% Each row of scheduleRows corresponds to ONE selected fixed-grid shock.
% True shock values and original X(t) are audit-only and NEVER input to EM.
% Columns: unit, shock time, INTERSHOCK gap hours, W, continuous X(t),
% cum W, total S(t), U from MOST RECENT 250h interval, Y=1,
% first-threshold crossing for diagnosis, all shocks USED for EM.
scheduleRows = zeros(0,11);
emRows       = zeros(0,6);
truthRows    = zeros(0,6);
gridRows     = zeros(0,6);
firstFailures = 0;
allExposureUnits = 0;
observedExposureUnits = 0; % all 240 recorded grid intervals: NO failure censoring

for iu = 1:numel(unitIDs)
    uid = unitIDs(iu);
    ii = find(idRaw == uid);
    hours = tRaw(ii); x = xRaw(ii);
    assert(numel(hours)==17 && hours(1)==0 && ...
        all(abs(diff(hours)-H)<1e-9) && ...
        x(1)==0 && all(diff(x)>0), ...
        'Unexpected source grid/degradation for unit %d.',uid);
    allExposureUnits = allExposureUnits + (hours(end)-hours(1))/H;

    % The complete injected path, used for visualization/audit.
    cumW = 0;
    failed = false;
    firstFailureTime = Inf;
    previousShockT = 0;
    gridRows(end+1,:) = [uid,hours(1),x(1),x(1),0,1]; %#ok<AGROW>

    for j = 2:numel(hours)
        t = hours(j);
        selected = shockMask(iu,j-1);  % fixed beforehand; independent of W size draws
        W = 0;
        if selected
            W = shockMagnitude(iu,j-1);
            cumW = cumW + W;
        end
        mixedLevel = x(j)+cumW;
        failureNow = ~failed && (mixedLevel >= L);
        isEligible = true;  % use ALL existing monitoring points even after L
        observedExposureUnits = observedExposureUnits + 1;

        if selected
            % Each selected shock belongs to the immediately preceding
            % observed 250-hour interval, not the interval since the last
            % shock. Thus the EM continuous component has fixed Y=1.
            % The inter-shock gap remains available for auditing only.
            previousMonitoringT = hours(j-1);
            Y = (t-previousMonitoringT)/H;
            deltaX = x(j)-x(j-1);
            U = deltaX + W;
            assert(abs(Y-1)<1e-12 && U>0, ...
                'Invalid 250-hour EM increment for unit %d.',uid);
            scheduleRows(end+1,:) = [uid,t,t-previousShockT,W,x(j), ...
                cumW,mixedLevel,U,Y,double(failureNow),double(isEligible)]; %#ok<AGROW>
            % No post-failure exclusion: EVERY injected shock is an EM sample.
            emRows(end+1,:) = [uid,previousMonitoringT,t,t-previousMonitoringT,Y,U]; %#ok<AGROW>
            truthRows(end+1,:) = [uid,t,deltaX,W,U,double(failureNow)]; %#ok<AGROW>
            previousShockT=t;
        end

        if failureNow
            failed = true;
            firstFailureTime = t;
            firstFailures = firstFailures+1;
        end
        % Failure flags are for reporting only. Even after L, readings
        % remain in the semi-synthetic sample and the full-trajectory plot.
        gridRows(end+1,:) = [uid,t,x(j),mixedLevel,double(selected), ...
                            double(t<=firstFailureTime)]; %#ok<AGROW>
    end
end

Schedule = array2table(scheduleRows,'VariableNames', ...
    {'UnitID','ShockTimeHours','InterarrivalHours','TrueShockDamage', ...
    'TrueBackgroundDegradation','CumulativeShockDamage', ...
    'ObservedDegradation','U_MixedIncrement','Y_TimeUnits', ...
    'FirstDetectedFailure','UsedForEM'});
EMInput = array2table(emRows,'VariableNames', ...
    {'UnitID','PreviousHours','ShockTimeHours', ...
     'IntervalHours','Y_TimeUnits','U_MixedIncrement'});
Truth = array2table(truthRows,'VariableNames', ...
    {'UnitID','ShockTimeHours','TrueBackgroundIncrement', ...
     'TrueShockDamage','U_MixedIncrement','FirstDetectedFailure'});
Grid = array2table(gridRows,'VariableNames', ...
    {'UnitID','Hours','OriginalDegradation','MixedDegradation', ...
     'ShockAtThisPoint','ValidUntilFirstDetectedFailure'});

writetable(Schedule,fullfile(rootDir,'GaAs_Grid_ShockSchedule_All.csv'));
writetable(EMInput,fullfile(rootDir,'GaAs_Grid_EM_Input_ONLY.csv'));
writetable(Truth,fullfile(rootDir,'GaAs_Grid_Truth_DO_NOT_FEED_EM.csv'));
writetable(Grid,fullfile(rootDir,'GaAs_Grid_GridMeasurements.csv'));

nAll = height(Schedule);
nEM = height(EMInput);
realizedLambda = nAll / allExposureUnits;  % full simulated horizon; audit only
% Estimate event probability from ALL 240 recorded 250h grid intervals.
% Under the CURRENT injection mechanism: Bernoulli event at 250h endpoint.
% Hence lambdaHatFrequency = nEM / observedExposureUnits is the empirical
% event FREQUENCY used to match the approximate downstream HPP policy model.
% With NO truncation, nEM=nAll and observedExposureUnits=allExposureUnits.
% If a HPP is fitted to the grid event-presence indicators instead,
% P(N(1)>=1)=1-exp(-lambda), giving lambdaHatHPP=-log(1-pHat).
% These are distinct estimands: do NOT label the Bernoulli injection as HPP.
assert(observedExposureUnits==allExposureUnits && ...
    nEM==nAll && observedExposureUnits==numel(unitIDs)*nIntervals, ...
    'NO-CENSORING INVARIANT FAILED: all selected shocks must enter EM.');
shockProbMLE = nEM / observedExposureUnits;
lambdaHatFrequency = shockProbMLE;
if shockProbMLE < 1
    lambdaHatHPP = -log1p(-shockProbMLE);
else
    lambdaHatHPP = Inf;
end
fprintf('Random grid shocks: %d across %.0f candidate 250-hour intervals.\n', ...
    nAll,allExposureUnits);
fprintf('Mean event frequency = %.6f vs Bernoulli probability = %.6f.\n', ...
    realizedLambda,shockProb);
fprintf('EM: %d selected 250-hour increments (FULL trajectory, NO truncation).\n',nEM);
fprintf('All %d unit trajectories retained to 4000 h; first threshold crossings: %d.\n', ...
    numel(unitIDs),firstFailures);
fprintf('Full-grid exposure: %g recorded 250-hour intervals.\n',observedExposureUnits);
fprintf('Estimated grid event probability: %d / %g = %.6f.\n', ...
    nEM,observedExposureUnits,shockProbMLE);
fprintf('Estimated frequency-matched lambda = %.6f per 250 h.\n', ...
    lambdaHatFrequency);
fprintf('HPP intensity implied by occupancy indicators = %.6f (alternative).\n', ...
    lambdaHatHPP);
if nEM<10
    error('Too few EM increments for stable estimation. Try another seed.');
end

%% EM: only U_i and Y_i are used, hidden true shock sizes are excluded.
U = double(EMInput.U_MixedIncrement(:));
Y = double(EMInput.Y_TimeUnits(:));
kX = alphaX*Y;
assert(all(U>0) && all(Y>0) && all(isfinite(U)), ...
    'EM observation input invalid.');

Z = U-kX/betaX;
mu0 = mean(Z);
var0 = var(Z,0)-mean(kX)/betaX^2;
if mu0>0 && var0>0
    alphaS = mu0^2/var0;
    betaS = mu0/var0;
    initMode = 'moment';
else
    warning('Moment starting values invalid; using positive fallback.');
    alphaS=1.0;
    betaS=alphaS/max(mean(Z),0.1);
    initMode = 'fallback';
end

maxIter = 800;
tol=1e-8;
history=zeros(maxIter+1,4);
history(1,:)=[0,alphaS,betaS,NaN];
converged=false;
for iter=1:maxIter
    [A,B] = localEStep(U,kX,alphaS,betaS,betaX);
    h = log(mean(A))-mean(B);
    if ~(h>0 && isfinite(h))
        error('EM M-step target invalid at iteration %d.',iter);
    end
    newAlpha = fzero(@(aa)log(aa)-psi(aa)-h,[1e-9,1e8]);
    newBeta = newAlpha/mean(A);
    relChange=max(abs([newAlpha-alphaS,newBeta-betaS]) ./ ...
       max(1,abs([alphaS,betaS])));
    alphaS=newAlpha; betaS=newBeta;
    history(iter+1,:)=[iter,alphaS,betaS,relChange];
    if iter<=5 || mod(iter,25)==0
        fprintf('EM %3d: shape=%.6f, rate=%.6f, change=%.3e\n', ...
            iter,alphaS,betaS,relChange);
    end
    if relChange<tol
        converged=true;
        break;
    end
end
history=history(1:iter+1,:);
if ~converged
    warning('EM did NOT meet tolerance after %d iterations; inspect estimates.',maxIter);
end

History= array2table(history,'VariableNames', ...
    {'Iteration','AlphaS_Shape','BetaS_Rate','RelativeParameterChange'});
writetable(History,fullfile(rootDir,'GaAs_Grid_EM_Iterations.csv'));

paramLabels={'AlphaS_Shape';'BetaS_Rate';'ShockMean'};
truthVals=[shapeS_True;rateS_True;shapeS_True/rateS_True];
estVals=[alphaS;betaS;alphaS/betaS];
relativeErr=100*abs(estVals-truthVals)./truthVals;
Compare=table(paramLabels,truthVals,estVals,relativeErr, ...
    'VariableNames',{'Parameter','TrueValue','Estimate','AbsoluteRelativeErrorPct'});
writetable(Compare,fullfile(rootDir,'GaAs_Grid_EM_Truth_vs_Estimate.csv'));

Summary=table(realizedLambda,shockProb,shockProbMLE, ...
    lambdaHatFrequency,lambdaHatHPP, ...
    H,alphaX,betaX,L,shapeS_True,rateS_True,alphaS,betaS, ...
    nAll,nEM,firstFailures,allExposureUnits,observedExposureUnits, ...
    iter,double(converged),rngSeed, ...
    'VariableNames',{'Lambda_Realized_FullHorizon', ...
    'ShockProbability_Set','ShockProbability_Estimated', ...
    'Lambda_Estimated_Frequency','Lambda_HPP_Implied_From_Indicators', ...
    'HoursPerUnit','AlphaX','BetaX_Rate','L', ...
    'AlphaS_True','BetaS_True','AlphaS_EM','BetaS_EM', ...
    'NTotalShocks','NEMRows','NDetectedFailures', ...
    'TotalExposureTimeUnits','ObservedExposureTimeUnits', ...
    'EMIterations','EMConverged','Seed'});
writetable(Summary,fullfile(rootDir,'GaAs_Grid_EM_Summary.csv'));

fprintf('\nEM finished: alphaS = %.8f, betaS = %.8f (rate), %d iterations.\n', ...
    alphaS,betaS,iter);
fprintf('Relative errors: shape %.2f%%, rate %.2f%%.\n',relativeErr(1),relativeErr(2));
fprintf('Both optimizer groups use the SAME estimated lambda; only shock magnitude shape/rate differ.\n');

%% Plot the ORIGINAL and randomly shock-injected grid trajectories for laser 101
uid = 101;
raw101All = Grid(Grid.UnitID==uid,:);
shock101 = Schedule(Schedule.UnitID==uid,:);
fig=figure('Color','w','Position',[100,100,980,530]);
plot(raw101All.Hours,raw101All.OriginalDegradation,'-o', ...
    'LineWidth',1.6,'MarkerSize',4); hold on;
plot(raw101All.Hours,raw101All.MixedDegradation,'-s', ...
    'LineWidth',1.8,'MarkerSize',4);
if ~isempty(shock101)
    plot(shock101.ShockTimeHours,shock101.ObservedDegradation, ...
        'ko','MarkerFaceColor','k','MarkerSize',6);
else
    plot(nan,nan,'ko','MarkerFaceColor','k');
end
xlabel('Operating time (hours)');
ylabel('Increase in operating current (%)');
title('GaAs laser unit 101: original and shock-injected degradation');
legend('Original experimental degradation','Shock-injected degradation', ...
       'Selected shock times','Location','northwest');
grid on; box on;
savefig(fig,fullfile(rootDir,'GaAs_Grid_Unit101_Comparison.fig'));
print(fig,fullfile(rootDir,'GaAs_Grid_Unit101_Comparison.png'),'-dpng','-r220');

f2=figure('Color','w','Position',[120,120,940,630]);
subplot(2,1,1);
plot(history(:,1),history(:,2),'-o','LineWidth',1.5,'MarkerSize',3);hold on;
yline(shapeS_True,'--');ylabel('Shock shape \alpha_S');grid on;
subplot(2,1,2);
plot(history(:,1),history(:,3),'-o','LineWidth',1.5,'MarkerSize',3);hold on;
yline(rateS_True,'--');ylabel('Shock rate \beta_S');
xlabel('EM iteration');grid on;
savefig(f2,fullfile(rootDir,'GaAs_Grid_EM_Convergence.fig'));
print(f2,fullfile(rootDir,'GaAs_Grid_EM_Convergence.png'),'-dpng','-r220');

fprintf('Random grid shock injection + EM outputs saved to %s\n',rootDir);

end

function [A,B] = localEStep(U,kX,aS,bS,betaX)
% Conditional moments under U=X_Y+W for Gamma(shape, RATE) components.
% z=w/U, conditional unnormalized weight is
% z^(aS-1)*(1-z)^(kX-1)*exp((betaX-bS)*U*z).
A=zeros(size(U)); B=zeros(size(U));
for i=1:numel(U)
    tilt=(betaX-bS)*U(i);
    % Center the integrand near its mode to avoid underflow. Without this,
    % the exponential factor for long inter-shock intervals can make the
    % unnormalized integral smaller than the absolute quadrature tolerance.
    zRef=min(max((aS/bS)/U(i),1e-5),1-1e-5);
    if aS>1 && kX(i)>1
        zRef=fzero(@(z) (aS-1)/z-(kX(i)-1)/(1-z)+tilt, ...
                   [1e-12,1-1e-12]);
    end
    logShift=(aS-1)*log(zRef)+(kX(i)-1)*log1p(-zRef)+tilt*zRef;
    ansVec=integral(@(z)emMomentIntegrand(z,aS,kX(i),tilt,logShift), ...
        0,1,'ArrayValued',true,'RelTol',1e-9,'AbsTol',1e-12);
    assert(numel(ansVec)==3 && all(isfinite(ansVec)) && ansVec(1)>0, ...
        'E-step integration failed at observation %d.',i);
    A(i)=U(i)*ansVec(2)/ansVec(1);
    B(i)=log(U(i))+ansVec(3)/ansVec(1);
end
end

function out = emMomentIntegrand(z,a,k,tilt,logShift)
z=max(1e-15,min(1-1e-15,z));
logw=(a-1)*log(z)+(k-1)*log1p(-z)+tilt*z-logShift;
w=exp(logw);
out=[w;z*w;log(z)*w];
end
