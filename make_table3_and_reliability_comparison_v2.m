%MAKE_TABLE3_AND_RELIABILITY_COMPARISON_V2
% Load the completed results of:
%   1) Policy 1,
%   2) Policy 2,
%   3) the proposed hybrid policy,
%
% and then:
%   - output Table 3;
%   - draw the stationary-start reliability comparison;
%   - add the Appendix-A reliability curve obtained by starting from the
%     stationary distribution of the proposed policy and applying no
%     inspection or maintenance afterwards.
%
% The proposed policy is NOT rerun by this script. Its optimal result and
% reliability curve are loaded directly from the existing MAT file.
%
% Required MAT files in the current MATLAB folder:
%   shock_cbm_results_v5_coarse_timed.mat
%   policy1_results_v2_coarse_timed.mat
%   policy2_results_v2_coarse_timed.mat
%
% How to run:
%   >> make_table3_and_reliability_comparison_v2
%
% Main output variables:
%   table3Numeric
%   table3Display
%   reliabilityComparison
%   comparisonResults
%   timingLog

clearvars;
clc;
close all;

scriptTimer = tic;
timingNames = {};
timingSeconds = [];

fprintf('\n============================================================\n');
fprintf(' Table 3 and stationary reliability comparison\n');
fprintf('============================================================\n');

%% --------------------- Input files ---------------------------
proposedFile = 'shock_cbm_results_v5_coarse_timed.mat';
policy1File  = 'policy1_results_v2_coarse_timed.mat';
policy2File  = 'policy2_results_v2_coarse_timed.mat';

assertFileExists(proposedFile);
assertFileExists(policy1File);
assertFileExists(policy2File);

%% --------------------- Load completed results ----------------
fprintf('\n[1] Loading the three completed policy results...\n');

functionTimer = tic;

SH = load(proposedFile,'results');
S1 = load(policy1File,'results');
S2 = load(policy2File,'results');

timingNames{end+1,1} = 'load completed policy results';
timingSeconds(end+1,1) = toc(functionTimer);

if ~isfield(SH,'results') || ~isfield(SH.results,'finalBest')
    error('%s does not contain results.finalBest.',proposedFile);
end

if ~isfield(S1,'results') || ~isfield(S1.results,'finalBest')
    error('%s does not contain results.finalBest.',policy1File);
end

if ~isfield(S2,'results') || ~isfield(S2.results,'finalBest')
    error('%s does not contain results.finalBest.',policy2File);
end

resultsProposed = SH.results;
resultsPolicy1  = S1.results;
resultsPolicy2  = S2.results;

repProposed = resultsProposed.finalBest;
repPolicy1  = resultsPolicy1.finalBest;
repPolicy2  = resultsPolicy2.finalBest;

checkResultStructure(repProposed,'Proposed policy');
checkResultStructure(repPolicy1,'Policy 1');
checkResultStructure(repPolicy2,'Policy 2');

checkParameterConsistency( ...
    repProposed.parameters, ...
    repPolicy1.parameters, ...
    repPolicy2.parameters);

par = repProposed.parameters;

fprintf('\nSaved reliability-tail diagnostics:\n');
fprintf('  Proposed policy: %.3e\n', ...
    repProposed.reliability.Rstationary(end));
fprintf('  Policy 1:        %.3e\n', ...
    repPolicy1.reliability.Rstationary(end));
fprintf('  Policy 2:        %.3e\n', ...
    repPolicy2.reliability.Rstationary(end));

if repProposed.reliability.Rstationary(end) > 1e-5
    warning(['The saved proposed-policy reliability curve ends at %.3g. ', ...
        'The comparison interval is limited by the existing result.'], ...
        repProposed.reliability.Rstationary(end));
end

%% --------------------- Table 3 -------------------------------
fprintf('\n[2] Constructing Table 3...\n');

functionTimer = tic;

Policy = [ ...
    "Policy 1"; ...
    "Policy 2"; ...
    "Proposed policy"];

RIStar = [ ...
    repPolicy1.policy.RI; ...
    NaN; ...
    repProposed.policy.RI];

MStar = [ ...
    repPolicy1.policy.M; ...
    repPolicy2.policy.M; ...
    repProposed.policy.M];

CostRate = [ ...
    repPolicy1.costRate; ...
    repPolicy2.costRate; ...
    repProposed.costRate];

AverageAvailability = [ ...
    repPolicy1.averageAvailability; ...
    repPolicy2.averageAvailability; ...
    repProposed.averageAvailability];

MTTF = [ ...
    repPolicy1.MTTF_stationary; ...
    repPolicy2.MTTF_stationary; ...
    repProposed.MTTF_stationary];

table3Numeric = table( ...
    Policy, RIStar, MStar, CostRate, AverageAvailability, MTTF, ...
    'VariableNames',{ ...
    'Policy','RIStar','MStar','CostRate', ...
    'AverageAvailability','MTTF'});

% Display table with "--" for the Policy 2 reliability threshold.
RIStarDisplay = strings(3,1);
RIStarDisplay(1) = sprintf('%.6f',RIStar(1));
RIStarDisplay(2) = "--";
RIStarDisplay(3) = sprintf('%.6f',RIStar(3));

table3Display = table( ...
    Policy, ...
    RIStarDisplay, ...
    compose('%.6f',MStar), ...
    compose('%.6f',CostRate), ...
    compose('%.6f',AverageAvailability), ...
    compose('%.6f',MTTF), ...
    'VariableNames',{ ...
    'Policy','RIStar','MStar','CostRate', ...
    'AverageAvailability','MTTF'});

timingNames{end+1,1} = 'construct Table 3';
timingSeconds(end+1,1) = toc(functionTimer);

fprintf('\n---------------- Table 3 ----------------\n');
disp(table3Display);

writetable(table3Display,'Table3_results_v2.csv');

%% --------------------- Common reliability grid --------------
fprintf('\n[3] Preparing a common reliability time grid...\n');

% User-adjustable plotting settings.
plotOpts.Nt = 401;

% Use [] to display the entire common interval available in all three
% completed policy results. A value such as 60 displays [0,60].
plotOpts.tMaximum = [];

% Clean very small numerical increases in reliability curves.
plotOpts.enforceMonotonicity = true;

tEndAvailable = min([ ...
    repProposed.reliability.t(end), ...
    repPolicy1.reliability.t(end), ...
    repPolicy2.reliability.t(end)]);

if isempty(plotOpts.tMaximum)
    tEnd = tEndAvailable;
else
    tEnd = min(plotOpts.tMaximum,tEndAvailable);
end

if ~(isfinite(tEnd) && tEnd > 0)
    error('The common reliability time interval is invalid.');
end

tCommon = linspace(0,tEnd,plotOpts.Nt)';

RProposed = interpolateReliabilityCurve( ...
    repProposed.reliability.t, ...
    repProposed.reliability.Rstationary, ...
    tCommon, ...
    plotOpts.enforceMonotonicity);

RPolicy1 = interpolateReliabilityCurve( ...
    repPolicy1.reliability.t, ...
    repPolicy1.reliability.Rstationary, ...
    tCommon, ...
    plotOpts.enforceMonotonicity);

RPolicy2 = interpolateReliabilityCurve( ...
    repPolicy2.reliability.t, ...
    repPolicy2.reliability.Rstationary, ...
    tCommon, ...
    plotOpts.enforceMonotonicity);

%% --------------------- Appendix-A reliability ---------------
fprintf('\n[4] Computing Appendix-A reliability without actions...\n');

% The initial state distribution is the stationary post-maintenance
% distribution of the proposed policy.
noActionOpts.NuCompound   = 1201;
noActionOpts.NfftCompound = 8192;

functionTimer = tic;

[RNoAction,RNoActionByState,noActionDiagnostics] = ...
    reliabilityWithoutInspectionAndMaintenance( ...
    tCommon, ...
    repProposed.stationary, ...
    par, ...
    noActionOpts);

timingNames{end+1,1} = ...
    'Appendix-A reliability without inspection and maintenance';

timingSeconds(end+1,1) = toc(functionTimer);

%% --------------------- Reliability comparison table --------
reliabilityComparison = table( ...
    tCommon, ...
    RProposed, ...
    RPolicy1, ...
    RPolicy2, ...
    RNoAction, ...
    'VariableNames',{ ...
    'Time','ProposedPolicy','Policy1','Policy2', ...
    'WithoutInspectionAndMaintenance'});

writetable( ...
    reliabilityComparison, ...
    'stationary_reliability_comparison_v2.csv');

%% --------------------- Reliability comparison plot ----------
fprintf('\n[5] Drawing the four reliability curves...\n');

functionTimer = tic;

fig = figure( ...
    'Name','Stationary reliability comparison', ...
    'Color','w', ...
    'Position',[100,80,860,620]);

ax = axes('Parent',fig);
hold(ax,'on');

plot(ax,tCommon,RProposed, ...
    '-','LineWidth',2.1);

plot(ax,tCommon,RPolicy1, ...
    '--','LineWidth',1.9);

plot(ax,tCommon,RPolicy2, ...
    '-.','LineWidth',1.9);

plot(ax,tCommon,RNoAction, ...
    ':','LineWidth',2.3);

xlabel(ax,'Time');
ylabel(ax,'Reliability');

title(ax, ...
    'Reliability functions starting from stationary distributions');

legend(ax, ...
    'Proposed policy', ...
    'Policy 1', ...
    'Policy 2', ...
    'Without inspection and maintenance', ...
    'Location','best');

xlim(ax,[0,tEnd]);
ylim(ax,[0,1]);

grid(ax,'on');
box(ax,'on');

set(ax, ...
    'FontSize',11, ...
    'LineWidth',1);

savefig(fig,'stationary_reliability_comparison_v2.fig');

try
    exportgraphics( ...
        fig, ...
        'stationary_reliability_comparison_v2.png', ...
        'Resolution',300);

    exportgraphics( ...
        fig, ...
        'stationary_reliability_comparison_v2.pdf', ...
        'ContentType','vector');
catch
    print( ...
        fig, ...
        'stationary_reliability_comparison_v2.png', ...
        '-dpng','-r300');

    print( ...
        fig, ...
        'stationary_reliability_comparison_v2.pdf', ...
        '-dpdf','-painters');
end

timingNames{end+1,1} = ...
    'plot stationary reliability comparison';

timingSeconds(end+1,1) = toc(functionTimer);

fprintf(['Saved reliability comparison to ', ...
    'stationary_reliability_comparison_v2.fig, .png, and .pdf\n']);

%% --------------------- Save comparison result ---------------
comparisonResults = struct();

comparisonResults.files = struct( ...
    'proposed',proposedFile, ...
    'policy1',policy1File, ...
    'policy2',policy2File);

comparisonResults.table3Numeric = table3Numeric;
comparisonResults.table3Display = table3Display;

comparisonResults.reliability = struct( ...
    't',tCommon, ...
    'proposedPolicy',RProposed, ...
    'policy1',RPolicy1, ...
    'policy2',RPolicy2, ...
    'withoutInspectionAndMaintenance',RNoAction, ...
    'withoutActionsByInitialState',RNoActionByState);

comparisonResults.noActionDiagnostics = ...
    noActionDiagnostics;

comparisonResults.parameters = par;
comparisonResults.plotOptions = plotOpts;

timingNames{end+1,1} = 'total script time before save';
timingSeconds(end+1,1) = toc(scriptTimer);

timingLog = table( ...
    timingNames,timingSeconds, ...
    'VariableNames',{'FunctionOrStep','Seconds'});

comparisonResults.timingLog = timingLog;

save( ...
    'table3_and_reliability_comparison_v2.mat', ...
    'comparisonResults', ...
    'table3Numeric', ...
    'table3Display', ...
    'reliabilityComparison', ...
    'timingLog');

fprintf('\n---------------- Runtime log ----------------\n');
disp(timingLog);

fprintf('\nSaved outputs:\n');
fprintf('  Table3_results_v2.csv\n');
fprintf('  stationary_reliability_comparison_v2.csv\n');
fprintf('  table3_and_reliability_comparison_v2.mat\n');
fprintf('  stationary_reliability_comparison_v2.fig\n');
fprintf('  stationary_reliability_comparison_v2.png\n');
fprintf('  stationary_reliability_comparison_v2.pdf\n');
fprintf('Done.\n');

%% ========================================================================
%  Appendix-A reliability without inspection and maintenance
% ========================================================================
function [Rstationary,RbyState,diagnostics] = ...
    reliabilityWithoutInspectionAndMaintenance( ...
    tGrid,stationary,par,opts)

if ~isfield(stationary,'z') || ...
        ~isfield(stationary,'piVec')

    error(['The proposed-policy stationary result must contain ', ...
        'stationary.z and stationary.piVec.']);
end

states = [0; stationary.z(:)];
piVec = stationary.piVec(:);

if numel(states) ~= numel(piVec)
    error('The proposed-policy states and piVec have different lengths.');
end

piVec = max(real(piVec),0);
piVec = piVec/sum(piVec);

compound = buildCompoundCdfAtTimes( ...
    tGrid,par,opts);

Ns = numel(states);
Nt = numel(tGrid);

RbyState = zeros(Nt,Ns);

for is = 1:Ns
    margin = par.L-states(is);

    RbyState(:,is) = interp1( ...
        compound.uGrid, ...
        compound.cdf.', ...
        margin, ...
        'linear').';
end

RbyState = min(max(real(RbyState),0),1);
RbyState = cummin(RbyState,1);
RbyState(1,:) = 1;

Rstationary = RbyState*piVec;
Rstationary = min(max(real(Rstationary),0),1);
Rstationary = cummin(Rstationary);
Rstationary(1) = 1;

diagnostics.fftUpperBound = compound.fftUpperBound;
diagnostics.fftStep = compound.fftStep;
diagnostics.minimumRetainedMass = ...
    min(compound.totalRetainedMass);
diagnostics.maximumRetainedMass = ...
    max(compound.totalRetainedMass);
diagnostics.maximumUpperTailMass = ...
    max(compound.upperTailMass);
end

function compound = buildCompoundCdfAtTimes(tGrid,par,opts)
tGrid = tGrid(:);

Nt = numel(tGrid);
Nu = opts.NuCompound;
Nfft = 2^nextpow2(opts.NfftCompound);

tMaximum = max(tGrid);

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

shockMass = diff(gammaCdfRate( ...
    fftEdges,par.alphaS,par.betaS));

shockMass = max(real(shockMass),0);
shockTransform = fft(shockMass,Nfft);

uGrid = linspace(0,par.L,Nu);
cdfTable = zeros(Nt,Nu);
retainedMass = zeros(Nt,1);
upperTailMass = zeros(Nt,1);

nTail = max(10,round(0.02*Nfft));

for it = 1:Nt
    t = tGrid(it);

    if t <= 0
        cdfTable(it,:) = 1;
        cdfTable(it,1) = 0;
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
    upperTailMass(it) = ...
        sum(totalMass(end-nTail+1:end));

    cdfAtEdges = [0,cumsum(totalMass)];

    cdfTable(it,:) = interp1( ...
        fftEdges,cdfAtEdges,uGrid,'linear',0);

    cdfTable(it,:) = ...
        min(max(cdfTable(it,:),0),1);
end

for iu = 1:Nu
    cdfTable(:,iu) = cummin(cdfTable(:,iu));
end

compound.tGrid = tGrid;
compound.uGrid = uGrid;
compound.cdf = cdfTable;
compound.fftUpperBound = fftUpperBound;
compound.fftStep = du;
compound.totalRetainedMass = retainedMass;
compound.upperTailMass = upperTailMass;
end

%% ========================================================================
%  Reliability interpolation
% ========================================================================
function Rcommon = interpolateReliabilityCurve( ...
    tOriginal,ROriginal,tCommon,enforceMonotonicity)

tOriginal = tOriginal(:);
ROriginal = ROriginal(:);

valid = ...
    isfinite(tOriginal) & ...
    isfinite(ROriginal);

tOriginal = tOriginal(valid);
ROriginal = ROriginal(valid);

if numel(tOriginal) < 2
    error('A saved reliability curve contains fewer than two valid points.');
end

[tOriginal,uniqueIndex] = unique(tOriginal,'stable');
ROriginal = ROriginal(uniqueIndex);

Rcommon = interp1( ...
    tOriginal,ROriginal,tCommon,'pchip');

Rcommon = min(max(real(Rcommon),0),1);

if enforceMonotonicity
    Rcommon = cummin(Rcommon);
end

Rcommon(1) = 1;
end

%% ========================================================================
%  Result and parameter validation
% ========================================================================
function checkResultStructure(rep,name)
required = { ...
    'policy', ...
    'costRate', ...
    'averageAvailability', ...
    'MTTF_stationary', ...
    'reliability', ...
    'stationary', ...
    'parameters'};

for k = 1:numel(required)
    if ~isfield(rep,required{k})
        error('%s result lacks field %s.',name,required{k});
    end
end

if ~isfield(rep.reliability,'t') || ...
        ~isfield(rep.reliability,'Rstationary') || ...
        isempty(rep.reliability.t) || ...
        isempty(rep.reliability.Rstationary)

    error('%s does not contain a saved stationary reliability curve.',name);
end

if ~isfield(rep.stationary,'z') || ...
        ~isfield(rep.stationary,'piVec')

    error('%s does not contain stationary.z and stationary.piVec.',name);
end
end

function checkParameterConsistency(pH,p1,p2)
names = { ...
    'alpha','beta','lambda','alphaS','betaS', ...
    'p','q','L','Ci','Cm','Cp','Cc','Cd'};

for k = 1:numel(names)
    name = names{k};

    if ~isfield(pH,name) || ...
            ~isfield(p1,name) || ...
            ~isfield(p2,name)

        error('Parameter %s is missing from at least one result.',name);
    end

    values = [pH.(name),p1.(name),p2.(name)];
    tolerance = 1e-10*max(1,max(abs(values)));

    if max(values)-min(values) > tolerance
        error(['Parameter mismatch for %s: proposed=%.12g, ', ...
            'Policy1=%.12g, Policy2=%.12g.'], ...
            name,values(1),values(2),values(3));
    end
end
end

function assertFileExists(fileName)
if exist(fileName,'file') ~= 2
    error('Cannot find required file: %s',fileName);
end
end

%% ========================================================================
%  Gamma CDF with rate parameterization
% ========================================================================
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
