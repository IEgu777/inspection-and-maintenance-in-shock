%PLOT_ALL_SENSITIVITY_50PCT
% Plot the completed +/-50% sensitivity-analysis results.
%
% This script does NOT rerun any optimization. It only reads the 12 saved
% sensitivity result files and produces a two-panel tornado-style figure.
%
% Expected result files:
%   sensitivity_Ci_50pct_results.mat
%   sensitivity_Cm_50pct_results.mat
%   sensitivity_Cp_50pct_results.mat
%   sensitivity_Cc_50pct_results.mat
%   sensitivity_Cd_50pct_results.mat
%   sensitivity_alpha_50pct_results.mat
%   sensitivity_beta_50pct_results.mat
%   sensitivity_lambda_50pct_results.mat
%   sensitivity_alphaS_50pct_results.mat
%   sensitivity_betaS_50pct_results.mat
%   sensitivity_p_50pct_results.mat
%   sensitivity_q_50pct_results.mat
%
% Baseline result file:
%   shock_cbm_results_v5_coarse_timed.mat
%
% Put this script in the same folder as the sensitivity result files, or
% change resultsDir below.

clearvars;
clc;
close all;

%% ----------------------- User settings -----------------------
thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

% Folder containing the 12 sensitivity result MAT files.
resultsDir = thisDir;

% If the baseline result MAT file is not available, enter the baseline
% optimized long-run cost rate here, e.g. baselineCostManual = 37.623;
% Keep NaN to read it automatically from the baseline MAT file.
baselineCostManual = NaN;

baselineFileName = 'shock_cbm_results_v5_coarse_timed.mat';

% Sort parameters within each panel by sensitivity magnitude.
sortByImpact = true;

% Figure size.
figurePosition = [2 2 23 10];   % cm

%% -------------------- Parameter definitions -----------------
% Panel (a): degradation, shock, and PM parameters.
modelFields = { ...
    'alpha', ...
    'beta', ...
    'lambda', ...
    'alphaS', ...
    'betaS', ...
    'p', ...
    'q'};

modelLabels = { ...
    '$\alpha$', ...
    '$\beta$', ...
    '$\lambda$', ...
    '$\alpha_S$', ...
    '$\beta_S$', ...
    '$p$', ...
    '$q$'};

% Panel (b): cost parameters.
costFields = { ...
    'Ci', ...
    'Cm', ...
    'Cp', ...
    'Cc', ...
    'Cd'};

costLabels = { ...
    '$C_i$', ...
    '$C_m$', ...
    '$C_p$', ...
    '$C_c$', ...
    '$C_d$'};

allFields = [modelFields, costFields];
allLabels = [modelLabels, costLabels];

%% --------------------- Read baseline cost --------------------
if isfinite(baselineCostManual)
    baselineCost = baselineCostManual;
else
    baselineCandidates = { ...
        fullfile(resultsDir,baselineFileName), ...
        fullfile(pwd,baselineFileName), ...
        fullfile(fileparts(resultsDir),baselineFileName)};

    baselineFile = '';

    for k = 1:numel(baselineCandidates)
        if isfile(baselineCandidates{k})
            baselineFile = baselineCandidates{k};
            break;
        end
    end

    if isempty(baselineFile)
        error([ ...
            'Cannot find "%s". Either place it in the results folder ' ...
            'or set baselineCostManual at the top of this script.'], ...
            baselineFileName);
    end

    S0 = load(baselineFile,'results');

    if ~isfield(S0,'results')
        error('The baseline MAT file does not contain variable "results".');
    end

    if isfield(S0.results,'finalBest') && ...
            isfield(S0.results.finalBest,'costRate')

        baselineCost = S0.results.finalBest.costRate;

    elseif isfield(S0.results,'bestPolicy') && ...
            isfield(S0.results.bestPolicy,'CostRate')

        baselineCost = S0.results.bestPolicy.CostRate;

    else
        error('Cannot find the optimized baseline cost rate.');
    end
end

fprintf('\nBaseline optimized long-run cost rate = %.8f\n',baselineCost);

%% ---------------- Read all sensitivity results ---------------
nPar = numel(allFields);

lowCost   = nan(nPar,1);
highCost  = nan(nPar,1);
lowRI     = nan(nPar,1);
highRI    = nan(nPar,1);
lowMOverL = nan(nPar,1);
highMOverL= nan(nPar,1);

for i = 1:nPar

    fieldName = allFields{i};

    resultFile = fullfile( ...
        resultsDir, ...
        sprintf('sensitivity_%s_50pct_results.mat',fieldName));

    if ~isfile(resultFile)
        error('Missing sensitivity result file:\n%s',resultFile);
    end

    S = load(resultFile);

    if isfield(S,'sensitivityTable')
        T = S.sensitivityTable;
    elseif isfield(S,'results') && isfield(S.results,'table')
        T = S.results.table;
    else
        error('Cannot find sensitivity table in:\n%s',resultFile);
    end

    idxLow  = find(abs(T.Multiplier-0.50) < 1e-10,1,'first');
    idxHigh = find(abs(T.Multiplier-1.50) < 1e-10,1,'first');

    if isempty(idxLow) || isempty(idxHigh)
        error('%s does not contain both -50%% and +50%% cases.',resultFile);
    end

    lowCost(i)  = T.OptimalCostRate(idxLow);
    highCost(i) = T.OptimalCostRate(idxHigh);

    if ismember('OptimalRI',T.Properties.VariableNames)
        lowRI(i)  = T.OptimalRI(idxLow);
        highRI(i) = T.OptimalRI(idxHigh);
    end

    if ismember('OptimalMOverL',T.Properties.VariableNames)
        lowMOverL(i)  = T.OptimalMOverL(idxLow);
        highMOverL(i) = T.OptimalMOverL(idxHigh);
    end
end

%% ------------------- Sensitivity magnitude -------------------
% Distance from the baseline cost. This is used only for ordering rows.
impact = max( ...
    abs([lowCost-baselineCost, highCost-baselineCost]), ...
    [],2);

modelIdx = 1:numel(modelFields);
costIdx  = numel(modelFields) + (1:numel(costFields));

if sortByImpact
    [~,orderModel] = sort(impact(modelIdx),'descend');
    [~,orderCost]  = sort(impact(costIdx),'descend');

    modelIdx = modelIdx(orderModel);
    costIdx  = costIdx(orderCost);
end

%% ----------------------- Common x-axis -----------------------
allDisplayedCosts = [lowCost; highCost; baselineCost];

xmin = min(allDisplayedCosts);
xmax = max(allDisplayedCosts);

xspan = xmax-xmin;

if xspan <= 0
    xspan = 1;
end

xpad = 0.08*xspan;

% Panel (a): data-driven x-axis.
xLimitsModel = [xmin-xpad, xmax+xpad];

% Panel (b): keep the same lower-side logic, but fix the maximum at 50.
costDisplayed = [lowCost(costIdx); highCost(costIdx); baselineCost];
xminCost = min(costDisplayed);
xspanCost = max(50-xminCost,1);
xpadCost = 0.08*xspanCost;
xLimitsCost = [xminCost-xpadCost, 45];

%% ------------------------- Plot figure -----------------------
fig = figure( ...
    'Name','Sensitivity analysis of long-run cost rate', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',figurePosition);

tl = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Colors.
lowColor  = [0.0000 0.4470 0.7410];
highColor = [0.8500 0.3250 0.0980];

% Panel (a).
ax1 = nexttile(tl,1);

[hLow1,hHigh1] = drawSensitivityPanel( ...
    ax1,modelIdx,allLabels,lowCost,highCost,baselineCost, ...
    xLimitsModel,lowColor,highColor);

xlabel(ax1,'(a)','Interpreter','latex');

text(ax1,0.5,-0.14,'(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'Interpreter','latex');

% Panel (b).
ax2 = nexttile(tl,2);

drawSensitivityPanel( ...
    ax2,costIdx,allLabels,lowCost,highCost,baselineCost, ...
    xLimitsCost,lowColor,highColor);

xlabel(ax2,'(b)','Interpreter','latex');

text(ax2,0.5,-0.14,'(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'Interpreter','latex');

% Shared legend.
lgd = legend( ...
    ax1,[hLow1,hHigh1], ...
    {'Parameter $-50\%$ (Low input)', ...
     'Parameter $+50\%$ (High input)'}, ...
    'Orientation','horizontal', ...
    'Interpreter','latex', ...
    'Location','northoutside');

lgd.Layout.Tile = 'north';

%% ---------------------- Print summary ------------------------
fprintf('\n============================================================\n');
fprintf(' Sensitivity summary\n');
fprintf('============================================================\n');

fprintf('%-10s  %-14s  %-14s  %-14s\n', ...
    'Parameter','Low (-50%)','Baseline','High (+50%)');

for i = 1:nPar
    fprintf('%-10s  %-14.6f  %-14.6f  %-14.6f\n', ...
        allFields{i},lowCost(i),baselineCost,highCost(i));
end

%% ----------------------- Save summary ------------------------
Parameter = string(allFields(:));
Label = string(allLabels(:));
Low50CostRate = lowCost;
BaselineCostRate = repmat(baselineCost,nPar,1);
High50CostRate = highCost;
SensitivityMagnitude = impact;
Low50OptimalRI = lowRI;
High50OptimalRI = highRI;
Low50OptimalMOverL = lowMOverL;
High50OptimalMOverL = highMOverL;

summaryTable = table( ...
    Parameter,Label, ...
    Low50CostRate,BaselineCostRate,High50CostRate, ...
    SensitivityMagnitude, ...
    Low50OptimalRI,High50OptimalRI, ...
    Low50OptimalMOverL,High50OptimalMOverL);

writetable( ...
    summaryTable, ...
    fullfile(resultsDir,'sensitivity_50pct_summary.csv'));

save( ...
    fullfile(resultsDir,'sensitivity_50pct_summary.mat'), ...
    'summaryTable','baselineCost');

%% ----------------------- Save figure -------------------------
savefig( ...
    fig, ...
    fullfile(resultsDir,'sensitivity_50pct_tornado_v2.fig'));

exportgraphics( ...
    fig, ...
    fullfile(resultsDir,'sensitivity_50pct_tornado_v2.png'), ...
    'Resolution',600);

exportgraphics( ...
    fig, ...
    fullfile(resultsDir,'sensitivity_50pct_tornado_v2.pdf'), ...
    'ContentType','vector');

fprintf('\nSaved:\n');
fprintf('  sensitivity_50pct_tornado_v2.fig\n');
fprintf('  sensitivity_50pct_tornado_v2.png\n');
fprintf('  sensitivity_50pct_tornado_v2.pdf\n');
fprintf('  sensitivity_50pct_summary.csv\n');
fprintf('  sensitivity_50pct_summary.mat\n');

%% ========================================================================
%  Local function
% ========================================================================
function [hLow,hHigh] = drawSensitivityPanel( ...
    ax,idx,labels,lowCost,highCost,baselineCost, ...
    xLimits,lowColor,highColor)

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

nRow = numel(idx);
y = 1:nRow;

% Draw the -50% and +50% bars on the same row for each parameter,
% as in a standard tornado plot. Both bars start from the baseline cost.
barHeight = 0.56;

hLow = barh( ...
    ax,y,lowCost(idx),barHeight, ...
    'BaseValue',baselineCost, ...
    'FaceColor',lowColor, ...
    'EdgeColor','none');

hHigh = barh( ...
    ax,y,highCost(idx),barHeight, ...
    'BaseValue',baselineCost, ...
    'FaceColor',highColor, ...
    'EdgeColor','none');

xline( ...
    ax,baselineCost,'--k', ...
    'LineWidth',1.1, ...
    'HandleVisibility','off');

ax.YDir = 'reverse';
ax.YTick = y;
ax.YTickLabel = labels(idx);
ax.TickLabelInterpreter = 'latex';

ax.XLim = xLimits;
ax.YLim = [0.45,nRow+0.55];

ax.FontName = 'Times New Roman';
ax.FontSize = 10;
ax.LineWidth = 0.8;

end
