function run_all_sensitivity_and_plot
%RUN_ALL_SENSITIVITY_AND_PLOT
% Run all twelve +/-50% sensitivity analyses, then create a two-panel
% tornado-style figure of the optimized long-run cost rate.
%
% The baseline cost is read from:
%   shock_cbm_results_v5_coarse_timed.mat
%
% Put that baseline MAT file either in this folder, the current folder,
% or the parent folder before running this function.

close all;
clc;

thisDir = fileparts(mfilename('fullpath'));
if isempty(thisDir)
    thisDir = pwd;
end

scripts = { ...
    'sensitivity_Ci_50pct.m', ...
    'sensitivity_Cm_50pct.m', ...
    'sensitivity_Cp_50pct.m', ...
    'sensitivity_Cc_50pct.m', ...
    'sensitivity_Cd_50pct.m', ...
    'sensitivity_alpha_50pct.m', ...
    'sensitivity_beta_50pct.m', ...
    'sensitivity_lambda_50pct.m', ...
    'sensitivity_alphaS_50pct.m', ...
    'sensitivity_betaS_50pct.m', ...
    'sensitivity_p_50pct.m', ...
    'sensitivity_q_50pct.m'};

fprintf('\n============================================================\n');
fprintf(' Running all 12 sensitivity analyses\n');
fprintf('============================================================\n');

for k = 1:numel(scripts)
    fprintf('\n\n########## %d/%d: %s ##########\n', ...
        k,numel(scripts),scripts{k});

    scriptPath = fullfile(thisDir,scripts{k});
    if ~isfile(scriptPath)
        error('Missing sensitivity script: %s',scriptPath);
    end

    % Run in the base workspace so that clearvars inside each script
    % does not remove variables local to this driver function.
    evalin('base',sprintf('run(''%s'');',strrep(scriptPath,,'')));
end

% fprintf('
% ============================================================
% ');
% fprintf(' All sensitivity runs completed. Creating figure...
% ');
% fprintf('============================================================
% ');

makeSensitivityTornadoFigure(thisDir);

end

%% ========================================================================
function makeSensitivityTornadoFigure(thisDir)

baselineFileName = 'shock_cbm_results_v5_coarse_timed.mat';

candidates = { ...
    fullfile(thisDir,baselineFileName), ...
    fullfile(pwd,baselineFileName), ...
    fullfile(fileparts(thisDir),baselineFileName)};

baselineFile = '';
for k = 1:numel(candidates)
    if isfile(candidates{k})
        baselineFile = candidates{k};
        break;
    end
end

if isempty(baselineFile)
    error(['Baseline result file "%s" was not found. ' ...
           'Copy it into the sensitivity folder, current folder, or parent folder.'], ...
           baselineFileName);
end

S0 = load(baselineFile,'results');
if ~isfield(S0,'results')
    error('The baseline MAT file does not contain variable "results".');
end

if isfield(S0.results,'finalBest') && isfield(S0.results.finalBest,'costRate')
    baselineCost = S0.results.finalBest.costRate;
elseif isfield(S0.results,'bestPolicy') && isfield(S0.results.bestPolicy,'CostRate')
    baselineCost = S0.results.bestPolicy.CostRate;
else
    error('Cannot find the baseline optimized cost rate in the baseline MAT file.');
end

defs = struct( ...
    'field',   {'lambda','beta','alpha','betaS','p','alphaS','q', ...
                'Ci','Cm','Cc','Cp','Cd'}, ...
    'label',   {'\lambda','eta','lpha','eta_S','p','lpha_S','q', ...
                'C_i','C_m','C_c','C_p','C_d'}, ...
    'group',   {'model','model','model','model','model','model','model', ...
                'cost','cost','cost','cost','cost'});

n = numel(defs);
lowCost  = nan(n,1);
highCost = nan(n,1);

for i = 1:n
    f = defs(i).field;
    resultFile = fullfile(thisDir,sprintf('sensitivity_%s_50pct_results.mat',f));

    if ~isfile(resultFile)
        error('Missing result file: %s',resultFile);
    end

    S = load(resultFile,'sensitivityTable');
    T = S.sensitivityTable;

    idxLow  = find(abs(T.Multiplier-0.5)<1e-12,1);
    idxHigh = find(abs(T.Multiplier-1.5)<1e-12,1);

    if isempty(idxLow) || isempty(idxHigh)
        error('Result file %s does not contain both -50%% and +50%% cases.',resultFile);
    end

    lowCost(i)  = T.OptimalCostRate(idxLow);
    highCost(i) = T.OptimalCostRate(idxHigh);
end

% Sensitivity magnitude used only to order rows within each panel.
impact = max(abs([lowCost-baselineCost, highCost-baselineCost]),[],2);

modelIdx = find(strcmp({defs.group},'model'));
costIdx  = find(strcmp({defs.group},'cost'));

[~,ordModel] = sort(impact(modelIdx),'descend');
[~,ordCost]  = sort(impact(costIdx),'descend');

modelIdx = modelIdx(ordModel);
costIdx  = costIdx(ordCost);

allCosts = [lowCost; highCost; baselineCost];
xmin = min(allCosts);
xmax = max(allCosts);
pad = 0.10*max(xmax-xmin,1);
xLimits = [xmin-pad, xmax+pad];

fig = figure( ...
    'Name','Sensitivity of optimized long-run cost rate', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 22 10]);

tl = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

ax1 = nexttile(tl,1);
drawPanel(ax1,modelIdx,defs,lowCost,highCost,baselineCost,xLimits);
xlabel(ax1,'Long-run Cost Rate','Interpreter','latex');
text(ax1,0.5,-0.14,'(a)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'Interpreter','latex');

ax2 = nexttile(tl,2);
drawPanel(ax2,costIdx,defs,lowCost,highCost,baselineCost,xLimits);
xlabel(ax2,'Long-run Cost Rate','Interpreter','latex');
text(ax2,0.5,-0.14,'(b)', ...
    'Units','normalized', ...
    'HorizontalAlignment','center', ...
    'FontWeight','bold', ...
    'Interpreter','latex');

% Shared legend above the two panels.
hLow  = plot(ax1,nan,nan,'s', ...
    'MarkerFaceColor',[0 0.4470 0.7410], ...
    'MarkerEdgeColor',[0 0.4470 0.7410], ...
    'MarkerSize',8);
hHigh = plot(ax1,nan,nan,'s', ...
    'MarkerFaceColor',[0.8500 0.3250 0.0980], ...
    'MarkerEdgeColor',[0.8500 0.3250 0.0980], ...
    'MarkerSize',8);

lgd = legend(ax1,[hLow,hHigh], ...
    {'Parameter -50\% (Low input)','Parameter +50\% (High input)'}, ...
    'Orientation','horizontal', ...
    'Interpreter','latex', ...
    'Location','northoutside');
lgd.Layout.Tile = 'north';

% Save summary data.
Parameter = string({defs.field})';
Label = string({defs.label})';
Low50CostRate = lowCost;
BaselineCostRate = repmat(baselineCost,n,1);
High50CostRate = highCost;
SensitivityMagnitude = impact;

summaryTable = table( ...
    Parameter,Label,Low50CostRate,BaselineCostRate, ...
    High50CostRate,SensitivityMagnitude);

writetable(summaryTable,fullfile(thisDir,'sensitivity_50pct_summary.csv'));
save(fullfile(thisDir,'sensitivity_50pct_summary.mat'), ...
    'summaryTable','baselineCost');

exportgraphics(fig, ...
    fullfile(thisDir,'sensitivity_50pct_tornado.png'), ...
    'Resolution',600);

exportgraphics(fig, ...
    fullfile(thisDir,'sensitivity_50pct_tornado.pdf'), ...
    'ContentType','vector');

savefig(fig,fullfile(thisDir,'sensitivity_50pct_tornado.fig'));

fprintf('
Baseline optimized cost rate = %.8f
',baselineCost);
fprintf('Saved sensitivity figure and summary files in:
  %s
',thisDir);

end

%% ========================================================================
function drawPanel(ax,idx,defs,lowCost,highCost,baselineCost,xLimits)

hold(ax,'on');
box(ax,'on');
grid(ax,'on');

nRow = numel(idx);
y = 1:nRow;

blue   = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];

% Use two thin horizontal bars per parameter. This remains readable even
% when the low and high cases move the cost in the same direction.
barOffset = 0.17;
barWidth  = 0.30;

barh(ax,y-barOffset,lowCost(idx),barWidth, ...
    'BaseValue',baselineCost, ...
    'FaceColor',blue, ...
    'EdgeColor','none');

barh(ax,y+barOffset,highCost(idx),barWidth, ...
    'BaseValue',baselineCost, ...
    'FaceColor',orange, ...
    'EdgeColor','none');

xline(ax,baselineCost,'--k','LineWidth',1.1, ...
    'HandleVisibility','off');

ax.YDir = 'reverse';
ax.YTick = y;
ax.YTickLabel = {defs(idx).label};
ax.TickLabelInterpreter = 'latex';
ax.FontName = 'Times New Roman';
ax.FontSize = 10;
ax.LineWidth = 0.8;
ax.XLim = xLimits;
ax.YLim = [0.45,nRow+0.55];

end
