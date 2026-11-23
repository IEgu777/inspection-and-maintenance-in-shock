%% EM initialization sensitivity: alpha_S and beta_S only
% Two subplots; each subplot has its own legend.
% Read the existing results without rerunning EM or Monte Carlo.
clearvars; close all; clc;

maxIterShown = 15;  % Horizontal axis: 0-20. Change to 15 or 25 if needed.

% Find the MAT file next to this script; otherwise let the user select it.
scriptDir = fileparts(mfilename('fullpath'));
matFile = fullfile(scriptDir, 'section61_all_results.mat');
if exist(matFile, 'file') ~= 2
    [file, folder] = uigetfile('*.mat', 'Select section61_all_results.mat');
    if isequal(file, 0)
        error('No MAT file selected.');
    end
    matFile = fullfile(folder, file);
end

S = load(matFile, 'results');
if ~isfield(S, 'results') || ~isfield(S.results, 'initialFits')
    error('The MAT file does not contain results.initialFits.');
end
fits = S.results.initialFits;
if ~iscell(fits) || numel(fits) ~= 6
    error('Expected 6 EM initialization fits in results.initialFits.');
end

labels = {'Moment-based', 'Fixed_1', 'Fixed_2', 'Fixed_3', 'Fixed_4', 'Fixed_5'};
colors = lines(6);

fig = figure('Color','w', 'Units','pixels', 'Position',[100 100 820 580]);
% Explicit axes avoid problems with MATLAB subplot listeners.
ax1 = axes('Parent',fig, 'Position',[0.12 0.57 0.84 0.37]);
ax2 = axes('Parent',fig, 'Position',[0.12 0.12 0.84 0.37]);
hold(ax1,'on'); hold(ax2,'on');

hAlpha = gobjects(1,6);
hBeta  = gobjects(1,6);
for j = 1:6
    H = fits{j}.history;
    % Columns: iteration, alpha_S, beta_S, log-likelihood, parameter change.
    if size(H,2) ~= 5
        if size(H,1) == 5
            H = H.';
        else
            error('Fit %d: history is %d-by-%d; expected N-by-5.', ...
                  j, size(H,1), size(H,2));
        end
    end
    iter = H(:,1);
    keep = isfinite(iter) & iter >= 0 & iter <= maxIterShown;
    if ~any(keep)
        error('Fit %d has no history points within the selected axis range.',j);
    end
    hAlpha(j) = plot(ax1, iter(keep), H(keep,2), '-', ...
                     'LineWidth',1.5, 'Color',colors(j,:));
    hBeta(j)  = plot(ax2, iter(keep), H(keep,3), '-', ...
                     'LineWidth',1.5, 'Color',colors(j,:));
end

for ax = [ax1 ax2]
    xlim(ax,[0 maxIterShown]);
    xticks(ax,0:5:maxIterShown);
    grid(ax,'on'); box(ax,'on');
    set(ax,'FontName','Times New Roman','FontSize',11, ...
        'LineWidth',0.8, 'GridAlpha',0.18);
    xlabel(ax,'EM iteration','FontSize',12);
end
ylabel(ax1,'\alpha_S','FontSize',13);
ylabel(ax2,'\beta_S','FontSize',13);
ylim(ax1,[0 8.2]);
ylim(ax2,[0 3.1]);

% Separate six-curve legends on both subplots.
legend(ax1,hAlpha,labels,'Interpreter','none', ...
       'Location','northeast','FontSize',9,'Box','on');
legend(ax2,hBeta,labels,'Interpreter','none', ...
       'Location','northeast','FontSize',9,'Box','on');

outDir = fileparts(matFile);
pngPath = fullfile(outDir,'Figure_EM_initialization_two_panels.png');
figPath = fullfile(outDir,'Figure_EM_initialization_two_panels.fig');
if exist('exportgraphics','file') == 2
    exportgraphics(fig,pngPath,'Resolution',300);
else
    print(fig,pngPath,'-dpng','-r300');
end
savefig(fig,figPath);
fprintf('Saved: %s\nSaved: %s\n',pngPath,figPath);
fprintf('Displayed iterations 0-%d. EM results are unchanged.\n',maxIterShown);
