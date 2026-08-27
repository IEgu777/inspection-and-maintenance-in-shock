%% redraw_coarse_cost_surface_projected.m
% Draw a three-dimensional long-run cost-rate surface with
% bottom and top heat-map projections.
%
% This script does NOT rerun the optimization.
%
% Required data:
%   results.coarseGridTable
%   results.parameters.L
%
% The coarse-grid table must contain:
%   RI, M, CostRate

fprintf('\nRedrawing projected three-dimensional cost-rate surface...\n');

%% ------------------------------------------------------------------------
%  1. Get data from workspace or saved MAT file
% -------------------------------------------------------------------------
if exist('results','var') && isstruct(results) && ...
        isfield(results,'coarseGridTable') && ...
        isfield(results,'parameters')

    tabLocal = results.coarseGridTable;
    parLocal = results.parameters;

elseif exist('coarseTab','var') && exist('par','var')

    tabLocal = coarseTab;
    parLocal = par;

elseif exist('ans','var') && isstruct(ans) && ...
        isfield(ans,'coarseGridTable') && ...
        isfield(ans,'parameters')

    results = ans; %#ok<NASGU>

    tabLocal = ans.coarseGridTable;
    parLocal = ans.parameters;

else

    matName = 'shock_cbm_results_v5_coarse_timed.mat';

    if exist(matName,'file') ~= 2
        error(['Cannot find workspace variables or ', ...
               matName, ...
               ' in the current MATLAB folder.']);
    end

    S = load(matName,'results');

    if ~isfield(S,'results')
        error('%s does not contain the variable results.',matName);
    end

    results = S.results; %#ok<NASGU>

    if ~isfield(results,'coarseGridTable') || ...
            ~isfield(results,'parameters')

        error(['The saved results structure must contain ', ...
               'coarseGridTable and parameters.']);
    end

    tabLocal = results.coarseGridTable;
    parLocal = results.parameters;
end

%% ------------------------------------------------------------------------
%  2. Validate input data
% -------------------------------------------------------------------------
if ~istable(tabLocal)
    error('The coarse-grid result must be a MATLAB table.');
end

requiredColumns = {'RI','M','CostRate'};

if ~all(ismember(requiredColumns, ...
        tabLocal.Properties.VariableNames))

    error('The coarse-grid table must contain RI, M, and CostRate.');
end

if ~isstruct(parLocal) || ...
        ~isfield(parLocal,'L') || ...
        ~isscalar(parLocal.L) || ...
        ~isfinite(parLocal.L) || ...
        parLocal.L <= 0

    error('The parameter structure must contain a positive scalar par.L.');
end

valid = ...
    isfinite(tabLocal.RI) & ...
    isfinite(tabLocal.M) & ...
    isfinite(tabLocal.CostRate);

tabLocal = tabLocal(valid,:);

if isempty(tabLocal)
    error('No finite coarse-grid cost values are available.');
end

%% ------------------------------------------------------------------------
%  3. Reconstruct the rectangular coarse grid
% -------------------------------------------------------------------------
RIvals = unique(tabLocal.RI,'sorted');
Mvals  = unique(tabLocal.M,'sorted');

% Normalize M using the failure threshold L.
MLvals = Mvals ./ parLocal.L;

Z = nan(numel(RIvals),numel(Mvals));

for i = 1:numel(RIvals)

    for j = 1:numel(Mvals)

        idx = ...
            abs(tabLocal.RI - RIvals(i)) < 1e-12 & ...
            abs(tabLocal.M  - Mvals(j))  < 1e-12;

        if any(idx)

            Z(i,j) = ...
                tabLocal.CostRate(find(idx,1,'first'));
        end
    end
end

if all(~isfinite(Z(:)))
    error('The reconstructed cost surface contains no finite values.');
end

%% ------------------------------------------------------------------------
%  4. User-adjustable plotting settings
% -------------------------------------------------------------------------

% x-axis range: M/L
% Use [] to show the complete optimization range.
xShow = [0.0 1.0];

% y-axis range: reliability threshold R_I
% Use [] to show the complete optimization range.
yShow = [0.0 1.0];

% z-axis range: long-run cost rate
%
% Example:
% zShow = [30 60];
%
% Use [] for automatic limits.
zShow = [];

% Whether to interpolate the coarse-grid data.
useSmoothInterpolation = true;

% Number of interpolation points.
nFineX = 180;
nFineY = 180;

% Number of discrete color bands.
nColorBands = 18;

% Number of contour levels on the top and bottom planes.
nContourLines = 14;

showCoarseMesh       = false;
showBottomProjection = true;
showTopProjection    = true;

surfaceAlpha = 1.00;
bottomAlpha  = 0.90;
topAlpha     = 0.15;

azimuth   = -50;
elevation = 18;

% Coordinate tick spacing.
% Set to [] to use MATLAB automatic ticks.
xTickStep = 0.05;
yTickStep = 0.05;
zTickStep = [];

%% ------------------------------------------------------------------------
%  5. Restrict data to selected x-y range
% -------------------------------------------------------------------------
ix = true(size(MLvals));
iy = true(size(RIvals));

if ~isempty(xShow)

    validateRange(xShow,'xShow');

    ix = ...
        MLvals >= xShow(1) & ...
        MLvals <= xShow(2);
end

if ~isempty(yShow)

    validateRange(yShow,'yShow');

    iy = ...
        RIvals >= yShow(1) & ...
        RIvals <= yShow(2);
end

xBase = MLvals(ix);
yBase = RIvals(iy);
ZBase = Z(iy,ix);

if numel(xBase) < 2 || numel(yBase) < 2

    error(['The selected x-y range contains fewer than two ', ...
           'coarse-grid points in at least one direction.']);
end

if all(~isfinite(ZBase(:)))

    error('No finite cost values exist in the selected x-y range.');
end

%% ------------------------------------------------------------------------
%  6. Interpolate the surface
% -------------------------------------------------------------------------
[Xbase,Ybase] = meshgrid(xBase,yBase);

if useSmoothInterpolation

    xPlot = linspace( ...
        min(xBase), ...
        max(xBase), ...
        nFineX);

    yPlot = linspace( ...
        min(yBase), ...
        max(yBase), ...
        nFineY);

    [Xplot,Yplot] = meshgrid(xPlot,yPlot);

    finiteBase = isfinite(ZBase);

    if all(finiteBase(:))

        % Use makima when supported.
        try

            Zplot = interp2( ...
                Xbase, ...
                Ybase, ...
                ZBase, ...
                Xplot, ...
                Yplot, ...
                'makima');

        catch

            % Fallback for older MATLAB versions.
            Zplot = interp2( ...
                Xbase, ...
                Ybase, ...
                ZBase, ...
                Xplot, ...
                Yplot, ...
                'spline');
        end

    else

        % More robust when a few grid values are missing.
        F = scatteredInterpolant( ...
            Xbase(finiteBase), ...
            Ybase(finiteBase), ...
            ZBase(finiteBase), ...
            'natural', ...
            'none');

        Zplot = F(Xplot,Yplot);
    end

else

    xPlot = xBase;
    yPlot = yBase;

    Xplot = Xbase;
    Yplot = Ybase;
    Zplot = ZBase;
end

finiteZ = Zplot(isfinite(Zplot));

if isempty(finiteZ)
    error('The plotted surface contains no finite cost values.');
end

zDataMin = min(finiteZ);
zDataMax = max(finiteZ);

zDataRange = max(zDataMax-zDataMin,eps);

%% ------------------------------------------------------------------------
%  7. Determine z-axis and projection-plane positions
% -------------------------------------------------------------------------
if isempty(zShow)

    % Automatically leave space below and above the real cost surface.
    zLower = zDataMin - 0.12*zDataRange;
    zUpper = zDataMax + 0.12*zDataRange;

else

    validateRange(zShow,'zShow');

    zLower = zShow(1);
    zUpper = zShow(2);
end

% Projection-plane locations.
zFloor   = zLower;
zCeiling = zUpper;

%% ------------------------------------------------------------------------
%  8. Convert continuous costs into discrete color bands
% -------------------------------------------------------------------------
if nColorBands < 2 || floor(nColorBands) ~= nColorBands
    error('nColorBands must be an integer greater than or equal to 2.');
end

colorEdges = linspace( ...
    zDataMin, ...
    zDataMax, ...
    nColorBands+1);

colorCenters = ...
    0.5*(colorEdges(1:end-1)+colorEdges(2:end));

Cband = nan(size(Zplot));

validZ = isfinite(Zplot);

bandIndex = floor( ...
    (Zplot(validZ)-zDataMin) ./ ...
    zDataRange * nColorBands) + 1;

bandIndex = min(max(bandIndex,1),nColorBands);

Cband(validZ) = colorCenters(bandIndex);

contourLevels = linspace( ...
    zDataMin, ...
    zDataMax, ...
    nContourLines);

%% ------------------------------------------------------------------------
%  9. Create figure
% -------------------------------------------------------------------------
fig = figure( ...
    'Name','Projected three-dimensional cost-rate surface', ...
    'Color','w', ...
    'Position',[80 55 1020 740]);

ax = axes('Parent',fig);
hold(ax,'on');

%% ------------------------------------------------------------------------
%  10. Bottom heat-map projection
% -------------------------------------------------------------------------
if showBottomProjection

    surf(ax, ...
        Xplot, ...
        Yplot, ...
        zFloor*ones(size(Zplot)), ...
        Cband, ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',bottomAlpha);

    % Add projected contour lines.
    plotProjectedContours( ...
        ax, ...
        xPlot, ...
        yPlot, ...
        Zplot, ...
        contourLevels, ...
        zFloor, ...
        [0.12 0.12 0.12], ...
        0.55);
end

%% ------------------------------------------------------------------------
%  11. Top heat-map projection
% -------------------------------------------------------------------------
if showTopProjection

    surf(ax, ...
        Xplot, ...
        Yplot, ...
        zCeiling*ones(size(Zplot)), ...
        Cband, ...
        'FaceColor','flat', ...
        'EdgeColor','none', ...
        'FaceAlpha',topAlpha);

    % Add projected contour lines.
    plotProjectedContours( ...
        ax, ...
        xPlot, ...
        yPlot, ...
        Zplot, ...
        contourLevels, ...
        zCeiling, ...
        [0.12 0.12 0.12], ...
        0.50);
end

%% ------------------------------------------------------------------------
%  12. Central three-dimensional cost surface
% -------------------------------------------------------------------------
hMainSurface = surf(ax, ...
    Xplot, ...
    Yplot, ...
    Zplot, ...
    Zplot, ...
    'FaceColor','interp', ...
    'EdgeColor',[0.15 0.15 0.15], ...
    'EdgeAlpha',0.28, ...
    'LineWidth',0.25, ...
    'FaceAlpha',surfaceAlpha);

%% ------------------------------------------------------------------------
%  13. Overlay original coarse-grid mesh
% -------------------------------------------------------------------------
if showCoarseMesh

    mesh(ax, ...
        Xbase, ...
        Ybase, ...
        ZBase, ...
        'FaceColor','none', ...
        'EdgeColor','k', ...
        'EdgeAlpha',0.62, ...
        'LineWidth',0.50);
end

%% ------------------------------------------------------------------------
%  14. Axis labels and title
% -------------------------------------------------------------------------
xlabel(ax, ...
    'Preventive maintenance threshold M/L');

ylabel(ax, ...
    'Reliability threshold R_I');

zlabel(ax, ...
    'Long-run cost rate');

title(ax, ...
    'Long-run cost-rate surface');

%% ------------------------------------------------------------------------
%  15. Colormap and colorbar
% -------------------------------------------------------------------------
try

    % Low values are blue and high values are red.
    colormap(ax,turbo(nColorBands));

catch

    % Fallback for older MATLAB versions.
    colormap(ax,jet(nColorBands));
end

caxis(ax,[zDataMin zDataMax]);

cb = colorbar(ax);
ylabel(cb,'Long-run cost rate');

cb.Ticks = linspace( ...
    zDataMin, ...
    zDataMax, ...
    min(9,nColorBands+1));

%% ------------------------------------------------------------------------
%  16. Axis ranges
% -------------------------------------------------------------------------
if ~isempty(xShow)

    xlim(ax,xShow);

else

    xlim(ax,[min(xBase) max(xBase)]);
end

if ~isempty(yShow)

    ylim(ax,yShow);

else

    ylim(ax,[min(yBase) max(yBase)]);
end

zlim(ax,[zLower zUpper]);

%% ------------------------------------------------------------------------
%  17. Axis ticks
% -------------------------------------------------------------------------
if ~isempty(xTickStep)

    xLimits = xlim(ax);

    xticks(ax, ...
        xLimits(1):xTickStep:xLimits(2));
end

if ~isempty(yTickStep)

    yLimits = ylim(ax);

    yticks(ax, ...
        yLimits(1):yTickStep:yLimits(2));
end

if ~isempty(zTickStep)

    zticks(ax, ...
        zLower:zTickStep:zUpper);
end

%% ------------------------------------------------------------------------
%  18. Appearance and view
% -------------------------------------------------------------------------
grid(ax,'on');
box(ax,'on');

set(ax, ...
    'Layer','top', ...
    'FontSize',10, ...
    'Projection','perspective');

view(ax,azimuth,elevation);

axis(ax,'vis3d');

% Allow interactive mouse rotation.
rotate3d(fig,'on');

%% ------------------------------------------------------------------------
%  19. Save figure
% -------------------------------------------------------------------------
outFig = 'coarse_cost_surface_projected.fig';
outPng = 'coarse_cost_surface_projected.png';

savefig(fig,outFig);

try

    exportgraphics( ...
        fig, ...
        outPng, ...
        'Resolution',300);

catch

    print( ...
        fig, ...
        outPng, ...
        '-dpng', ...
        '-r300');
end

fprintf('Saved projected three-dimensional cost surface to:\n');
fprintf('  %s\n',outFig);
fprintf('  %s\n',outPng);


%% ========================================================================
%  Local functions
% ========================================================================

function validateRange(v,name)
% Check a two-element plotting range.

    if ~isnumeric(v) || ...
            numel(v) ~= 2 || ...
            any(~isfinite(v)) || ...
            v(1) >= v(2)

        error('%s must have the form [minimum maximum].',name);
    end
end


function plotProjectedContours( ...
    ax,xVals,yVals,Z,levels,zPlane,lineColor,lineWidth)
% Calculate contours from Z and draw them on a fixed horizontal plane.
%
% This avoids the version-dependent ZLocation property of contour objects.

    C = contourc(xVals,yVals,Z,levels);

    column = 1;

    while column < size(C,2)

        % First column of each segment contains:
        % C(1,column): contour level
        % C(2,column): number of following points

        nPoints = C(2,column);

        firstPoint = column + 1;
        lastPoint  = column + nPoints;

        if lastPoint > size(C,2)
            break;
        end

        xSegment = C(1,firstPoint:lastPoint);
        ySegment = C(2,firstPoint:lastPoint);
        zSegment = zPlane*ones(size(xSegment));

        plot3(ax, ...
            xSegment, ...
            ySegment, ...
            zSegment, ...
            '-', ...
            'Color',lineColor, ...
            'LineWidth',lineWidth);

        column = lastPoint + 1;
    end
end