function T = export_section65_grid_reliability_values()
% EXPORT_SECTION65_GRID_RELIABILITY_VALUES
% Report actual R_m(t|0) and R_m(t) from Coarse, Baseline, and Fine grids.
% Reuses the SAVED baseline analytical calculation from section65_all_results.mat.
% Only the Coarse/Fine analytical calculations are rerun; NO Monte Carlo.
%
% Put this file, section65_analytical_reference.m, and
% section65_all_results.mat in the same directory, or keep the MAT file
% in the 'section65_results' subfolder.
%
% Run:
%   T = export_section65_grid_reliability_values();
%
% Outputs in section65_results:
%   Table_section65_reliability_values.csv
%   Data_section65_grid_curves.mat
%
root = fileparts(mfilename('fullpath'));
matFile = fullfile(root, 'section65_results', 'section65_all_results.mat');
if ~isfile(matFile)
    matFile = fullfile(root, 'section65_all_results.mat');
end
if ~isfile(matFile)
    error('Cannot find section65_all_results.mat. Place it in this folder or section65_results.');
end
L = load(matFile, 'results');
res = L.results;
par = res.parameters;
policy = res.policy;
opts = res.analyticalOptions;
baseline = res.analytical;

coarse = opts;
coarse.Nz = 30; coarse.Nv = 40; coarse.Ny = 10;
coarse.NtInt = 40; coarse.Nconv = 30;
coarse.NuRel = 10; coarse.NtRel = 40;
coarse.warnOnMass = false;

fine = opts;
fine.Nz = 90; fine.Nv = 120; fine.Ny = 28;
fine.NtInt = 120; fine.Nconv = 90;
fine.NuRel = 28; fine.NtRel = 90;
fine.warnOnMass = false;

fprintf('Calculating Coarse grid reliability ...\n');
[repCoarse, ~] = section65_analytical_reference( ...
    policy.RI, policy.M, par, coarse, [], true);
fprintf('Reusing saved Baseline grid reliability ...\n');
repBaseline = baseline;
fprintf('Calculating Fine grid reliability ...\n');
[repFine, ~] = section65_analytical_reference( ...
    policy.RI, policy.M, par, fine, [], true);

curves = struct();
curves.Coarse = repCoarse.reliability;
curves.Baseline = repBaseline.reliability;
curves.Fine = repFine.reliability;

% Evaluate ALL three curves at identical physical times, not at
% grid indices (different discretizations use different time nodes).
tmax = min([curves.Coarse.t(end), curves.Baseline.t(end), curves.Fine.t(end)]);
tValues = [0; 5; 10; 20; 30; 40; 50];
tValues = tValues(tValues <= tmax + 1e-10);
assert(~isempty(tValues), 'No common time points in the computed horizon.');

names = {'Coarse', 'Baseline', 'Fine'};
newValues = zeros(numel(tValues), 3);
statValues = zeros(numel(tValues), 3);
for k=1:3
    r = curves.(names{k});
    newValues(:,k) = interp1(r.t(:), r.RfromNew(:), tValues, 'linear');
    statValues(:,k) = interp1(r.t(:), r.Rstationary(:), tValues, 'linear');
end

T = table(tValues, newValues(:,1), newValues(:,2), newValues(:,3), ...
    statValues(:,1), statValues(:,2), statValues(:,3), ...
    'VariableNames', {'Time', 'RnewCoarse', 'RnewBaseline', 'RnewFine', ...
    'RstatCoarse', 'RstatBaseline', 'RstatFine'});

outDir = fullfile(root, 'section65_results');
if ~exist(outDir, 'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir, 'Table_section65_reliability_values.csv'));
save(fullfile(outDir, 'Data_section65_grid_curves.mat'), 'curves', '-v7.3');
disp(T);
fprintf('Saved all three grid values to: %s\n', outDir);
end
