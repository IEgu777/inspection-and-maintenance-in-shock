function par = GaAs_case_parameters()
% Shared physical and COST parameters for all four policies.
% IMPORTANT: Shock magnitude estimates come from the NEW random-grid-shock
% EM output, never from the obsolete every-250h injection dataset.
% ALL time units: 1 model time unit = 250 hours.
%
% Illustrative cost scenario, not measured GaAs maintenance costs.
% Preventive replacement cost Cp=50 uses the literature reference scale.
% All four strategies and both shock-magnitude parameter groups use these costs.
rootDir=fileparts(mfilename('fullpath'));
summaryFile=fullfile(rootDir,'GaAs_Grid_EM_Summary.csv');
if ~isfile(summaryFile)
    error(['GaAs_Grid_EM_Summary.csv is missing. ', ...
        'Run Run_GaAs_Grid_EM or Run_GaAs_Grid_All first.']);
end
T=readtable(summaryFile);
assert(height(T)==1 && all(ismember( ...
    {'Lambda_Estimated_Frequency', ...
    'AlphaS_True','BetaS_True','AlphaS_EM','BetaS_EM'}, ...
    T.Properties.VariableNames)), 'Invalid EM summary schema.');

% Set automatically by Run_GaAs_Compare_True_Estimated.m.
% Single-policy manual execution defaults to estimated parameters.
source = getenv('GAAS_POLICY_PARAMETER_SOURCE');
if isempty(source)
    source = 'estimated';
end
par.alpha = 7.426;
par.beta = 1/0.0674;
% SHARED between TRUE shock-magnitude and ESTIMATED shock-magnitude runs.
% Calculated from observed grid shocks / observed eligible 250-h intervals.
% It matches mean event frequency to the HPP optimization model, and is
% NOT a known true continuous-time HPP intensity.
par.lambda = double(T.Lambda_Estimated_Frequency(1));
if strcmpi(source,'true')
    % Known injected shock-MAGNITUDE parameters only; lambda remains estimated.
    par.alphaS = double(T.AlphaS_True(1));
    par.betaS = double(T.BetaS_True(1));
elseif strcmpi(source,'estimated')
    % EM-estimated shock-MAGNITUDE parameters; same lambda as TRUE case.
    par.alphaS = double(T.AlphaS_EM(1));
    par.betaS = double(T.BetaS_EM(1));
else
    error('GAAS_POLICY_PARAMETER_SOURCE must be true or estimated.');
end
par.p = 6;
par.q = 5;
par.L = 10;

% Illustrative scenario, not empirically measured maintenance costs.
% Results must be recomputed after any change in cost or search grid.
% The four strategies and both magnitude-parameter groups share identical costs.
par.Ci = 3;
par.Cm = 20;
par.Cp = 50;
par.Cc = 90;
par.Cd = 40; % loss rate per 250-hour time unit

assert(par.lambda>0 && isfinite(par.lambda), ...
    'Estimated lambda must be finite and positive.');
assert(all(isfinite([par.alphaS,par.betaS])) && ...
    par.alphaS>0 && par.betaS>0,'EM estimates must be positive.');
end
