function Check_RI_Grid_Before_Run()
% Quick checks of the actual source files before expensive optimization.
rootDir = fileparts(mfilename('fullpath'));
expected = '[0.01:0.01:0.98, 0.981:0.001:0.999]';
list = {'GaAs_hybrid_fast.m','GaAs_policy1_fast.m'};
for k = 1:numel(list)
    name = list{k};
    p = fullfile(rootDir,name);
    assert(isfile(p),'Missing source file: %s',p);
    txt = fileread(p);
    assert(contains(txt,['coarse.RIGrid = ',expected,';']), ...
        'RI dense grid missing from: %s',name);
    fprintf('SOURCE VERIFIED: %s | 117 RI points | max RI = 0.999\n',p);
end
T = 0.01:0.01:0.98;
R = [T,0.981:0.001:0.999];
assert(numel(R)==117 && any(abs(R-0.995)<1e-12));
% We run from this exact directory; otherwise old same-name functions can win.
previousFolder=pwd;
cleaner=onCleanup(@()cd(previousFolder)); %#ok<NASGU>
cd(rootDir);
for k = 1:numel(list)
    [~,base] = fileparts(list{k});
    actual=which(base);
    assert(strcmpi(actual,fullfile(rootDir,list{k})), ...
        'MATLAB function-path collision for %s: %s',base,actual);
end
P=GaAs_case_parameters();
actualCosts=[P.Ci,P.Cm,P.Cp,P.Cc,P.Cd];
assert(isequal(actualCosts,[3,20,50,90,40]), ...
    'Incorrect costs; need (3,20,50,90,40).');
fprintf('COST VERIFIED: (Ci,Cm,Cp,Cc,Cd)=(3,20,50,90,40)\n');
fprintf('RI DENSE GRID VERIFIED: FULL | 117 RI x 99 M = 11583 candidates\n');
end
