function results = run_section61_EM(mode, outDir, overrides)
% RUN_SECTION61_EM  论文第6.1节：冲击Gamma分布参数估计及统计性能实验。
% 用法：
%   results = run_section61_EM('pilot');  % 小规模试运行
%   results = run_section61_EM('paper');  % 正式规模（计算量较大）
%   results = run_section61_EM('paper','results61',struct('useParallel',true));
%   results = run_section61_EM('pilot','test61',struct('R',3,'B',10));
%   % 【新增】直接使用已有的 checkpoint_MC.mat，复用 Table 2 中 n=300 的 200 组数据：
%   results = run_section61_EM('tworeuse','上次实验的结果文件夹');
%   % 仅重算两阶段实验，不重复运行 Table 2 的 500 次 Bootstrap。
%   % 新结果存到结果文件夹下 two_stage_reuse_n300 子文件夹。
%   % 可在 tworeuse 模式覆盖 mList，例如：
%   results = run_section61_EM('tworeuse','结果文件夹',struct('mList',[20,50,100,300,500]));
% 不需要 Statistics and Machine Learning Toolbox；并行可选（需要 PCT）。
% Gamma(shape,rate)：randg(shape)/rate；与论文第5节的密度一致。

% 输出 CSV、MAT、PNG；Table2/3/4 分别对应三个数值实验。
totalTimer = tic;  % 记录整个程序的开始时间
if nargin<1 || isempty(mode), mode='pilot'; end
if nargin<2 || isempty(outDir), outDir=fullfile(pwd,['section61_' char(mode)]); end
if nargin<3, overrides=struct(); end
if ~exist(outDir,'dir'), mkdir(outDir); end

% 【修复新增】plotonly 模式：复用上次保存的实验结果，只重新绘图。
% 用法：run_section61_EM('plotonly','上次结果文件夹');
% 不会再次运行 Monte Carlo、Bootstrap 或两阶段估计。
if strcmpi(char(mode),'plotonly')
    matFile=fullfile(outDir,'section61_all_results.mat');
    assert(exist(matFile,'file')==2, ...
        '找不到已有结果 %s。请确认 outDir 是上次实验的输出文件夹。',matFile);
    S=load(matFile,'results');
    results=S.results;
    makeFigures(results.Table2,results.Table3,results.initialFits,results.Table4,outDir);
    fprintf('绘图结束；实验估计没有重新运行。结果文件夹：%s\n',outDir);
    return;
end

% 【新增】只运行 Table 4：利用 Table 2 的固定随机种子恢复 n=300 的原始(U,Y)，
% 并直接复用 Table 2 保存的 Known EM 估计；不再计算 Table 2/Bootstrap。
% 第一阶段依然需要单独的无冲击退化数据，它不包含在冲击观测 U=X+W 中。
if strcmpi(char(mode),'tworeuse')
    results=rerunTwoStageFromSavedMC(outDir,overrides);
    return;
end

%% 1. 模型参数、实验参数（论文表1及第6.1节）
par=struct('alpha',1.02,'beta',1.2,'lambda',0.25,'L',12.5);
truth=[2.5,1];                         % 只用于生成数据与评价误差
cfg.nList=[50,100,150,200,250,300];             % 冲击触发观测数，不是轨迹数
cfg.nInit=300;                         % 初值敏感性：同一组数据
cfg.nTwo=300;                          % 两阶段实验：复用Table 2中n=300的冲击观测
cfg.mList=[20,50,100,500];              % 第一阶段无冲击观测数
cfg.dtClean=1;                         % 第一阶段固定观测间隔（时间单位）
cfg.reuseTwoFromMC=true;              % 【新增】Table 4使用Table 2相同的冲击数据及Known估计
cfg.seed=20261116;
cfg.useParallel=false;
cfg.minBootFraction=0.95;              % 有效bootstrap比例不足则不报告该CI
cfg.makePlots=true;
cfg.runSelfTest=true;                 % 验证E-step数值积分及Gamma更新方程
% EM停止准则采用绝对参数范数，与论文公式(69)一致。
% maxIter为计算保护，不等同于达到收敛；绝不把超限结果当作收敛。
cfg.em=struct('tol',1e-6,'maxIter',2000,'maxTerms',5000, ...
              'maxShape',1e5,'monotoneTol',1e-9,'fallback',[1,1]);
switch lower(char(mode))
    case 'pilot', cfg.R=20; cfg.B=40; cfg.Rtwo=30;
    case 'paper', cfg.R=200; cfg.B=500; cfg.Rtwo=200;
    otherwise, error('mode must be pilot or paper');
end
% 允许覆盖实验参数，例如增大最大迭代次数或启用并行。
fn=fieldnames(overrides);
for j=1:numel(fn)
    if strcmp(fn{j},'em')
        sub=fieldnames(overrides.em);
        for k=1:numel(sub), cfg.em.(sub{k})=overrides.em.(sub{k}); end
    else
        if ~isfield(cfg,fn{j}), error('Unknown option: %s',fn{j}); end
        cfg.(fn{j})=overrides.(fn{j});
    end
end
assert(cfg.R>=2 && cfg.B>=2 && cfg.Rtwo>=2 && all(cfg.nList>=2));
assert(all(cfg.mList>=2) && cfg.dtClean>0 && cfg.minBootFraction>0 && cfg.minBootFraction<=1);
assert(cfg.em.maxIter>=1 && cfg.em.tol>0 && cfg.em.maxTerms>=1);
assert(cfg.nInit>=2 && cfg.nTwo>=2);
if cfg.reuseTwoFromMC && cfg.Rtwo>cfg.R
    warning('两阶段实验计划 %d 次，但Table 2仅有 %d 次；使用已有的 %d 次进行配对比较。', ...
        cfg.Rtwo,cfg.R,cfg.R);
    cfg.Rtwo=cfg.R;
end
if cfg.reuseTwoFromMC
    assert(any(cfg.nList==cfg.nTwo), ...
        'cfg.nTwo=%d 不在 cfg.nList 中。如需独立生成，设置 reuseTwoFromMC=false。',cfg.nTwo);
end
usePar=cfg.useParallel && license('test','Distrib_Computing_Toolbox');
if cfg.useParallel && ~usePar, warning('PCT unavailable; using serial execution.'); end
if usePar && isempty(gcp('nocreate')), parpool; end
fprintf('Mode=%s, R=%d, B=%d, Rtwo=%d, EM tol=%g, maxIter=%d\n', ...
    mode,cfg.R,cfg.B,cfg.Rtwo,cfg.em.tol,cfg.em.maxIter);

if cfg.runSelfTest
    check=validateEstimator(par,cfg.em);
    writetable(check,fullfile(outDir,'Numerical_selftest.csv'));
    fprintf('E-step numerical self-test passed.\n');
end

%% 2. Monte Carlo + 参数bootstrap：生成 Table 2
% 每个外层实验独立生成 n 对 (U,Y)。每次bootstrap固定该次实验的Y，
% 用该次估计参数生成新U，而非用真实参数。
% 为方便复现，每个任务使用独立、固定的随机数种子，串并行结果一致。
t2=table(); t2diag=table(); t2raw=table(); mc=cell(numel(cfg.nList),1);
for j=1:numel(cfg.nList)
    n=cfg.nList(j); fprintf('\nTable 2: n=%d\n',n); timer=tic;
    runs=cell(cfg.R,1);
    if usePar
        parfor r=1:cfg.R
            runs{r}=oneMC(n,par,truth,cfg,seedFor(cfg.seed,1,j,r));
        end
    else
        for r=1:cfg.R
            runs{r}=oneMC(n,par,truth,cfg,seedFor(cfg.seed,1,j,r));
            if mod(r,max(1,floor(cfg.R/5)))==0 || r==cfg.R
                fprintf('  %d/%d, %.1f s\n',r,cfg.R,toc(timer));
            end
        end
    end
    mc{j}=vertcat(runs{:});
    [part,diagPart]=summarizeMC(mc{j},n,truth,cfg.B);
    t2=appendTable(t2,part); t2diag=appendTable(t2diag,diagPart);
    t2raw=appendTable(t2raw,replicationDetails(mc{j},n));
    writetable(t2,fullfile(outDir,'Table2_statistical_performance.csv'));
    writetable(t2diag,fullfile(outDir,'Table2_diagnostics.csv'));
    writetable(t2raw,fullfile(outDir,'Table2_replication_details.csv'));
    save(fullfile(outDir,'checkpoint_MC.mat'),'mc','t2','t2diag','cfg','par','truth','-v7.3');
    fprintf('  n=%d finished in %.1f s\n',n,toc(timer));
end

%% 3. 初值敏感性：同一组数据，不随初值重新生成样本（Table 3）
rng(seedFor(cfg.seed,2,1,1),'twister');
[U0,Y0,path0]=simulateData(cfg.nInit,par,truth);
[initMoM,rawMom,fb]=momentInit(U0,Y0,par,cfg.em);
starts=[initMoM; .5,.5; 1,1; 2,.5; 5,2; 8,3];
labels={'MoM';'Fixed_1';'Fixed_2';'Fixed_3';'Fixed_4';'Fixed_5'};
if fb, labels{1}='MoM_fallback'; end
initFits=cell(size(starts,1),1); t3=table();
for j=1:size(starts,1)
    f=fitEM(U0,Y0,par,cfg.em,starts(j,:),true);
    initFits{j}=f;
    row=table(labels(j),starts(j,1),starts(j,2),f.theta(1),f.theta(2), ...
        f.logLik,f.iter,f.status==0,{statusText(f.status)},f.lastChange, ...
        'VariableNames',{'Initialization','InitialAlpha','InitialBeta','AlphaHat','BetaHat', ...
        'LogLikelihood','Iterations','Converged','Status','LastChange'});
    t3=appendTable(t3,row);
end
writetable(t3,fullfile(outDir,'Table3_initialization_sensitivity.csv'));
writetable(table(Y0,U0,path0,'VariableNames',{'Y','U','TrajectoryID'}), ...
    fullfile(outDir,'Table3_fixed_dataset.csv'));
% 预先指定此数据集用于第6.2节（如果决定更新优化结果，需重新求解）。
refFit=initFits{1};
writetable(table(refFit.theta(1),refFit.theta(2),refFit.status==0,{statusText(refFit.status)}, ...
    'VariableNames',{'AlphaHat','BetaHat','Converged','Status'}), ...
    fullfile(outDir,'Reference_estimate_for_section62.csv'));
% 保存每次迭代的参数与对应的观测数据对数似然，而不是上一轮似然。
histTable=table();
for j=1:numel(initFits)
    H=initFits{j}.history;
    rows=table(repmat(labels(j),size(H,1),1),H(:,1),H(:,2),H(:,3),H(:,4),H(:,5), ...
        'VariableNames',{'Initialization','Iteration','Alpha','Beta','LogLikelihood','ParameterChange'});
    histTable=appendTable(histTable,rows);
end
writetable(histTable,fullfile(outDir,'Table3_iteration_history.csv'));

%% 4. 两阶段参数估计：同一冲击数据比较 Known / Plug-in（Table 4）
% 【修改】默认复用 Table 2 中 n=cfg.nTwo 的同一批冲击数据(U,Y)。
% Table 2 仅保留了估计值、未直接保存U和Y，所以借助相同种子恢复U和Y。
% Known估计直接复用mc{jTwo}(r).theta；只有Plug-in进行新EM。
% 第一阶段另外生成独立无冲击固定间隔样本，并用Gamma MLE估计(alpha,beta)。
% 各m使用同一组无冲击数据的前m个，形成嵌套配对实验。
twoRuns=cell(cfg.Rtwo,1); fprintf('\nTable 4: two-stage experiment\n'); timer=tic;
if cfg.reuseTwoFromMC
    jTwo=find(cfg.nList==cfg.nTwo,1,'first');
    refRuns=mc{jTwo};
    assert(numel(refRuns)>=cfg.Rtwo,'Table 2 可复用的重复次数不足。');
    verifyTable2Replay(refRuns,par,truth,cfg,jTwo,cfg.Rtwo);
    fprintf('  复用 Table 2: n=%d, Rtwo=%d；Known结果直接读取，新增独立无冲击样本。\n', ...
        cfg.nTwo,cfg.Rtwo);
    if usePar
        parfor r=1:cfg.Rtwo
            twoRuns{r}=oneTwoStageReused(par,truth,cfg,jTwo,r,refRuns(r));
        end
    else
        for r=1:cfg.Rtwo
            twoRuns{r}=oneTwoStageReused(par,truth,cfg,jTwo,r,refRuns(r));
            if mod(r,max(1,floor(cfg.Rtwo/5)))==0 || r==cfg.Rtwo
                fprintf('  %d/%d, %.1f s\n',r,cfg.Rtwo,toc(timer));
            end
        end
    end
else
    % 可选：保留原先独立生成 Table 4 冲击数据的实验方式。
    if usePar
        parfor r=1:cfg.Rtwo
            twoRuns{r}=oneTwoStage(par,truth,cfg,seedFor(cfg.seed,3,1,r));
        end
    else
        for r=1:cfg.Rtwo
            twoRuns{r}=oneTwoStage(par,truth,cfg,seedFor(cfg.seed,3,1,r));
            if mod(r,max(1,floor(cfg.Rtwo/5)))==0 || r==cfg.Rtwo
                fprintf('  %d/%d, %.1f s\n',r,cfg.Rtwo,toc(timer));
            end
        end
    end
end
[t4,t4diag,t4first]=summarizeTwo(twoRuns,cfg,truth,par);
writetable(t4,fullfile(outDir,'Table4_two_stage_performance.csv'));
writetable(t4diag,fullfile(outDir,'Table4_diagnostics.csv'));
writetable(t4first,fullfile(outDir,'Table4_first_stage_estimates.csv'));

%% 5. 数值结果和图形输出
results=struct('config',cfg,'parameters',par,'truth',truth,'Table2',t2, ...
    'Table2Diagnostics',t2diag,'Table2Replications',t2raw,'Table3',t3,'Table4',t4,'Table4Diagnostics',t4diag, ...
    'Table4FirstStage',t4first,'MC',{mc},'initialFits',{initFits}, ...
    'initialData',struct('U',U0,'Y',Y0,'rawMoments',rawMom,'fallback',fb), ...
    'twoStageRuns',{twoRuns},'referenceFit',refFit);
save(fullfile(outDir,'section61_all_results.mat'),'results','-v7.3');
if cfg.makePlots, makeFigures(t2,t3,initFits,t4,outDir); end
fprintf('\n========== Table 2 ==========\n'); disp(t2);
fprintf('\n========== Table 3 ==========\n'); disp(t3);
fprintf('\n========== Table 4 ==========\n'); disp(t4);
fprintf('Results written to %s\n',outDir);
% 计算并输出整个程序的运行时间
totalTime = toc(totalTimer);

fprintf('\n========== 运行时间 ==========\n');
fprintf('总运行时间：%.2f 秒\n', totalTime);
fprintf('总运行时间：%.2f 分钟\n', totalTime/60);
fprintf('总运行时间：%.2f 小时\n', totalTime/3600);

% 将运行时间保存为CSV文件
TimeTable = table(totalTime, totalTime/60, totalTime/3600, ...
    'VariableNames', {'Seconds','Minutes','Hours'});

writetable(TimeTable, fullfile(outDir,'Computational_time.csv'));
end

%% ======================= 数据生成、初始化 =========================
function s=seedFor(base,block,j,r)
% 互不相同的实验种子；随机数调用顺序不依赖其他实验是否失败。
s=mod(double(base)+10000000*block+100000*j+r,2^32-1);
end

function [U,Y,id]=simulateData(n,par,truth)
% 保留原始代码的轨迹生成机制：只在冲击到达时观测。
% 固定n个连续收集的观测；跨越L的那个观测也纳入，然后新设备从0开始。
% 因独立增量和替换后重新开始，冲击区间增量可以一次向量化生成。
% 不模拟完整混合检测政策；这与第6.1节用于验证EM的仿真设计一致。
Y=-log(max(rand(n,1),realmin))/par.lambda;
X=randg(par.alpha.*Y)/par.beta;
W=randg(truth(1)*ones(n,1))/truth(2);
U=X+W;
if nargout>2
    id=zeros(n,1); S=0; k=1;
    for i=1:n
        id(i)=k; S=S+U(i);
        if S>=par.L, S=0; k=k+1; end
    end
end
end

function [theta,raw,fallback]=momentInit(U,Y,par,opt)
% 论文公式(58)-(59)：不能将负残差截断为正，也不能漏掉连续退化方差。
Z=U(:)-par.alpha/par.beta.*Y(:);
mu=mean(Z); v=var(Z,0)-par.alpha/par.beta^2*mean(Y);
raw=[mu,v];
fallback=~(isfinite(mu)&&isfinite(v)&&mu>0&&v>0);
if ~fallback
    theta=[mu^2/v,mu/v];
    fallback=any(~isfinite(theta)) || any(theta<=0);
end
if fallback, theta=opt.fallback; end
end

%% ============================ EM 主算法 ============================
function fit=fitEM(U,Y,par,opt,initial,keepHistory)
% status: 0=收敛, 1=迭代上限, 2=数值/输入失败, 3=似然下降,
%         4=超出数值保护范围（不是人为截断参数的估计值）。
U=U(:); Y=Y(:); n=numel(U);
if isempty(initial), [theta,~,fb]=momentInit(U,Y,par,opt);
else, theta=initial; fb=false; end
fit=struct('theta',theta,'logLik',NaN,'iter',0,'status',2, ...
    'lastChange',NaN,'fallback',fb,'history',zeros(0,5));
if n<2 || numel(Y)~=n || any(~isfinite(U)|U<=0|~isfinite(Y)|Y<=0) || ...
        any(~isfinite(theta)|theta<=0) || any(~isfinite([par.alpha,par.beta])) || ...
        par.alpha<=0 || par.beta<=0, return; end
% 缓存仅依赖于数据和已知连续退化参数的部分，避免每次重复计算。
d=struct('U',U,'logU',log(U),'aX',par.alpha.*Y,'beta',par.beta);
d.base=d.aX.*log(par.beta)+(d.aX-1).*d.logU;
[A,B,ll,code]=estep(d,theta,opt);
fit.status=code; fit.logLik=ll;
if code~=0, return; end
if keepHistory, H=NaN(opt.maxIter+1,5); H(1,:)=[0,theta,ll,NaN]; end
fit.status=1;
for r=1:opt.maxIter
    % M-step: 论文公式(67)-(68)。唯一的一维形状方程；rate解析更新。
    ma=mean(A); c=log(ma)-mean(B);
    [a,ok]=gammaShape(c,theta(1));
    if ~ok || ~isfinite(ma) || ma<=0, fit.status=2; break; end
    next=[a,a/ma];
    % 先计算新参数的E-step，令似然和参数历史一一对应。
    % 这一结果下一轮直接复用，不需要重复积分。
    [An,Bn,lln,code]=estep(d,next,opt);
    if code~=0, fit.status=code; break; end
    if lln < ll-opt.monotoneTol*(1+abs(ll))
        fit.status=3; break;
    end
    change=norm(next-theta);
    theta=next; A=An; B=Bn; ll=lln;
    fit.theta=theta; fit.logLik=ll; fit.iter=r; fit.lastChange=change;
    if keepHistory, H(r+1,:)=[r,theta,ll,change]; end
    if change<opt.tol, fit.status=0; break; end
end
if keepHistory, fit.history=H(1:fit.iter+1,:); end
end

function [A,B,ll,status]=estep(d,theta,opt)
% E-step的高效数值实现：将w=u*z，利用指数倾斜Beta分布的正项级数。
% 这与原程序的卷积积分是同一个积分，不改变统计模型或EM更新。
% 若beta>=betaS, z=w/u；否则 z=x/u，保证指数级数的参数t>=0。
% 积分端点即使有奇异性（alpha*y<1或alphaS<1），也被Beta积分解析吸收。
% posterior = Beta(p+k,q)的混合，权重正比 t^k/k! * (p)_k/(p+q)_k。
% 因此可直接求 E[W]、E[log W]及卷积密度，无需固定Gauss-Legendre节点。
A=[]; B=[]; ll=NaN; status=2;
a=theta(1); b=theta(2);
if any(~isfinite(theta)|theta<=0), return; end
if a>opt.maxShape, status=4; return; end
n=numel(d.U); s=d.aX+a;
t=abs(d.beta-b).*d.U;
% 级数截断：泊松尾部的保守上界；混合项的尾部不大于Poisson(t)的尾部。
% 先检查最大项数，超出时明确标记，而不是悄悄修改估计值。
K=ceil(max(t)+12*sqrt(max(t)+1)+12);
if ~isfinite(K)||K>opt.maxTerms, status=4; return; end
if d.beta>=b, p=a*ones(n,1); q=d.aX; wIsZ=true;
else, p=d.aX; q=a*ones(n,1); wIsZ=false; end
% 向量化构造权重的对数，避免指数或Gamma函数溢出。
k=0:K-1;
lr=log(t)+log(p+k)-log(s+k)-log(k+1);
logw=[zeros(n,1),cumsum(lr,2)];
m=max(logw,[],2); weights=exp(logw-m); den=sum(weights,2);
prob=weights./den;
% digamma递推 psi(x+k)=psi(x)+sum_{j=0}^{k-1}1/(x+j)。
psiS=psi(s)+[zeros(n,1),cumsum(1./(s+k),2)];
if wIsZ
    zmean=(p+(0:K))./(s+(0:K));
    zlog=psi(p)+[zeros(n,1),cumsum(1./(p+k),2)]-psiS;
else
    zmean=q./(s+(0:K));
    zlog=psi(q)-psiS;
end
A=d.U.*sum(prob.*zmean,2);
B=d.logU+sum(prob.*zlog,2);
% 这里的 ll 是观测数据对数似然，不是Q函数。
lg=d.base+a*log(b)+a.*d.logU-gammaln(s) ...
    -max(d.beta,b).*d.U+m+log(den);
ll=sum(lg);
if all(isfinite(A)&A>0&isfinite(B))&&isfinite(ll), status=0; end
end

function T=validateEstimator(par,opt)
% 数值自检：等rate时 W/U 的后验为Beta(aS,aX)，有解析解；
% 不等rate时，使用MATLAB自带的自适应积分独立核对。
% 选取积分端点平滑的参数，避免自适应积分本身的端点误差。
U=3.7; Y=2.4; aX=par.alpha*Y; a=2.6;
d=struct('U',U,'logU',log(U),'aX',aX,'beta',par.beta);
d.base=aX*log(par.beta)+(aX-1)*log(U);
rates=[par.beta,.6,3.0]; errors=zeros(numel(rates),3);
for k=1:numel(rates)
    b=rates(k); [A,B,ll,st]=estep(d,[a,b],opt); assert(st==0);
    if b==par.beta
        Ae=U*a/(a+aX); Be=log(U)+psi(a)-psi(a+aX);
        lle=(a+aX)*log(b)-gammaln(a+aX)+(a+aX-1)*log(U)-b*U;
    else
        f=@(z) exp((a-1).*log(max(z,realmin)) + ...
            (aX-1).*log(max(1-z,realmin)) + (par.beta-b)*U.*z);
        den=integral(f,0,1,'AbsTol',1e-12,'RelTol',1e-11);
        Ae=U*integral(@(z)z.*f(z),0,1,'AbsTol',1e-12,'RelTol',1e-11)/den;
        Be=log(U)+integral(@(z)log(max(z,realmin)).*f(z),0,1, ...
            'AbsTol',1e-12,'RelTol',1e-11)/den;
        lle=aX*log(par.beta)+a*log(b)+(a+aX-1)*log(U)- ...
            gammaln(aX)-gammaln(a)-par.beta*U+log(den);
    end
    errors(k,:)=abs([A-Ae,B-Be,ll-lle]);
end
assert(max(errors(:))<1e-7,'E-step self-test failed.');
c=log(a)-psi(a); [ar,ok]=gammaShape(c,1); assert(ok&&abs(ar-a)<1e-8);
T=table(rates(:),errors(:,1),errors(:,2),errors(:,3), ...
    'VariableNames',{'BetaS','AbsErrorA','AbsErrorB','AbsErrorLogDensity'});
end

function [a,ok]=gammaShape(c,a0)
% 求解 log(a)-psi(a)=c>0。单调严格递减，唯一正根。
% 使用带区间保护的Newton迭代；避免负形状参数和无约束Newton跳出定义域。
a=NaN; ok=false;
if ~isfinite(c)||c<=0, return; end
lo=0.5/c; hi=1/c;  % 1/(2a)<log(a)-psi(a)<1/a，故根落在此区间
if ~isfinite(hi)||lo<=0, return; end
a=min(max(a0,lo),hi);
for k=1:100
    [h,hp]=shapeFun(a); f=h-c;
    if abs(f)<=1e-12*c, ok=true; return; end
    if f>0, lo=a; else, hi=a; end
    next=a-f/hp;
    if ~isfinite(next)||next<=lo||next>=hi, next=lo+(hi-lo)/2; end
    if abs(next-a)<=1e-12*max(a,realmin), a=next; ok=true; return; end
    a=next;
end
end

function [h,hp]=shapeFun(a)
% 大形状参数下 log(a)-psi(a) 存在消去误差，使用渐近式提高精度。
if a>1e4
    x=1/a; h=x/2+x^2/12-x^4/120+x^6/252;
    hp=-x^2/2-x^3/6+x^5/30-x^7/42;
else
    h=log(a)-psi(a); hp=1/a-psi(1,a);
end
end

%% ====================== 外层MC、内层Bootstrap =====================
function run=oneMC(n,par,truth,cfg,seed)
rng(seed,'twister');
[U,Y]=simulateData(n,par,truth);
f=fitEM(U,Y,par,cfg.em,[],false);
run=struct('theta',f.theta,'status',f.status,'iter',f.iter,'logLik',f.logLik, ...
    'fallback',f.fallback,'se',[NaN,NaN],'ci',NaN(2,2),'bootOK',0, ...
    'bootStatusCounts',zeros(1,5),'bootTheta',NaN(cfg.B,2),'bootStatus',NaN(cfg.B,1));
if f.status~=0, return; end
% 这里固定原始Y，并用估计出的冲击Gamma分布生成B组数据。
% 不能在每次bootstrap中使用真实参数truth。
for b=1:cfg.B
    Xb=randg(par.alpha.*Y)/par.beta;
    Wb=randg(f.theta(1)*ones(n,1))/f.theta(2);
    fb=fitEM(Xb+Wb,Y,par,cfg.em,[],false);
    run.bootTheta(b,:)=fb.theta; run.bootStatus(b)=fb.status;
end
valid=run.bootStatus==0;
run.bootOK=sum(valid);
for k=0:4, run.bootStatusCounts(k+1)=sum(run.bootStatus==k); end
% 不将失败的拟合当成成功，也不自动重抽直到成功（会掩盖失败比例）。
% 若有失败，SE和CI是基于有效重复的条件估计；必须报告有效率。
if run.bootOK>=max(2,ceil(cfg.minBootFraction*cfg.B))
    v=run.bootTheta(valid,:);
    run.se=std(v,0,1);
    run.ci=[empQuantile(v,.025);empQuantile(v,.975)];
end
end

function q=empQuantile(x,p)
% 明确使用线性插值分位数：位置 1+(N-1)*p；无需额外工具箱。
x=sort(x,1); h=1+(size(x,1)-1)*p; j=floor(h); k=ceil(h);
q=x(j,:)+(h-j).*(x(k,:)-x(j,:));
end

function [T,D]=summarizeMC(runs,n,truth,B)
R=numel(runs); est=vertcat(runs.theta); ok=[runs.status]'==0;
se=vertcat(runs.se); ciOK=all(isfinite(se),2); N=sum(ok); C=sum(ciOK);
T=table();
for j=1:2
    v=est(ok,j); [avg,bias,variance,rmse,sd]=stats(v,truth(j));
    cover=false(R,1); widths=NaN(R,1);
    for r=find(ciOK)'
        bounds=runs(r).ci(:,j);
        cover(r)=bounds(1)<=truth(j)&&truth(j)<=bounds(2);
        widths(r)=bounds(2)-bounds(1);
    end
    CP=NaN; if C>0, CP=100*sum(cover)/C; end
    % CP以有效CI数为分母；同时输出全R分母的保守覆盖率及CI有效率。
    name={'alphaS','betaS'};
    row=table(n,name(j),truth(j),avg,bias,variance,sd,rmse,mean(se(ciOK,j)), ...
        mean(widths(ciOK)),CP,R,N,C,100*N/R,100*C/R,100*sum(cover)/R, ...
        'VariableNames',{'n','Parameter','TrueValue','MeanEstimate','Bias','Variance','EmpiricalSD', ...
        'RMSE','MeanSE','MeanCIWidth','CP_pct','R','NConverged','NValidCI','Convergence_pct', ...
        'ValidCI_pct','CoverageAllR_pct'});
    T=appendTable(T,row);
end
sc=[runs.status]'; bc=vertcat(runs.bootStatusCounts); nBoot=sum(bc,1);
D=table(n,R,N,sum(sc==1),sum(sc==2),sum(sc==3),sum(sc==4),sum([runs.fallback]), ...
    N*B,nBoot(1),nBoot(2),nBoot(3),nBoot(4),nBoot(5),C,mean([runs(ok).iter]), ...
    'VariableNames',{'n','R','Converged','MaxIter','NumericalFailure','LikelihoodDecrease','GuardLimit', ...
    'MoMFallback','BootstrapAttempted','BootstrapConverged','BootstrapMaxIter','BootstrapNumericalFailure', ...
    'BootstrapLikelihoodDecrease','BootstrapGuardLimit','ValidCI','MeanIterationsConverged'});
end

function T=replicationDetails(runs,n)
% 逐次输出原始估计、置信区间和收敛状态，便于检查异常样本。
R=numel(runs); est=vertcat(runs.theta); se=vertcat(runs.se);
lo=NaN(R,2); hi=lo;
for r=1:R, lo(r,:)=runs(r).ci(1,:); hi(r,:)=runs(r).ci(2,:); end
T=table(repmat(n,R,1),(1:R)',est(:,1),est(:,2),[runs.logLik]', ...
    [runs.iter]',[runs.status]',[runs.fallback]',se(:,1),se(:,2), ...
    lo(:,1),lo(:,2),hi(:,1),hi(:,2),[runs.bootOK]', ...
    'VariableNames',{'n','Rep','AlphaHat','BetaHat','LogLikelihood','Iterations','Status', ...
    'MoMFallback','SE_alpha','SE_beta','Lo_alpha','Lo_beta','Hi_alpha','Hi_beta','BootstrapOK'});
end

function [avg,bias,v,rmse,sd]=stats(x,t)
% 只接受成功估计；若样本不足明确输出NaN。
avg=NaN; bias=NaN; v=NaN; rmse=NaN; sd=NaN;
if isempty(x), return; end
avg=mean(x); bias=avg-t; rmse=sqrt(mean((x-t).^2));
if numel(x)>=2, v=var(x,0); sd=sqrt(v); end
end

%% ========================== 两阶段实验 =============================
function run=oneTwoStage(par,truth,cfg,seed)
% 保留旧模式：独立生成新的冲击观测，并从中计算Known估计。
rng(seed,'twister');
[U,Y]=simulateData(cfg.nTwo,par,truth);
f=fitEM(U,Y,par,cfg.em,[],false);
% 使用独立随机种子生成无冲击样本，避免依赖fitEM的随机数调用。
rng(seedFor(seed,4,1,1),'twister');
run=fitTwoStageOnData(U,Y,f.theta,f.status,par,cfg);
end

function run=oneTwoStageReused(par,truth,cfg,jTwo,r,ref)
% 【新增】重现Table 2第r组n=cfg.nTwo的原始冲击数据(U,Y)，
% Known参数估计直接来自已保存的MC结果，不再运行Known EM。
rng(seedFor(cfg.seed,1,jTwo,r),'twister');
[U,Y]=simulateData(cfg.nTwo,par,truth);
% 生成一组完全独立的无冲击固定间隔退化增量，用于第一阶段MLE。
% 注意：绝对不能把U=X+W直接当成连续退化X使用。
rng(seedFor(cfg.seed,4,jTwo,r),'twister');
run=fitTwoStageOnData(U,Y,ref.theta,ref.status,par,cfg);
end

function run=fitTwoStageOnData(U,Y,knownTheta,knownStatus,par,cfg)
% Known和Plug-in严格共用U,Y；同一重复内不同m采用嵌套无冲击样本。
clean=randg(par.alpha*cfg.dtClean*ones(max(cfg.mList),1))/par.beta;
mCount=numel(cfg.mList);
% 保存恢复出的冲击数据与新增的无冲击数据，方便核验及后续复用。
run=struct('knownTheta',knownTheta,'knownStatus',knownStatus, ...
    'U',U,'Y',Y,'cleanData',clean, ...
    'continuous',NaN(mCount,2),'plugTheta',NaN(mCount,2),'plugStatus',2*ones(mCount,1));
for j=1:mCount
    x=clean(1:cfg.mList(j));
    % 无冲击退化数据的Gamma MLE：shape=alpha*dt, rate=beta。
    mx=mean(x); c=log(mx)-mean(log(x));
    [shape,ok]=gammaShape(c,max(mean(x)^2/var(x,0),.01));
    if ~ok || ~isfinite(mx) || mx<=0, continue; end
    run.continuous(j,:)=[shape/cfg.dtClean,shape/mx];
    pp=par; pp.alpha=run.continuous(j,1); pp.beta=run.continuous(j,2);
    g=fitEM(U,Y,pp,cfg.em,[],false);
    run.plugTheta(j,:)=g.theta; run.plugStatus(j)=g.status;
end
end

function verifyTable2Replay(refRuns,par,truth,cfg,jTwo,Rtwo)
% 为防止模型、随机种子或simulateData的实现变化造成错误复用，
% 抽查首/中/末三组数据：重新计算Known并对比Table 2已保存估计。
testReps=unique([1,ceil(Rtwo/2),Rtwo]);
for h=1:numel(testReps)
    r=testReps(h);
    rng(seedFor(cfg.seed,1,jTwo,r),'twister');
    [U,Y]=simulateData(cfg.nTwo,par,truth);
    f=fitEM(U,Y,par,cfg.em,[],false);
    saved=refRuns(r);
    err=norm(f.theta(:)-saved.theta(:));
    scale=max(1,norm(saved.theta(:)));
    if f.status~=saved.status || ~(isfinite(err)&&err<=1e-7*scale)
        error(['无法可靠复现Table 2的第%d组数据。请检查保存的checkpoint_MC.mat、', ...
               '随机种子、MATLAB版本及simulateData是否与原运行一致。'],r);
    end
end
fprintf('  Table 2随机种子复现检查通过（抽查 %d 组）。\n',numel(testReps));
end

function results=rerunTwoStageFromSavedMC(outDir,overrides)
% 直接从已完成的Table 2检查点执行两阶段实验，不重新执行Table 2与Bootstrap。
checkpoint=fullfile(outDir,'checkpoint_MC.mat');
assert(exist(checkpoint,'file')==2, ...
    '找不到 %s。请把outDir设置为原来包含checkpoint_MC.mat的结果文件夹。',checkpoint);
S=load(checkpoint,'mc','cfg','par','truth');
cfg=S.cfg; par=S.par; truth=S.truth; mc=S.mc;
if ~isfield(cfg,'reuseTwoFromMC'), cfg.reuseTwoFromMC=true; end
cfg.reuseTwoFromMC=true;
% 旧checkpoint可能来自nTwo=150的实验；本任务默认要复用Table 2的n=300。
cfg.nTwo=300;
% 允许修改两阶段实验参数，但必须保留原Table 2的数据生成设置与EM方法。
allowed={'nTwo','Rtwo','mList','dtClean','useParallel','makePlots'};
fn=fieldnames(overrides);
for k=1:numel(fn)
    if ~any(strcmp(fn{k},allowed))
        error('tworeuse模式不允许修改%s（需保持原Table 2的模型、种子和EM设置）。',fn{k});
    end
    cfg.(fn{k})=overrides.(fn{k});
end
jTwo=find(cfg.nList==cfg.nTwo,1,'first');
assert(~isempty(jTwo), 'checkpoint_MC.mat中没有n=%d的Table 2结果。',cfg.nTwo);
refRuns=mc{jTwo};
if ~isfield(overrides,'Rtwo'), cfg.Rtwo=numel(refRuns); end
assert(cfg.Rtwo>=2 && cfg.Rtwo<=numel(refRuns),'Rtwo必须介于2和现有MC重复次数之间。');
assert(all(cfg.mList>=2) && all(mod(cfg.mList,1)==0) && cfg.dtClean>0);
usePar=cfg.useParallel && license('test','Distrib_Computing_Toolbox');
if cfg.useParallel && ~usePar, warning('PCT不可用，改用串行。'); end
if usePar && isempty(gcp('nocreate')), parpool; end
verifyTable2Replay(refRuns,par,truth,cfg,jTwo,cfg.Rtwo);
reuseDir=fullfile(outDir,sprintf('two_stage_reuse_n%d',cfg.nTwo));
if ~exist(reuseDir,'dir'), mkdir(reuseDir); end
fprintf('\n只重算Table 4：复用Table 2 n=%d的%d组冲击数据，无需Bootstrap。\n', ...
    cfg.nTwo,cfg.Rtwo);
timer=tic; twoRuns=cell(cfg.Rtwo,1);
if usePar
    parfor r=1:cfg.Rtwo
        twoRuns{r}=oneTwoStageReused(par,truth,cfg,jTwo,r,refRuns(r));
    end
else
    for r=1:cfg.Rtwo
        twoRuns{r}=oneTwoStageReused(par,truth,cfg,jTwo,r,refRuns(r));
        if mod(r,max(1,floor(cfg.Rtwo/5)))==0 || r==cfg.Rtwo
            fprintf('  %d/%d, %.1f s\n',r,cfg.Rtwo,toc(timer));
        end
    end
end
[t4,t4diag,t4first]=summarizeTwo(twoRuns,cfg,truth,par);
writetable(t4,fullfile(reuseDir,'Table4_two_stage_performance.csv'));
writetable(t4diag,fullfile(reuseDir,'Table4_diagnostics.csv'));
writetable(t4first,fullfile(reuseDir,'Table4_first_stage_estimates.csv'));
results=struct('config',cfg,'parameters',par,'truth',truth, ...
    'Table4',t4,'Table4Diagnostics',t4diag,'Table4FirstStage',t4first, ...
    'twoStageRuns',{twoRuns},'sourceCheckpoint',checkpoint,'reusedN',cfg.nTwo);
save(fullfile(reuseDir,'two_stage_reuse_results.mat'),'results','-v7.3');
if cfg.makePlots, makeTwoStageFigure(t4,reuseDir); end
totalTime=toc(timer);
writetable(table(totalTime,totalTime/60, ...
    'VariableNames',{'Seconds','Minutes'}), ...
    fullfile(reuseDir,'Computational_time_two_stage.csv'));
fprintf('\n新两阶段结果保存于：%s\n',reuseDir);
disp(t4);
end

function [T,D,F]=summarizeTwo(runs,cfg,truth,par)
R=numel(runs); known=vertcatRuns(runs,'knownTheta'); knownOK=cellfun(@(x)x.knownStatus==0,runs);
T=table(); D=table(); F=table(); names={'alphaS','betaS'}; cnames={'alpha','beta'};
for j=1:numel(cfg.mList)
    m=cfg.mList(j); plug=NaN(R,2); cont=NaN(R,2); st=NaN(R,1);
    for r=1:R
        plug(r,:)=runs{r}.plugTheta(j,:); cont(r,:)=runs{r}.continuous(j,:); st(r)=runs{r}.plugStatus(j);
    end
    paired=knownOK & st==0; N=sum(paired);
    % 用两种方法都收敛的相同重复进行配对比较，避免比较不同数据子集。
    for k=1:2
        for method=1:2
            if method==1, v=known(paired,k); label={'Known'};
            else, v=plug(paired,k); label={'PlugIn'}; end
            [avg,bias,variance,rmse,sd]=stats(v,truth(k));
            row=table(m,cfg.nTwo,label,names(k),truth(k),avg,bias,variance,sd,rmse,R,N, ...
                'VariableNames',{'m','n','Method','Parameter','TrueValue','MeanEstimate','Bias', ...
                'Variance','EmpiricalSD','RMSE','R','NPaired'});
            T=appendTable(T,row);
        end
        t=[par.alpha,par.beta]; good=all(isfinite(cont),2);
        [avg,bias,variance,rmse,sd]=stats(cont(good,k),t(k));
        F=appendTable(F,table(m,cnames(k),t(k),avg,bias,variance,sd,rmse,sum(good), ...
            'VariableNames',{'m','Parameter','TrueValue','MeanEstimate','Bias','Variance','EmpiricalSD','RMSE','NValid'}));
    end
    D=appendTable(D,table(m,R,sum(knownOK),sum(st==0),N,sum(st==1),sum(st==2),sum(st==3),sum(st==4), ...
        'VariableNames',{'m','R','KnownConverged','PlugInConverged','BothConverged','PlugInMaxIter', ...
        'PlugInNumericalFailure','PlugInLikelihoodDecrease','PlugInGuardLimit'}));
end
end

function x=vertcatRuns(runs,field)
x=NaN(numel(runs),2); for r=1:numel(runs), x(r,:)=runs{r}.(field); end
end

function T=appendTable(T,row)
if isempty(T), T=row; else, T=[T;row]; end
end

function s=statusText(code)
labels={'Converged','MaxIter','NumericalFailure','LikelihoodDecrease','GuardLimit'};
s=labels{code+1};
end

%% =============================== 绘图 ===============================
function makeFigures(t2,t3,fits,t4,outDir)
% 绘图兼容性修复：
% 1) 不再使用 subplot()：部分 MATLAB 版本的 subplot 监听器会报告“对象无效或已删除”。
% 2) 不再直接调用 print()：避免 alternatePrintPath 图形输出异常。
% 3) 使用显式 axes + exportgraphics()，同时保留可编辑的 .fig 文件。
% 4) 三张图独立处理：某张图保存失败，不影响其他图或已保存的统计表。
% 5) 如在图形界面中运行，可将 'Visible','off' 改为 'Visible','on'。

% -------- 图 1：不同初值下的 EM 迭代轨迹 --------
fig=[];
try
    fig=figure('Visible','off','Color','w','Position',[100 100 950 760]);
    labels={'\alpha_S','\beta_S','Observed log-likelihood'};
    for k=1:3
        % 显式设置坐标轴位置，不调用 subplot()。
        ax=axes('Parent',fig,'Units','normalized', ...
                'Position',[0.12,0.715-(k-1)*0.316,0.82,0.245]);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');
        for j=1:numel(fits)
            H=fits{j}.history;
            if ~isempty(H)
                plot(ax,H(:,1),H(:,k+1),'LineWidth',1.1);
            end
        end
        ylabel(ax,labels{k}); xlabel(ax,'EM iteration');
        if k==1
            legend(ax,t3.Initialization,'Location','best','Interpreter','none');
        end
    end
    saveBothFormats(fig,outDir,'Figure_EM_initialization');
catch ME
    warning('EM 迭代图保存失败：%s',ME.message);
end
closeIfValid(fig);

% -------- 图 2：不同样本量的 EM 估计 RMSE --------
fig=[];
try
    fig=figure('Visible','off','Color','w','Position',[100 100 920 400]);
    names={'alphaS','betaS'};
    for k=1:2
        ax=axes('Parent',fig,'Units','normalized', ...
                'Position',[0.10+(k-1)*0.50,0.17,0.36,0.73]);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');
        ix=strcmp(t2.Parameter,names{k});
        plot(ax,t2.n(ix),t2.RMSE(ix),'-o','LineWidth',1.4);
        xlabel(ax,'Sample size n'); ylabel(ax,'RMSE'); title(ax,names{k});
    end
    saveBothFormats(fig,outDir,'Figure_EM_RMSE');
catch ME
    warning('EM RMSE 图保存失败：%s',ME.message);
end
closeIfValid(fig);

% -------- 图 3：两阶段估计的 RMSE 对比 --------
makeTwoStageFigure(t4,outDir);
end

function makeTwoStageFigure(t4,outDir)
% 新增的tworeuse模式可直接调用，不必再次绘制Table 2和初值敏感性图。
fig=[];
try
    fig=figure('Visible','off','Color','w','Position',[100 100 920 400]);
    names={'alphaS','betaS'};
    for k=1:2
        ax=axes('Parent',fig,'Units','normalized', ...
                'Position',[0.10+(k-1)*0.50,0.17,0.36,0.73]);
        hold(ax,'on'); box(ax,'on'); grid(ax,'on');
        ixKnown=strcmp(t4.Parameter,names{k}) & strcmp(t4.Method,'Known');
        ixPlug=strcmp(t4.Parameter,names{k}) & strcmp(t4.Method,'PlugIn');
        plot(ax,t4.m(ixKnown),t4.RMSE(ixKnown),'-o','LineWidth',1.4);
        plot(ax,t4.m(ixPlug),t4.RMSE(ixPlug),'-o','LineWidth',1.4);
        xlabel(ax,'Shock-free sample size m');
        ylabel(ax,'RMSE'); title(ax,names{k});
        legend(ax,{'Known','Plug-in'},'Location','best');
    end
    saveBothFormats(fig,outDir,'Figure_two_stage_RMSE');
catch ME
    warning('两阶段 RMSE 图保存失败：%s',ME.message);
end
closeIfValid(fig);
end

function saveBothFormats(fig,outDir,baseName)
% 先保存 MATLAB 可编辑 FIG 文件，再导出 200 dpi PNG。
% exportgraphics 为 R2020a 起提供的官方图形导出函数。
figFile=fullfile(outDir,[baseName,'.fig']);
pngFile=fullfile(outDir,[baseName,'.png']);
savefig(fig,figFile);
if exist('exportgraphics','file')==2
    exportgraphics(fig,pngFile,'Resolution',200);
else
    % 旧版 MATLAB 不支持 exportgraphics 时备用。
    % 尽量显式指定 OpenGL 渲染；如本机图形驱动仍报错，可仅保留 FIG 输出。
    warning('当前 MATLAB 不支持 exportgraphics，尝试使用 saveas 导出 PNG。');
    saveas(fig,pngFile);
end
fprintf('已保存 %s.fig 和 %s.png\n',baseName,baseName);
end

function closeIfValid(fig)
% MATLAB 不同版本对已失效图形句柄的关闭操作不完全一致。
if isempty(fig), return; end
try
    if isgraphics(fig,'figure')
        close(fig);
    end
catch
    % 对图形监听器的清理异常不继续抛错，以免掩盖计算结果。
end
end
