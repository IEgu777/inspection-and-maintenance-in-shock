clear; clc; close all;
rng(1);

%% ============================================================
%  1. Parameters from Section 6.1
%  注意：本文公式采用的是 Gamma(shape, rate) 参数化：
%  f(x)= beta^alpha/Gamma(alpha) * x^(alpha-1)*exp(-beta*x)
%  MATLAB 的 gamrnd/随机Gamma通常采用 scale 参数，所以生成时要用 1/rate
%% ============================================================

par.alpha  = 1.02;     % Gamma process shape rate per unit time
par.beta   = 1.20;     % Gamma process rate parameter
par.lambda = 2.50;     % shock arrival intensity
par.L      = 12.5;     % failure threshold

% p 和 q 是PM效果的Beta分布参数，第五章EM估计冲击幅值时暂时不用
par.p = 6;
par.q = 5;

%% ============================================================
%  2. True shock magnitude parameters used only for simulation
%  第6.1节没有给 alphaS 和 betaS，因为它们是待估计参数。
%  所以仿真时必须人为给定一组真值。
%% ============================================================

truePar.alphaS = 2.5;   % true shock Gamma shape
truePar.betaS  = 1;   % true shock Gamma rate

%% ============================================================
%  3. Simulate 100 degradation trajectories
%% ============================================================

nTraj = 200;

[traj, U, Y, Wtrue] = simulateShockDegradationTrajectories(nTraj, par, truePar);

fprintf('Number of trajectories: %d\n', nTraj);
fprintf('Number of shock-triggered observations: %d\n', numel(U));
fprintf('True alphaS = %.4f, true betaS = %.4f\n\n', truePar.alphaS, truePar.betaS);

%% ============================================================
%  4. EM estimation of shock magnitude distribution
%  Known: alpha, beta of continuous degradation
%  Unknown: alphaS, betaS of shock magnitude
%% ============================================================

opt.maxIter = 150;
opt.tol     = 1e-8;
opt.nQuad   = 200;       % Gauss-Legendre quadrature nodes

% Initial values
% 用一个简单的初值：U 减去连续退化期望 E[X_Y]=alpha/beta*Y
W0 = U - (par.alpha / par.beta) .* Y;
W0(W0 <= 0) = 0.1 * mean(U);

m0 = mean(W0);
v0 = var(W0, 1);

init.alphaS = max(m0^2 / v0, 0.1);
init.betaS  = max(m0 / v0,  0.1);

fprintf('Initial alphaS = %.4f, initial betaS = %.4f\n', init.alphaS, init.betaS);

[alphaHat, betaHat, history] = EM_estimateShockGamma(U, Y, par, init, opt);

fprintf('\n========== EM estimation result ==========\n');
fprintf('True alphaS      = %.6f\n', truePar.alphaS);
fprintf('Estimated alphaS = %.6f\n', alphaHat);
fprintf('True betaS       = %.6f\n', truePar.betaS);
fprintf('Estimated betaS  = %.6f\n', betaHat);
fprintf('Iterations       = %d\n', size(history,1));

%% ============================================================
%  5. Plot simulated degradation trajectories
%% ============================================================

figure; hold on; box on;

for k = 1:nTraj
    plot(traj(k).t, traj(k).S);
end

yline(par.L, 'r--', 'L');

xlabel('Time');
ylabel('Cumulative degradation level');
title('Simulated degradation trajectories');
grid on;

%% ============================================================
%  6. Plot EM convergence
%% ============================================================

figure; box on; hold on;
plot(history(:,1), history(:,2), '-o', 'LineWidth', 1.2);
plot(history(:,1), history(:,3), '-s', 'LineWidth', 1.2);
xlabel('Iteration');
ylabel('Parameter estimate');
legend('\alpha_S', '\beta_S', 'Location', 'best');
title('EM convergence');
grid on;


%% ============================================================
%  单独画一条：实际退化轨道 vs 检测到的退化轨道
%% ============================================================

dtPlot = 0.002;   % 用于显示真实轨道的时间细分步长，可适当调大或调小

oneTraj = simulateOneTrajectoryForPlot(par, truePar, dtPlot);

figure; hold on; box on;

% 实际退化轨道：包含连续退化 + 冲击瞬时跳变
plot(oneTraj.tActual, oneTraj.SActual, 'k-', 'LineWidth', 1.5);

% 检测到的退化轨道：只在冲击触发检查时刻观测
plot(oneTraj.tDetect, oneTraj.SDetect, 'ro-', ...
    'LineWidth', 1.2, ...
    'MarkerSize', 5, ...
    'MarkerFaceColor', 'r');

% 失效阈值
yline(par.L, 'r--', 'L', 'LineWidth', 1.2);

xlabel('Time');
ylabel('Cumulative degradation level');
legend('Actual degradation trajectory', ...
       'Detected degradation trajectory', ...
       'Failure threshold', ...
       'Location', 'best');

grid on;

%% ========================================================================
%                          Local functions
%% ========================================================================

function [traj, U_all, Y_all, W_all] = simulateShockDegradationTrajectories(nTraj, par, truePar)
% Simulate degradation trajectories observed at shock-triggered inspection points.
%
% For each interval:
%   Y_i ~ Exp(lambda)
%   X_{Y_i} ~ Gamma(alpha*Y_i, rate beta)
%   W_i ~ Gamma(alphaS, rate betaS)
%   U_i = X_{Y_i} + W_i
%
% The cumulative degradation is updated as:
%   S <- S + U_i
%
% The trajectory stops when S >= L.

    traj = struct([]);
    U_all = [];
    Y_all = [];
    W_all = [];

    for k = 1:nTraj

        t = 0;
        S = 0;

        t_path = 0;
        S_path = 0;

        while S < par.L

            % shock inter-arrival time
            y = -log(rand) / par.lambda;

            % continuous degradation increment over interval y
            % MATLAB uses scale parameter, so Gamma(shape, rate beta) = randg(shape)/beta
            x_inc = randg(par.alpha * y) / par.beta;

            % shock magnitude
            w_inc = randg(truePar.alphaS) / truePar.betaS;

            % total observed increment at shock-triggered inspection
            u = x_inc + w_inc;

            % update cumulative time and degradation
            t = t + y;
            S = S + u;

            % store trajectory
            t_path(end+1,1) = t;
            S_path(end+1,1) = S;

            % store shock-triggered observation
            Y_all(end+1,1) = y;
            U_all(end+1,1) = u;
            W_all(end+1,1) = w_inc;
        end

        traj(k).t = t_path;
        traj(k).S = S_path;
    end
end


function [alphaHat, betaHat, history] = EM_estimateShockGamma(U, Y, par, init, opt)
% EM estimation for shock magnitude W ~ Gamma(alphaS, rate betaS)
%
% Observed:
%   U_i = X_{Y_i} + W_i
%
% Known:
%   X_{Y_i} ~ Gamma(alpha*Y_i, rate beta)
%
% Unknown:
%   W_i ~ Gamma(alphaS, rate betaS)

    U = U(:);
    Y = Y(:);

    alphaS = init.alphaS;
    betaS  = init.betaS;

    [nodes, weights] = gaussLegendre01(opt.nQuad);

    history = zeros(opt.maxIter, 5);

    for r = 1:opt.maxIter

        % -------------------------
        % E-step
        % -------------------------
        [A, B, logLik] = Estep_compute_AB(U, Y, par, alphaS, betaS, nodes, weights);

        Abar = mean(A);
        Bbar = mean(B);

        % -------------------------
        % M-step
        % -------------------------
        % Solve:
        %   log(alphaS) - psi(alphaS) = log(Abar) - Bbar
        c = log(Abar) - Bbar;

        alphaNew = solveGammaShapeFromLogMoment(c, alphaS);
        betaNew  = alphaNew / Abar;

        relChange = norm([alphaNew, betaNew] - [alphaS, betaS]) / ...
                    (1 + norm([alphaS, betaS]));

        history(r,:) = [r, alphaNew, betaNew, logLik, relChange];

        alphaS = alphaNew;
        betaS  = betaNew;

        if relChange < opt.tol
            history = history(1:r,:);
            break;
        end
    end

    alphaHat = alphaS;
    betaHat  = betaS;
end


function [A, B, logLik] = Estep_compute_AB(U, Y, par, alphaS, betaS, nodes, weights)
% Compute:
%   A_i = E[W_i | U_i=u_i, Y_i=y_i]
%   B_i = E[log W_i | U_i=u_i, Y_i=y_i]
%
% Using Gauss-Legendre quadrature after transformation:
%   w = u * z, z in (0,1)
%
% The denominator is:
%   int_0^u f_X(u-w) f_W(w) dw
%
% The numerator for A_i is:
%   int_0^u w f_X(u-w) f_W(w) dw
%
% The numerator for B_i is:
%   int_0^u log(w) f_X(u-w) f_W(w) dw

    U = U(:);
    Y = Y(:);

    z  = nodes(:)';     % 1 x Q
    wz = weights(:)';   % 1 x Q

    Wpart = U .* z;          % n x Q
    Xpart = U .* (1 - z);    % n x Q

    aX = par.alpha .* Y;     % n x 1

    % log f_X(x)
    logFX = aX .* log(par.beta) ...
          - gammaln(aX) ...
          + (aX - 1) .* log(Xpart) ...
          - par.beta .* Xpart;

    % log f_W(w)
    logFW = alphaS .* log(betaS) ...
          - gammaln(alphaS) ...
          + (alphaS - 1) .* log(Wpart) ...
          - betaS .* Wpart;

    % because dw = U dz
    logIntegrand = log(U) + logFX + logFW;

    % log-sum-exp stabilization
    m = max(logIntegrand, [], 2);
    R = exp(logIntegrand - m);

    denom = sum(R .* wz, 2);

    A = sum(Wpart .* R .* wz, 2) ./ denom;
    B = sum(log(Wpart) .* R .* wz, 2) ./ denom;

    logLik = sum(m + log(denom));
end


function alpha = solveGammaShapeFromLogMoment(c, alpha0)
% Solve:
%   log(alpha) - psi(alpha) = c
%
% Newton iteration:
%   derivative = 1/alpha - psi(1, alpha)
%
% where psi(1, alpha) is the trigamma function.

    alpha = max(alpha0, 1e-4);

    for iter = 1:100

        f  = log(alpha) - psi(alpha) - c;
        fp = 1 / alpha - psi(1, alpha);

        step = f / fp;
        alphaNew = alpha - step;

        % keep positive
        count = 0;
        while (~isfinite(alphaNew) || alphaNew <= 0) && count < 20
            step = step / 2;
            alphaNew = alpha - step;
            count = count + 1;
        end

        if ~isfinite(alphaNew) || alphaNew <= 0
            alphaNew = alpha / 2;
        end

        if abs(alphaNew - alpha) / (1 + alpha) < 1e-12
            alpha = alphaNew;
            return;
        end

        alpha = alphaNew;
    end
end


function [x, w] = gaussLegendre01(N)
% Gauss-Legendre quadrature nodes and weights on [0,1]

    i = 1:N-1;
    a = i ./ sqrt(4*i.^2 - 1);

    CM = diag(a,1) + diag(a,-1);
    [V,D] = eig(CM);

    x = diag(D);
    [x,idx] = sort(x);

    w = 2 * (V(1,idx).^2)';

    % transform from [-1,1] to [0,1]
    x = (x + 1) / 2;
    w = w / 2;

    x = x(:);
    w = w(:);
end

function oneTraj = simulateOneTrajectoryForPlot(par, truePar, dtPlot)
% Simulate one degradation trajectory for visualization.
%
% Actual trajectory:
%   continuous Gamma degradation between shocks
%   + instantaneous shock jumps at shock arrival times.
%
% Detected trajectory:
%   only observed at shock-triggered inspection times.

    t = 0;
    S = 0;

    tActual = 0;
    SActual = 0;

    tDetect = 0;
    SDetect = 0;

    while S < par.L

        % 下一次冲击到达时间间隔
        y = -log(rand) / par.lambda;

        % 将冲击间隔细分，用于画出“实际连续退化轨道”
        nSub = max(2, ceil(y / dtPlot));
        dt = y / nSub;

        for j = 1:nSub

            % Gamma process increment:
            % X_{t+dt}-X_t ~ Gamma(alpha*dt, rate beta)
            x_inc = randg(par.alpha * dt) / par.beta;

            t = t + dt;
            S = S + x_inc;

            tActual(end+1,1) = t;
            SActual(end+1,1) = S;

            % 如果连续退化已经达到失效阈值，也可以停止
            if S >= par.L
                tDetect(end+1,1) = t;
                SDetect(end+1,1) = S;
                break;
            end
        end

        if S >= par.L
            break;
        end

        % 冲击瞬时退化增量
        w_inc = randg(truePar.alphaS) / truePar.betaS;

        % 为了在图中显示竖直跳变，先记录冲击前状态，再记录冲击后状态
        tActual(end+1,1) = t;
        SActual(end+1,1) = S;

        S = S + w_inc;

        tActual(end+1,1) = t;
        SActual(end+1,1) = S;

        % 冲击触发检查，因此检测点是冲击后的累计退化水平
        tDetect(end+1,1) = t;
        SDetect(end+1,1) = S;
    end

    oneTraj.tActual = tActual;
    oneTraj.SActual = SActual;
    oneTraj.tDetect = tDetect;
    oneTraj.SDetect = SDetect;
end