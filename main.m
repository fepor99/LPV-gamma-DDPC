%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                                                               %%%%%%
%%%%%%      Author: Federico Porcari                                 %%%%%%
%%%%%%      Code of the paper "A subspace approach to data-driven    %%%%%%
%%%%%%      predictive control for linear parameter-varying systems" %%%%%%
%%%%%%                                                               %%%%%%
%%%%%%      Required software:                                       %%%%%%
%%%%%%          - Signal Processing Toolbox                          %%%%%%
%%%%%%          - MATLAB Coder (optional but helps)                  %%%%%%
%%%%%%          - YALMIP                                             %%%%%%
%%%%%%          - GUROBI Solver                                      %%%%%%
%%%%%%                                                               %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc
clearvars
close all

rng(5)
addpath("LPV_gamma-DDPC/")
addpath("LPV-IO-DPC/")
addpath("LPV-MPC")
addpath("utils/")


%% Design parameters
% Monte Carlo runs
runs = 100;

% Noise standard deviation
noise.offline   = 0.01;                 % noise std during data collection
noise.online    = 0.01;                 % noise std during control testing
noise.process   = false;                % if true, apply process noise
noise.K{1}      = [0.77, 0.44]';        % or [0.93; 0.36] if K{2} is set to 0
noise.K{2}      = [-0.02, 0.76]';

% Data collection parameters
Ndata   = 89;                   % number of data collected
T       = 20;                   % future horizon length
M       = 2;                    % past horizon length
N       = Ndata - M - T + 1;

% Matrix reduction parameters
reduction.rows.z             = 10;          % rows of reduced ZP
reduction.rows.u             = 28;          % rows of reduced UF
reduction.mult_max_ZP        = 3;           % don't consider rows of ZP/UF with at least 
reduction.mult_max_UF        = 3;           % mult_max multiplications with scheduling signals
reduction.residual_method_ZP = "exp";       % ["lin" or "exp"] method to drop high-residual rows
reduction.residual_method_UF = "exp";       % ["lin" or "exp"] method to drop high-residual rows

% Optimization parameters
optim.nsim   = 100;                 % simulation length
optim.Q      = 16;                  % tracking weight
optim.R      = 0.01;                % control input weight
optim.beta2  = 3;                   % regularization parameter for gamma2
optim.beta3  = 0;                   % regularization parameter for gamma2
optim.Q_KF   = diag([1e-4, 10]);    % Kalman filter Q matrix
optim.R_KF   = 0.01;                % Kalman filter R matrix

% LPV-IO-DPC optimization parameters
lambda_sigma = 1e9; 
lambda_g     = 0.6;

% Reference signal
y_ref = zeros(optim.nsim+T, 1);
optim.y_ref = y_ref;       

% Initial conditions
xini    = [-pi/4, 5]';          % data collection
x0      = [-pi/2, 0]';          % control deployment

% Plot metadata
plt.normalize_p = true;        % if true, plot p normalized in [-1, 1]
plt.padd        = 0.3;
plt.font        = 20;

% Force C code generation of reduced LQ factorization
force_codegen = 0;


%% Initialize disk system
% Construct system
[params, bounds, sys, dim] = init_disk(noise);

% Extend dim struct with horizons
dim.T   = T;
dim.N   = N;
dim.M   = M;

% Check assumptions of the script
reduction.rows.z_full = (dim.nu * (dim.np * (dim.np + 1) / 2) + dim.ny * dim.np) * ...
                        dim.np^T * polyval(ones(1, M), dim.np * (dim.np + 1) / 2);
reduction.rows.u_full = dim.nu * polyval([ones(1, T), 0], dim.np);
check_assumptions(dim, reduction, optim, bounds)

% Instantiate signals
u           = zeros(dim.nu, Ndata, runs);           % offline dataset
x           = zeros(dim.nx, Ndata+1, runs);
y           = zeros(dim.ny, Ndata, runs);
p           = zeros(dim.np-1, Ndata, runs);
x_MPC       = zeros(dim.nx, optim.nsim+1, runs);    % MPC simulation
y_MPC       = zeros(dim.ny, optim.nsim, runs);
u_MPC       = zeros(dim.nu, optim.nsim, runs);
p_MPC       = zeros(dim.np, optim.nsim, runs);
x_gamma_DPC = zeros(dim.nx, optim.nsim+1, runs);    % gamma-DDPC simulation
y_gamma_DPC = zeros(dim.ny, optim.nsim, runs);
u_gamma_DPC = zeros(dim.nu, optim.nsim, runs);
p_gamma_DPC = zeros(dim.np, optim.nsim, runs);
y_IO_DPC    = zeros(dim.ny, optim.nsim+T, runs);    % LPV-IO-DPC simulation


%% Monte Carlo simulations
% Instantiate vector to save unstable realizations
unstable_IO_DPC     = [];
unstable_gamma_DPC  = [];

% Instantiate vector to save computational time
time_IO_DPC     = zeros(1, runs);
time_gamma_DPC  = zeros(1, runs);

% Perform Monte Carlo simulation
for i = 1:runs
    fprintf('\n\nSTARTING MONTE CARLO SIMULATION: %d / %d\n\n', i, runs)
    
    
            %%%%% DATA COLLECTION %%%%%
    u(:,:,i) = bounds.u_max * (1 - 2*rand(1, Ndata));
    [x(:,:,i), y(:,:,i), p(:,:,i)] = collect_data(u(:,:,i), xini, dim, params, sys, noise);
    plt_MC.collection.u = u(:,:,i);
    plt_MC.collection.y = y(:,:,i);
    plt_MC.collection.p = p(:,:,i);

    
            %%%%% LPV GAMMA-DDPC PREPROCESSING %%%%%
    % Perform LQ decomposition
    [IZ, IU, LQ, res] = setup_data_matrices(plt_MC, dim, reduction, force_codegen);
    
    % Plot residuals
    if i == 1
        plot_residuals(res, IZ, dim, reduction, plt)
    end
    
    % Generate noise
    e = noise.online * randn(dim.ny, optim.nsim);

            %%%%% CONTROL SIMULATION %%%%%
    % Compute MPC + KF trajectory
    [x_MPC(:,:,i), y_MPC(:,:,i), u_MPC(:,:,i), p_MPC(:,:,i)] = simulate_LPV_MPC(x0, dim, sys, params, optim, bounds, e);

    % Compute LPV gamma-DDPC trajectory
    tic
    try
        [x_gamma_DPC(:,:,i), y_gamma_DPC(:,:,i), u_gamma_DPC(:,:,i), p_gamma_DPC(:,:,i)] = simulate_DPC(x0, dim, params, optim, IZ, LQ, bounds, e);
    catch
        unstable_gamma_DPC = [unstable_gamma_DPC, i];
    end
    time_gamma_DPC(i) = toc;
    
    % Compute LPV-IO-DPC trajectory
    tic
    try
        y_IO_DPC(:,:,i) = LPV_IO_DPC(params, dim, bounds, sys, plt_MC, optim, lambda_sigma, lambda_g, e);
    catch
        unstable_IO_DPC = [unstable_IO_DPC, i];
    end
    time_IO_DPC(i) = toc;
end


%% Store all signals into struct for plotter
% Simulation details
plt.y_ref       = optim.y_ref(1:end-T+1);
plt.runs        = runs;
plt.noise       = noise;
plt.dim         = dim;
plt.reduction   = reduction;
plt.optim       = optim;

% Data collection
plt.collection.x = x;
plt.collection.y = y;
plt.collection.u = u;
plt.collection.p = p;

% MPC + KF
plt.MPC.x = x_MPC(:, 1:end-1, :);
plt.MPC.y = y_MPC;
plt.MPC.u = u_MPC;
plt.MPC.p = p_MPC(2, :, :);

% LPV gamma-DDPC
plt.gamma_DPC.x         = x_gamma_DPC(:, 1:end-1, :);
plt.gamma_DPC.y         = y_gamma_DPC;
plt.gamma_DPC.u         = u_gamma_DPC;
plt.gamma_DPC.p         = p_gamma_DPC(2, :, :);
plt.gamma_DPC.time      = time_gamma_DPC;
plt.gamma_DPC.unstable  = unstable_gamma_DPC;

% LPV-IO-DPC
plt.IO_DPC.y        = y_IO_DPC;
plt.IO_DPC.time     = time_IO_DPC;
plt.IO_DPC.unstable = unstable_IO_DPC;

% Set trajectories to remove from plotting
plt.IO_DPC.rem      = unstable_IO_DPC;
plt.gamma_DPC.rem   = unstable_gamma_DPC;


%% Plotting
plot_data_mean_std(plt)


    


