function [x, y, u, p] = simulate_DPC(x0, dim, params, optim, IZ, LQ, bounds, e)
    %   Simulates the subspace data-driven predictive control given the
    %   reduced LQ decomposition matrices, the relevant indices of ZP and
    %   a reference signal
    %
    %   Inputs
    %           x0:         state initial condition
    %           dim:        system dimensions and horizons
    %           params:     struct containing system parameters
    %           optim:      struct containing optimization parameters
    %           IZ:         ordered selected indices of ZP
    %           LQ:         struct containing LQ decomposition matrices
    %           bounds:     struct containing saturation bounds
    %           e:          measurement noise vector
    % 
    %   Outputs
    %           x:          state of the system
    %           y:          measured output
    %           u:          control input
    %           p:          measured scheduling signal

    
    % Extract parameters
    pmaxo   = params.pmaxo;
    pmino   = params.pmino;
    m       = params.m;
    g       = params.g;
    l       = params.l;
    J       = params.J;
    taum    = params.taum;
    Km      = params.Km;
    Ts      = params.Ts;    
    
    % Extract dimensions
    M   = dim.M;
    T   = dim.T;
    nu  = dim.nu;
    ny  = dim.ny;
    nx  = dim.nx;
    np  = dim.np;
    
    % Extract simulation length
    nsim = optim.nsim;    
    
    
    %% Initialize signals
    p0 = 2 / (pmaxo - pmino) * (sinc(x0(1) / pi) - (pmaxo + pmino) / 2);

    % Boundary conditions
    y_init  = kron(x0(1), ones(1,M));
    u_init  = -(m*g*l / J * sin(y_init)) * (taum / Km);
    p_init  = [ones(1, M); 2 / (pmaxo-pmino) * (sinc(y_init/pi) - (pmaxo+pmino) / 2)];
    p_next  = repmat(p_init(:,end), 1, T);

    % Initialize function output
    u   = zeros(nu, nsim);
    y   = zeros(ny, nsim);
    p   = [ones(1, nsim); zeros(np-1, nsim)];
    x   = zeros(nx, nsim+1);

    p(2:np, 1)  = p0;
    x(:, 1)     = x0;

    
    %% Run simulation
    fprintf('\nStarting simulation with LPV DPC\n')
    for k = 1:nsim  
        % Initial conditions
        z_init = build_init_reduced(dim, u_init, y_init, p_init, p_next, IZ);

        % Compute non-variable elements
        gamma1  = LQ.L11_inv * z_init;
        u1      = LQ.L21 * gamma1;
        y1      = LQ.L31 * gamma1;

        % Set up and solve optimization problem
        fprintf('Solving optimization at time step %-10d', k)
        [u_opt, cost, info] = LPVgammaDDPC_reduced(dim, u1, y1, LQ.L22, LQ.L3, ...
                                            optim.y_ref(k:k+T-1), optim, bounds);
        fprintf('Done!\n')

        if info.problem ~= 0
            error('Problem status %d, issues arose during optimization', info.problem);
        end

        % Simulate system evolution
        x(:, k+1) = [x(1,k) + x(2,k)*Ts;
                     x(2,k) + (m*g*l/J*sin(x(1,k)) - 1/taum*x(2,k) + Km/taum*u_opt)*Ts];

        y(:,k)  = x(1,k) + e(k);
        p_curr  = [1; 2 / (pmaxo-pmino) * (sinc(y(:,k)/pi) - (pmaxo+pmino) / 2)];
        p(:,k)  = p_curr;
        u(:,k)  = u_opt;

        % Update scheduling parameters
        p_init  = [p_init(:, 2:end), p_curr];
        p_next  = repmat(p_curr, 1, T);

        % Update initial condition parameters
        u_init = [u_init(:, 2:end), u_opt];
        y_init = [y_init(:, 2:end), y(:,k)];
    end

end


