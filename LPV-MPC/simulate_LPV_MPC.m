function [x, y, u, p] = simulate_LPV_MPC(x0, dim, sys, params, optim, bounds, e)
    %   Simulates the LPV-MPC controller given an initial measurement of
    %   the disk position and a reference to be tracked
    %
    %   Inputs
    %           x0:         state initial condition
    %           dim:        system dimensions and horizons
    %           sys:        struct containing system matrices
    %           params:     struct containing system parameters
    %           optim:      struct containing optimization parameters
    %           bounds:     struct containing saturation bounds
    %           e:          measurement noise vector
    % 
    %   Outputs
    %           x:          state of the system
    %           y:          measured state
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
    T   = dim.T;
    nu  = dim.nu;
    nx  = dim.nx;
    ny  = dim.ny;
    np  = dim.np;
    
    % Extract simulation length and KF parameters
    nsim = optim.nsim;    
    Qk   = optim.Q_KF;
    Rk   = optim.R_KF;
    
    
    %% Initialize signals
    % Initialize function output
    u       = zeros(nu, nsim);
    y       = zeros(ny, nsim);
    p       = [ones(1, nsim); zeros(np-1, nsim)];
    x       = zeros(nx, nsim+1);
    x(:, 1) = x0;
    
    % Kalman Filter initialization
    x_hat   = x0;               % initial state estimate
    P_hat   = eye(nx);          % initial estimation error covariance

    
    %% Run simulation
    fprintf('\nStarting simulation with LPV MPC\n')
    for k = 1:nsim  
        
        % Get current scheduling matrices
        Ap = sys.A{1} + sys.A{2} * p(1,k);
        Bp = sys.B{1} + sys.B{2} * p(1,k);
        Cp = sys.C{1} + sys.C{2} * p(1,k);
        
        % Get noisy measurement
        y(:, k) = x(1, k) + e(:, k);
        
        % Kalman filter update
        if k > 1
            % KF prediction step
            x_hat_pred  = Ap * x_hat + Bp * u(:, k-1);
            P_pred      = Ap * P_hat * Ap' + Qk;

            % KF update step
            S       = Cp * P_pred * Cp' + Rk;
            Kf      = P_pred * Cp' / S;
            x_hat   = x_hat_pred + Kf * (y(:,k) - Cp * x_hat_pred);
            P_hat   = (eye(nx) - Kf * Cp) * P_pred;
        end
        
        % Update current noisy scheduling signal
        p_curr  = [1; 2 / (pmaxo-pmino) * (sinc(y(1, k)/pi) - (pmaxo+pmino) / 2)];
        p(:, k) = p_curr;
        
        % Update future scheduling parameters
        p_next = repmat(p_curr, 1, T);

        % Set up and solve optimization problem
        fprintf('Solving optimization at time step %-10d', k)
        [u_opt, ~, cost, info] = LPV_MPC(x_hat, dim, sys, optim.y_ref(k+1:k+T), p_next, optim, bounds);
        fprintf('Done!\n')

        if info.problem ~= 0
            fprintf('Problem status %d, issues arose during optimization\n', info.problem);
        end

        % Simulate system evolution
        x(:, k+1) = [x(1,k) + x(2,k)*Ts;
                     x(2,k) + (m*g*l/J*sin(x(1,k)) - 1/taum*x(2,k) + Km/taum*u_opt)*Ts];
        u(:,k)  = u_opt;

    end

end


