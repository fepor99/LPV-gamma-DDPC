function [x, y, p] = collect_data(u, xini, dim, params, sys, noise)
    %   Performs the data collection experiment given an initial condition
    %   for the unbalanced disk system and a specified control input for
    %   the whole simulation
    %
    %   Inputs
    %           u:          process noise state-space matrix
    %           xini:       state initial condition
    %           dim:        system dimensions
    %           params:     struct containing system parameters
    %           sys:        struct containing disk state-space matrices
    %           noise:      struct containing noise parameters
    % 
    %   Outputs
    %           x:          state of the system
    %           y:          measured output
    %           p:          measured scheduling signal


    % Extract scheduling bound (not normalized)
    pmaxo = params.pmaxo;
    pmino = params.pmino;
    
    % Extract system dimensions and horizons
    nx  = dim.nx;
    ny  = dim.ny;
    T   = dim.T;
    M   = dim.M;
    N   = dim.N;
    
    % Initialize LPV system
    fprintf('Generating data...')
    x       = zeros(nx, N+M+T);
    x(:,1)  = xini;

    % Noise realization
    e = noise.offline * randn(ny, N+M+T-1);
    
    % Simulate system
    for k = 1:N+M+T-1
        p_real      = 2/(pmaxo - pmino) * (sinc(x(1,k)/pi) - (pmaxo + pmino)/2);
        x(:, k+1)   = (sys.A{1} + sys.A{2} * p_real(1)) * x(:,k) + ...
                      (sys.B{1} + sys.B{2} * p_real(1)) * u(1,k) + ...
                      (sys.K{1} + sys.K{2} * p_real(1)) * e(:,k);
    end
    y = x(1,1:end-1) + e;
    p = 2 / (pmaxo - pmino) * (sinc(y / pi) - (pmaxo + pmino) / 2);
    fprintf('\tDone!\n\n')
    
    % Print signal-to-noise ratio
    fprintf('Signal-to-noise ratio: %.2fdB\n\n', snr(x(1,1:end-1), e))
end


