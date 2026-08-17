function [u_opt, cost, info] = LPVgammaDDPC_reduced(dim, u1, y1, L22, L3, y_ref, optim, bounds)
    %   Solves the optimization problem of the LPV gamma-DDPC
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           u1:         solution to u1 = L21 * gamma1
    %           y1:         solution to y1 = L32 * gamma1
    %           L22:        L22 matrix from LQ factorization
    %           L3:         [L32, L33] matrix from LQ factorization
    %           y_ref:      reference signal
    %           optim:      struct containing optimization parameters
    %           bounds:     input/output saturations
    %
    %   Outputs
    %           u_opt:      optimal first-step control input
    %           cost:       cost function value
    %           info:       optimization solution info

    
    % System dimensions
    nu  = dim.nu;
    T   = dim.T;
    
    % Optimization parameters
    Q       = optim.Q;
    R       = optim.R;
    beta2   = optim.beta2;
    beta3   = optim.beta3;
    
    % Signal values
    u_max   = bounds.u_max;
    u_min   = bounds.u_min;
    y_max   = bounds.y_max;
    y_min   = bounds.y_min;
    
    % Set optimizer options
    yalmip('clear');
    options = sdpsettings('solver', 'gurobi', 'verbose', 0);
                      
    % Optimization variables
    if optim.beta2 || ~optim.beta3
        gamma2 = sdpvar(height(L22), 1);
    else
        gamma2 = zeros(height(L22), 1);
    end
    
    if optim.beta3
        gamma3 = sdpvar(T, 1);
    else
        gamma3 = zeros(T, 1);
    end
    
    % Initialize cost with regularization (if present)
    Jcost = beta2 * (gamma2.' * gamma2) + beta3 * (gamma3.' * gamma3);
    
    % Get future inputs and outputs
	uf  = u1 + L22 * gamma2;
    yf  = y1 + L3 * [gamma2; gamma3];

    % Output saturation constraints
    cons = [yf <= repmat(y_max, T, 1)] + ...
           [yf >= repmat(y_min, T, 1)]; 
    
    for i = 1:T        
        % Update cost function
        Jcost = Jcost + (yf(i) - y_ref(i))' * Q * (yf(i) - y_ref(i)) ...
                      + uf((i-1)*nu+1:i*nu)' * R * uf((i-1)*nu+1:i*nu);
                  
        % Input saturation constraints
        cons = cons + [uf((i-1)*nu+1:i*nu) <= u_max] ...
                    + [uf((i-1)*nu+1:i*nu) >= u_min];
    end
    
    % Optimization
    info    = optimize(cons, Jcost, options);
    u_opt   = value(uf(1:nu));
    cost    = double(Jcost);   
    
end




