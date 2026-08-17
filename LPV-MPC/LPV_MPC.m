function [u_opt, y_fut, cost, info] = LPV_MPC(x0, dim, sys, y_ref, p_fut, optim, bounds)
    %   Solves the optimization problem of the LPV-MPC controller
    % 
    %   Inputs
    %           x0:         current measured state
    %           dim:        dimensions and horizons
    %           sys:        struct containing system matrices
    %           y_ref:      reference signal
    %           p_fut:      estimated future scheduling
    %           optim:      struct containing optimization parameters
    %           bounds:     input/output saturations
    %
    %   Outputs
    %           u_opt:      optimal first-step control input
    %           y_fut:      predicted output sequence
    %           cost:       cost function value
    %           info:       optimization solution info

    
    % System dimensions
    nu  = dim.nu;
    nx  = dim.nx;
    ny  = dim.ny;
    np  = dim.np;
    T   = dim.T;
    
    % Extract system matrices
    A   = reshape(cell2mat(sys.A), dim.nx, dim.nx, []);
    B   = reshape(cell2mat(sys.B), dim.nx, dim.nu, []);
    C   = reshape(cell2mat(sys.C), dim.ny, dim.nx, []);
    
    % Generate A, B, C matrices for each value of future scheduling
    A_fut = zeros(nx, nx, T+1);
    B_fut = zeros(nx, nu, T+1);
    C_fut = zeros(ny, nx, T+1);
    
    for i = 1:T
        A_fut(:,:,i) = sum(pagemtimes(A, reshape(kron(p_fut(:,i)', eye(nx)), nx, nx, np)), 3);
        B_fut(:,:,i) = sum(pagemtimes(B, reshape(kron(p_fut(:,i)', eye(nu)), nu, nu, np)), 3);
        C_fut(:,:,i) = sum(pagemtimes(C, reshape(kron(p_fut(:,i)', eye(nx)), nx, nx, np)), 3);
    end
    
    % Optimization parameters
    Q   = optim.Q;
    R   = optim.R;
    
    % Signal values
    u_max   = bounds.u_max;
    u_min   = bounds.u_min;
    y_max   = bounds.y_max;
    y_min   = bounds.y_min;
    
    % Set optimizer options
    yalmip('clear');
    options = sdpsettings('solver', 'gurobi', 'verbose', 0);
                      
    % Optimization variables
    U = sdpvar(nu*T, 1);            % containing values from t to t+T-1
    X = sdpvar(nx*(T+1), 1);        % containing values from t to t+T
    Y = sdpvar(ny*T, 1);            % containing values from t+1 to t+T
    
    % Initial conditions
    X(1:nx) = x0;
    
    % System dynamics
    for i = 1:T
        X(i*nx+1:(i+1)*nx) = A_fut(:,:,i) * X((i-1)*nx+1:i*nx) + B_fut(:,:,i) * U((i-1)*nu+1:i*nu);
        Y((i-1)*ny+1:i*ny) = C_fut(:,:,i) * X(i*nx+1:(i+1)*nx);
    end
    
    % Output saturation constraints
    cons = [Y <= repmat(y_max, T, 1)] + ...
           [Y >= repmat(y_min, T, 1)]; 
              
	% Input saturation constraints
    cons = cons + [U <= repmat(u_max, T, 1)] + ...
                  [U >= repmat(u_min, T, 1)];
    
	% Cost function
    Jcost = 0;
    for i = 1:T
        Jcost = Jcost + U((i-1)*nu+1:i*nu)' * R * U((i-1)*nu+1:i*nu) ...
                      + (Y((i-1)*ny+1:i*ny) - y_ref(i,:))' * Q * (Y((i-1)*ny+1:i*ny) - y_ref(i,:));
    end
    
    % Optimization
    info    = optimize(cons, Jcost, options);
    y_fut   = value(Y)';
    u_opt   = value(U(1:nu));
    cost    = double(Jcost);   
    
end




