function z_init = build_init_reduced(dim, u_past, y_past, p_past, p_fut, IZ)
    %   Generate the initial condition vector z_init for the system given 
    %   the system dimensions, the past and future horizon, past 
    %   input/output data, estimated future scheduling variable and a set
    %   of relevant rows to build
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           u_past:     past input
    %           y_past:     past output
    %           p_past:     past scheduling
    %           p_fut:      future scheduling
    %           IZ:         array of selected row indices to construct
    %
    %   Outputs
    %           z_init:     reduced initial conditions of the system


    % Get system and Hankel matrix dimensions
    ny  = dim.ny;
    nu  = dim.nu;
    np  = dim.np;
    M   = dim.M;
    T   = dim.T;
    
    % Scheduling Kronecker product with itself
    neta = np*(np+1)/2;
    eta = zeros(np*(np+1)/2, width(p_past));
    for i = 1:width(p_past)
        prod = p_past(:,i) * p_past(:,i).';
        eta(:,i) = prod(tril(true(size(prod))));
    end

    % Initialize function output
    z_init = ones(length(IZ), 1);
    
    % Generate one row at a time
    for j = 1:length(IZ)
        
        % Get row index
        idx = IZ(j); 
        
        % Determine the combination of future scheduling components starting from the leftmost scheduling variable
        for i = 1:T

            % Height of p(t+T-i) (x) ... (x) p(t) (x) [up(t); yp(t)]
            kron_height = np^(T - i + 1) * ...
                          (nu*(np*(np+1)/2) + ny*np) * polyval(ones(1,M), np*(np+1)/2);

            % Component of p(t+T-i) to consider to generate the requested row
            eta_comp  = floor(np * (idx-1) / kron_height) + 1;

            % Update scheduling product and position index for p(t+T-i-1) (x) ... (x) p(t+time_idx-1)
            z_init(j) = z_init(j) * p_fut(eta_comp, T-i+1);
            idx = idx - (eta_comp-1) * kron_height / np;
        end

        % Determine whether to use input or output data and its corresponding time index
        if idx < nu*polyval([ones(1,M), 0], neta) + 1

            % Use control input
            switch_idx  = [1, nu * neta.^(M:-1:2)] * triu(ones(M));
            time_idx    = sum(idx - switch_idx >= 0);
            component   = mod(idx-1, nu) + 1;
            data        = u_past(component, time_idx);
            idx         = idx - switch_idx(time_idx) + 1;

            % Determine the combination of past scheduling components starting from the leftmost scheduling variable
            for i = 1:M+1-time_idx

                % Height of eta(t-i) (x) ... (x) eta(t-M+time_idx-1) (x) u(t-M+time_idx-1)
                kron_height = nu*neta^(M - time_idx - i + 2);

                % Component of eta(t-i) to consider to generate the requested row
                eta_comp  = floor(neta * (idx-1) / kron_height) + 1;

                % Update scheduling product and position index for eta(t-i-1) (x) ... (x) eta(t-M+time_idx-1)
                z_init(j) = z_init(j) * eta(eta_comp, M-i+1);
                idx = idx - (eta_comp-1) * kron_height / neta;
            end
        else

            % Use output
            switch_idx  = [1, ny * np * neta.^(M-1:-1:1)] * triu(ones(M));
            idx         = idx - nu*polyval([ones(1,M), 0], neta);
            time_idx    = sum(idx - switch_idx >= 0);
            component   = mod(idx-1, ny) + 1;
            data        = y_past(component, time_idx);   
            idx         = idx - switch_idx(time_idx) + 1;

            % Determine the combination of past scheduling components starting from the leftmost scheduling variable
            for i = 1:M-time_idx

                % Height of eta(t-i) (x) ... (x) eta(t-M+time_idx) (x) p(t-M+time_idx-1) (x) y(t-M+time_idx-1)
                kron_height = ny*np*neta^(M - time_idx - i + 1);

                % Component of eta(t-i) to consider to generate the requested row
                eta_comp  = floor(neta * (idx-1) / kron_height) + 1;

                % Update scheduling product and position index for eta(t-i-1) (x) ... (x) eta(t-M+time_idx)
                z_init(j) = z_init(j) * eta(eta_comp, M-i+1);
                idx = idx - (eta_comp-1) * kron_height / neta;
            end

            % Add last Kronecker product with p_past
            kron_height = np*ny;
            p_comp      = floor(np * (idx-1) / kron_height) + 1;
            z_init(j)	= z_init(j) * p_past(p_comp, time_idx);
            idx         = idx - (p_comp-1) * kron_height / np;
        end

        % Add multiplication by control input or output
        z_init(j) = z_init(j) * data;
        
        % Check for correct construction
        assert(idx == component, 'Something went wrong in the construction of the product series')
    end
    
end



