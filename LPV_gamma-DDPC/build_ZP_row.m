function row = build_ZP_row(dim, u_past, y_past, p_past, eta, p_fut, idx)
    %   Construct only one row of ZP given the required row index. Note
    %   that this implementation works only for nu = 1 and ny = 1. The main
    %   rationale of the function is to first compute the required
    %   components of the Kronecker products between future scheduling
    %   signals, to then follow the same procedure described in the
    %   function build_UF_row()
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           u_past:     past input
    %           y_past:     past output
    %           p_past:     past scheduling
    %           eta:        past scheduling (in the form of p (x) p)
    %           p_fut:      future scheduling
    %           idx:        number of row to be constructed
    %
    %   Outputs
    %           row:        ZP(idx, :) row of matrix ZP

    
    % Get system dimensions and check assumptions
    ny  = dim.ny;
    nu  = dim.nu;
    np  = dim.np;
    T   = dim.T;
    N   = dim.N;
    M   = dim.M;

	% Initialize row
    row  = ones(1, N);
	neta = np*(np+1)/2;
    
    % Determine the combination of future scheduling components starting from the leftmost scheduling variable
    for i = 1:T
        
        % Height of p(t+T-1) (x) ... (x) p(t) (x) [up(t); yp(t)]
        kron_height = np^(T - i + 1) * ...
                      (nu*(np*(np+1)/2) + ny*np) * sum(cumprod([1, np*(np+1)/2 * ones(1,M-1)]));
        
        % Component of p(t+T-i) to consider to generate the requested row
        eta_comp  = floor(np * (idx-1) / kron_height) + 1;
        
        % Update scheduling product and position index for p(t+T-i-1) (x) ... (x) p(t+time_idx-1)
        row = row .* p_fut(eta_comp, T-i+1:T-i+N);
        idx = idx - (eta_comp-1) * kron_height / np;
    end
           
	% Determine whether to use input or output data and its corresponding time index
    if idx < nu*polyval([ones(1,M), 0], neta) + 1

    	% Use control input
        switch_idx  = [1, nu * neta.^(M:-1:2)] * triu(ones(M));
        time_idx    = sum(idx - switch_idx >= 0);
        component   = mod(idx-1, nu) + 1;
        data        = u_past(component, time_idx:time_idx+N-1);
        idx         = idx - switch_idx(time_idx) + 1;

        % Determine the combination of scheduling components starting from the leftmost scheduling variable
        for i = 1:M+1-time_idx

            % Height of eta(t-i) (x) ... (x) eta(t-M+time_idx-1) (x) u(t-M+time_idx-1)
            kron_height = nu*neta^(M - time_idx - i + 2);

            % Component of eta(t-i) to consider to generate the requested row
            eta_comp  = floor(neta * (idx-1) / kron_height) + 1;

            % Update scheduling product and position index for eta(t-i-1) (x) ... (x) eta(t-M+time_idx-1)
            row = row .* eta(eta_comp, M-i+1:M-i+N);
            idx = idx - (eta_comp-1) * kron_height / neta;
        end
    else

    	% Use output
        switch_idx  = [1, ny * np * neta.^(M-1:-1:1)] * triu(ones(M));
        idx         = idx - nu*polyval([ones(1,M), 0], neta);
        time_idx    = sum(idx - switch_idx >= 0);
        component   = mod(idx-1, ny) + 1;
        data        = y_past(component, time_idx:time_idx+N-1);   
        idx         = idx - switch_idx(time_idx) + 1;

        % Determine the combination of scheduling components starting from the leftmost scheduling variable
        for i = 1:M-time_idx

        	% Height of eta(t-i) (x) ... (x) eta(t-M+time_idx) (x) p(t-M+time_idx-1) (x) y(t-M+time_idx-1)
            kron_height = ny*np*neta^(M - time_idx - i + 1);

            % Component of eta(t-i) to consider to generate the requested row
            eta_comp  = floor(neta * (idx-1) / kron_height) + 1;
                
            % Update scheduling product and position index for eta(t-i-1) (x) ... (x) eta(t-M+time_idx)
            row = row .* eta(eta_comp, M-i+1:M-i+N);
            idx = idx - (eta_comp-1) * kron_height / neta;
        end

        % Add last Kronecker product with p_past
        kron_height = np*ny;
        p_comp      = floor(np * (idx-1) / kron_height) + 1;
        row         = row .* p_past(p_comp, time_idx:time_idx+N-1);
        idx         = idx - (p_comp-1) * kron_height / np;
    end

	% Add multiplication by control input or output and divide by samples length
    row = row .* data / sqrt(N);

    % Check for correct construction
	assert(idx == component, 'Something went wrong in the construction of the product series')
    
end


