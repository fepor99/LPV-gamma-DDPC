function row = build_UF_row(dim, u_fut, p_fut, idx)
    %   Construct only one row of UF given the required row index. Note
    %   that this implementation works only for nu = 1. The main rationale
    %   of the function is to first select the control input corresponding
    %   to the required row index, then the components of the scheduling 
    %   variable from time t+T-1 to time t are computed one at a time
    %   starting from time t+T-1. This procedure leverages the fact that
    %   the contribution of the leftmost Kronecker product can be easily
    %   distinguished if the Kronecker product matrix is divided into np
    %   equal parts
    %
    %   Inputs
    %           dim:        dimensions and horizons
    %           u_fut:      future input
    %           p_fut:      future scheduling
    %           idx:        number of row to be constructed
    %
    %   Outputs
    %           row:        UF(idx, :) row of matrix UF

    
    % Get system dimensions and check assumptions
    nu  = dim.nu;
    np  = dim.np;
    T   = dim.T;
    N   = dim.N;

	% Initialize row
    row = ones(1, N);
    
    % Divide index by the number of inputs to determine which scheduling element to use while building the row
    component = mod(idx-1, nu) + 1;
	
	% Determine time index of control input (of the first column)
    switch_idx  = [1, nu * np.^(T:-1:2)] * triu(ones(T));
    time_idx    = sum(idx - switch_idx >= 0);
    idx         = idx - switch_idx(time_idx) + 1;
    
    % Determine the combination of scheduling components starting from the leftmost scheduling variable
    for i = 1:T+1-time_idx
        
        % Height of p(t+T-i) (x) ... (x) p(t+time_idx-1) (x) u(t+time_idx-1)
        kron_height = nu * np^(T - time_idx - i + 2);
        
        % Component of p(t+T-i) to consider to generate the requested row
        p_comp  = floor(np * (idx-1) / kron_height) + 1;
        
        % Update scheduling product and position index for p(t+T-i-1) (x) ... (x) p(t+time_idx-1)
        row = row .* p_fut(p_comp, T-i+1:T-i+N);
        idx = idx - (p_comp-1) * kron_height / np;
    end
    
    % Add multiplication by control input and divide by samples length
    row = row .* u_fut(component, time_idx:time_idx+N-1) / sqrt(N);
    
    % Check for correct construction
    assert(idx == component, 'Something went wrong in the construction of the product series')
end


