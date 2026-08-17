function [ZP_occ, UF_occ, ETA_occ] = get_mult_frequency(dim)
	%   This function computes the occurrence for every number of
	%   multiplications with the scheduling signal. For example, given
    %   p = [1, p1, p2]', if the algorithm has to compute the occurrence 
    %   for the vector 
    %
    %      kron(p, p) = [1, p1, p2, p1, p1^2, p1*p2, p2, p1*p2, p2^2]'
    %
    %   then it will return the vector [4, 4, 1], i.e. 4 double products, 4
    %   single products, 1 element without products with scheduling 
    %   signals. This function does not give information about the position
    %   of the products, just their multiplicity.
    % 
    %   To better understand the structure of the following code, remember
    %   that
    %                                 | etaM (x) ... (x) eta1 (x) u1 |
    %                                 | etaM (x) ... (x) eta2 (x) u2 |
    %                                 |              ...             |
    %      ZP = pT (x) ... (x) p1 (x) | etaM (x) uM                  |
    %                                 | etaM (x) ... (x) p1 (x) y1   |
    %                                 |              ...             |
    %                                 | etaM (x) ... (x) p2 (x) y2   |
    %                                 | pM (x) yM                    |
    %
    %           | pT (x) ... (x) p1 (x) u1 |
    %      UF = | pT (x) ... (x) p2 (x) u2 |
    %           |            ...           |
    %           | pT (x) uM                |
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           
    %   Outputs
    %           ZP_occ:     number of occurences in ZP
    %           UF_ovv:     number of occurences in UF
    %           ETA_occ:    number of occurences in eta^M
    

    % Extract scheduling dimension, past and future horizons
    np  = dim.np;
    nu  = dim.nu;
    ny  = dim.ny;
    M   = dim.M;
    T   = dim.T;

    % eta contains 1 element with no multiplications with the scheduling, 
    % np-1 with 1 multiplications, and np*(np-1)/2 with 2 multiplications
	eta_mult = [np * (np - 1)/2, np-1, 1];
    
    % Compute occurrence in eta^(M-1) and of the stacked matrix 
    % eta_mat = [eta^(M-1); eta^(M-2); ...; eta^2; eta; 1]
    eta_mat_occ = [zeros(1, 2*(M-1)), 1];       % stacked matrix
    eta_exp_occ = [zeros(1, 2*(M-1)), 1];       % eta^(M-1)
    for i = 1:M-1
        eta_exp_occ(2*(M-i)-1:end) = conv(eta_exp_occ(2*(M-i)+1:end), eta_mult);
        eta_mat_occ(2*(M-i)-1:end) = eta_mat_occ(2*(M-i)-1:end) + eta_exp_occ(2*(M-i)-1:end);
    end
    
    % Compute occurrence in eta^M
    ETA_occ = conv(eta_exp_occ, eta_mult);
    
    % Combine occurrence of past inputs and outputs
    rows_u = ETA_occ + [zeros(1,2), eta_mat_occ];
    rows_y = [0, conv([np-1, 1], eta_mat_occ)];
    eta_mat_occ = rows_u*nu + rows_y*ny;
    eta_mat_occ(end) = eta_mat_occ(end) - nu;
    
    % Add occurrences of the future scheduling matrix p_mat = [pT (x) ... (x) p1]
    % p contains np-1 elements with 1 multiplication and 1 element with no multiplications
    p_mat = [zeros(1, T), 1];
    for i = 1:T
        p_mat(T-i+1:end) = conv(p_mat(T-i+2:end), [np-1, 1]);
    end
    ZP_occ = conv(p_mat, eta_mat_occ);
    
    % Compute occurrence for UF in a similar way to eta_mat_occ
	UF_occ = zeros(1, T+1);
    eta_exp_occ = [zeros(1, T), 1];
    for i = 1:T
        eta_exp_occ(T-i+1:end) = conv(eta_exp_occ(T-i+2:end), [np-1, 1]);
        UF_occ(T-i+1:end) = UF_occ(T-i+1:end) + eta_exp_occ(T-i+1:end);
    end  
    
    % Account for input dimensions
    UF_occ = UF_occ * dim.nu;
    
end


