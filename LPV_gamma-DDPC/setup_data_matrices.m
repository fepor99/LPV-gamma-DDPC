function [IZ, IU, LQ, res] = setup_data_matrices(plt, dim, reduction, force)
    %   Generates the library of past/future collected data and computes
    %   the LQ decomposition of the reduced [ZP; UF; YF] matrix
    %
    %   Inputs
    %           plt:        struct containing data collection signals
    %           dim:        struct containing system horizons
    %           reduction:  struct containing all the information for the
    %                       reduction of the ZP and UF matrices
    %           force:      force the C code generation if different from 0           
    % 
    %   Outputs
    %           IZ:         ordered indices of the selected rows of ZP
    %           IU:         ordered indices of the selected rows of UF
    %           LQ:         struct containing LQ decomposition matrices
    %           res:        row selection residuals
    

    % Extract horizons
    T = dim.T;
    M = dim.M;
    N = dim.N;  
    
    % Extract signals
    u = plt.collection.u;
    y = plt.collection.y;
    p = plt.collection.p;
    
    % Build past data
    u_past = u(:, 1:M+N-1);
    y_past = y(:, 1:M+N-1);
    p_past = [ones(1, width(y_past)); p(:, 1:M+N-1)];

    % Build future data
    u_fut  = u(:, M+1:end);
    y_fut  = y(:, M+1:end);
    p_fut  = [ones(1, width(y_fut)); p(:, M+1:end)];

    % Hankel matrices
    temp   = hankel(1:T, T:T+N-1);
    YF     = reshape(y_fut(:, temp), height(y_fut)*T, []) / sqrt(N);

    % Compute LQ factorization with reduced size versions of ZP and UF
    try 
        reduced_LQ_codegen(dim, YF, u_past, y_past, u_fut, p_past, p_fut, reduction, force);
        [H, IZ, IU, res] = reduced_LQ_mex(dim, YF, u_past, y_past, u_fut, p_past, p_fut, reduction);
    catch
        [H, IZ, IU, res] = reduced_LQ(dim, YF, u_past, y_past, u_fut, p_past, p_fut, reduction);
    end
    L = tril(H(:, 1:height(H)));
    Q = House_apply(eye(N), H);     % Q matrix of LQ decomposition is not used
    Q = Q(:, 1:height(H)).';
    
    % Add initial norm to residuals
    res = [norm(YF, 'fro')^2, res];

    % Splice L
    rows = reduction.rows;
    
    L11         = L(1:rows.z, 1:rows.z);
    LQ.L21      = L(rows.z+1:rows.z+rows.u, 1:rows.z);
    LQ.L22      = L(rows.z+1:rows.z+rows.u, rows.z+1:rows.z+rows.u);
    LQ.L31      = L(rows.z+rows.u+1:end, 1:rows.z);
    LQ.L3       = L(rows.z+rows.u+1:end, rows.z+1:end);
    LQ.L11_inv  = L11 \ eye(rows.z);
    fprintf('Computed the reduced LQ decomposition!\n\n')

end


