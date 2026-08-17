function [H, IZ, IU, res_old] = reduced_LQ(dim, YF, u_past, y_past, u_fut, p_past, p_fut, reduction)
    %   Computes the reduced version of the LQ decomposition. The
    %   algorithm to select the most relevant rows follows as similar 
    %   procedure to what is described in V. Verdult, M. Verhaegen,
    %   "Subspace identification of multivariable linear parameter-varying 
    %   systems", 2001
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           YF:         Hankel matrix of future outputs
    %           u_past:     past input
    %           y_past:     past output
    %           u_fut:      future input
    %           p_past:     past scheduling
    %           p_fut:      future scheduling
    %           reduction:  struct containing all the information for the
    %                       reduction of the ZP and UF matrices
    %
    %   Outputs
    %           H:          matrix containing both the triangular matrix
    %                       and the Householder vectors of the LQ decomposition
    %           IZ:         ordered indices of the selected rows of ZP
    %           IU:         ordered indices of the selected rows of UF
    %           res_old:    residual after selecting each row
    
    
    %% Algorithm initialization   
    % System and Hankel matrix dimensions
    ny   = dim.ny;
    nu   = dim.nu;
    np   = dim.np;
    T    = dim.T;
    N    = dim.N;
    M    = dim.M;
    rows = reduction.rows;
    
    % Get indices of rows without too many multiplications
    [IZ_no_p, IU_no_p, IZ_with_p, IU_with_p] = reduce_indices(dim, reduction);
    
    fprintf('Original row size of ZP: %.3e\n', cast(rows.z_full, "double"));
    fprintf('Size neglecting rows with more than %d multiplications with scheduling in ZP: %d\n\n', ...
                     cast(reduction.mult_max_ZP, "int16"), cast(length(IZ_no_p) + length(IZ_with_p), "int32"));    
    fprintf('Original row size of UF: %.3e\n', cast(rows.u_full, "double"));
    fprintf('Size neglecting rows with more than %d multiplications with scheduling in UF: %d\n', ...
                     cast(reduction.mult_max_UF, "int16"), cast(length(IU_no_p) + length(IU_with_p), "int32")); 
    
    % Initialize matrices
    IZ  = zeros(rows.z, 1);                     % Selected rows of ZP
    IU  = zeros(rows.u, 1);                     % Selected rows of UF
    H   = zeros(rows.z + rows.u + ny*T, N);     % Householder vector storage
    
    % Initialize residual vector
    res_old = zeros(1, rows.z + rows.u);
    
    
    %% Select rows of ZP without products with scheduling
    fprintf('\nStarting selection of %d rows for ZP\n\n', cast(rows.z, "int32"));
    fprintf('+-------------------+----------------+----------------------------------------------------------+\n');
    fprintf('|%-3sSelected rows%-3s|%-4sResidual%-4s|%-17sEstimated remaining time%-17s|\n','','','','','','');
    fprintf('+-------------------+----------------+-----------------+------------------+---------------------+\n');
    
    % Scheduling Kronecker product with itself
    eta = zeros(np*(np+1)/2, M+N-1);
    for i = 1:M+N-1
        prod = p_past(:,i) * p_past(:,i).';
        eta(:,i) = prod(tril(true(size(prod))));
    end

    for j = 1:length(IZ_no_p)
        
        % Construct j-th row and apply Householder transformation
        Z_best = build_ZP_row(dim, u_past, y_past, p_past, eta, p_fut, IZ_no_p(j));
        Z_best = House_apply(Z_best, H(1:j-1, :));

        % Compute and update residual
        [v, b]  = House_vec(Z_best(:, j:end));
        if j == 1
        	res_val = norm(YF(:, j:end) * ([zeros(1,N-j); eye(N-1,N-j)] - b*(v.'*v(2:end-j+1))), 'fro')^2;
        else
            res_val = res_old(j-1) - norm(YF(:, j:end) * (eye(N-j+1,1) - b*(v.'*v(1))))^2;
        end
        
        % Update residual and save index
        res_old(j)  = res_val;
        IZ(j)       = IZ_no_p(j);
        
        % Compute Householder transformation
        [v, b] = House_vec(Z_best(:, j:end));
        H(j, :) = Z_best * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        if j < N
            H(j, j+1:end) = v(2:end);
        end
        
        % Apply transformation
        YF = YF * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        
        % Print status
        fprintf('|%8s / %-8s|', sprintf('%d', cast(j, "int32")), sprintf('%d', cast(rows.z, "int32")));
        fprintf('%-4s%.2e%-4s|', '', cast(res_old(j), "single"), '');
        fprintf('%6s- hours%-4s|', '','');
        fprintf('%5s- minutes%-4s|', '','');
        fprintf('%8s- seconds%-4s|\n', '','');
    end
        
    
    %% Determine number of rows to eliminate at each iteration
    drop_rate = 0;      % initialize for C code generation
    switch reduction.residual_method_ZP
        
        % Linear decrease in rows at each iteration
        case "lin"
            row_elim = floor((length(IZ_with_p) - rows.z + (nu+ny)*M) / (rows.z - (nu+ny)*M)) * ones(rows.z - (nu+ny)*M, 1);   
            fprintf('\nAt each iteration %d rows are dropped\n\n', cast(row_elim(1), "int16"));
    
        % Exponential decrease 
        case "exp"
            row_elim = zeros(rows.z - (nu+ny)*M, 1);
            if length(IZ_with_p) - rows.z + (nu+ny)*M > 0
                drop_rate = 1 - (1 / (length(IZ_with_p) - rows.z + (nu+ny)*M))^(1 / (rows.z - (nu+ny)*M));

                extra_rows = length(IZ_with_p) - rows.z + (nu+ny)*M;
                for j = 1:rows.z - (nu+ny)*M
                    row_elim(j, :) = floor(drop_rate*extra_rows);
                    extra_rows = extra_rows - row_elim(j, :);
                end
                fprintf('\nAt each iteration %.3f%% of the remaining rows are dropped\n\n', cast(100*drop_rate, "single"));
            else
                fprintf('\n0%% of remaining rows are dropped at each iteration\n\n');
            end
                
        % Not a valid method is chosen
        otherwise
            row_elim = zeros(rows.z - (nu+ny)*M, 1);
            fprintf('\nNo large residual rows are dropped at each iteration\n\n');
    end
    
    
    %% Select most relevant rows of ZP with products with scheduling
    for j = length(IZ_no_p)+1:rows.z
        
        tic 
        % Select j-th row (1)
        res_val = Inf;
        res_idx = 0;
        Z_best  = zeros(1, N);
        elim_res = zeros(row_elim(j-(nu+ny)*M), 1);
        elim_idx = zeros(row_elim(j-(nu+ny)*M), 1);
            
        for k = 1:length(IZ_with_p)
        	% Construct k-th row and apply Householder transformation
            Z = build_ZP_row(dim, u_past, y_past, p_past, eta, p_fut, IZ_with_p(k));
            Z = House_apply(Z, H(1:j-1, :));

            % Compute and update residual
            [v, b] = House_vec(Z(:, j:end));
            res = res_old(j-1) - norm(YF(:, j:end) * (eye(N-j+1,1) - b*(v.'*v(1))))^2;

            % Select current most relevant row
            if res < res_val
            	Z_best  = Z;                % most relevant row
                res_val = res;              % residual value
                res_idx = k;                % residual index
            end
            
            % Populate worst residual vector
            if ~isempty(elim_res) && res > elim_res(1)
                elim_res(1) = res;
                elim_idx(1) = k;
                [elim_res, order] = sort(elim_res);
                elim_idx = elim_idx(order);
            end
        end
        
        % Update residual, save index and drop chosen index
        res_old(j)          = res_val;
        IZ(j)               = IZ_with_p(res_idx);
        IZ_with_p(res_idx)  = NaN;
        IZ_with_p(elim_idx) = NaN;
        
        IZ_with_p(isnan(IZ_with_p)) = [];
        
        % Compute Householder transformation
        [v, b] = House_vec(Z_best(:, j:end));
        H(j, :) = Z_best * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        if j < N
            H(j, j+1:end) = v(2:end);
        end
        
        % Apply rotation
        YF = YF * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        
        % Print status
        if strcmp(reduction.residual_method_ZP, "exp")
            time_row = toc * polyval([ones(1, rows.z - j), 0], 1-drop_rate);
        else
            time_row = toc * (rows.z - j);
        end
        fprintf('|%8s / %-8s|', sprintf('%d', cast(j, "int32")), sprintf('%d', cast(rows.z, "int32")));
        fprintf('%-4s%.2e%-4s|', '', cast(res_old(j), "single"), '');
        fprintf('%-4s%3d hours%-4s|', '', cast(floor(time_row/3600), "int16"), '');
        fprintf('%-4s%2d minutes%-4s|', '', cast(floor(rem(time_row, 3600)/60), "int8"), '');
        fprintf('%-4s%5.2f seconds%-4s|\n', '', rem(rem(time_row, 3600), 60), '');
    end  
    fprintf('+-------------------+----------------+-----------------+------------------+---------------------+\n');
    
    
    %% Select rows of UF without products with scheduling
    fprintf('\nStarting selection of %d rows for UF\n\n', cast(rows.u, "int32"));
    fprintf('+-------------------+----------------+----------------------------------------------------------+\n');
    fprintf('|%-3sSelected rows%-3s|%-4sResidual%-4s|%-17sEstimated remaining time%-17s|\n','','','','','','');
    fprintf('+-------------------+----------------+-----------------+------------------+---------------------+\n');
    
    for j = rows.z+1:rows.z+length(IU_no_p)
        
        % Construct j-th row and apply Householder transformation
        U_best = build_UF_row(dim, u_fut, p_fut, IU_no_p(j - rows.z));
        U_best = House_apply(U_best, H(1:j-1, :));

        % Compute and update residual
        [v, b]  = House_vec(U_best(:, j:end));
        res_val = res_old(j-1) - norm(YF(:, j:end) * (eye(N-j+1,1) - b*(v.'*v(1))))^2;
        
        % Update residual and save index
        res_old(j)   = res_val;
        IU(j-rows.z) = IU_no_p(j - rows.z);
        
        % Compute Householder transformation
        [v, b] = House_vec(U_best(:, j:end));
        H(j, :) = U_best * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        if j < N
            H(j, j+1:end) = v(2:end);
        end
        
        % Apply rotation
        YF = YF * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        
        % Print status
        fprintf('|%8s / %-8s|', sprintf('%d', cast(j-rows.z, "int32")), sprintf('%d', cast(rows.u, "int32")));
        fprintf('%-4s%.2e%-4s|', '', cast(res_old(j), "single"), '');
        fprintf('%6s- hours%-4s|', '','');
        fprintf('%5s- minutes%-4s|', '','');
        fprintf('%8s- seconds%-4s|\n', '','');
    end
    
    
    %% Determine number of rows to eliminate at each iteration
    switch reduction.residual_method_UF
        
        % Linear decrease in rows at each iteration
        case "lin"
            row_elim = floor((length(IU_with_p) - rows.u + nu*T) / (rows.u - nu*T)) * ones(rows.u - nu*T, 1);   
            fprintf('\nAt each iteration %d rows are dropped\n\n', cast(row_elim(1), "int16"));
    
        % Exponential decrease 
        case "exp"
            row_elim = zeros(rows.u - nu*T, 1);
            if length(IU_with_p) - rows.u + nu*T > 0
                drop_rate = 1 - (1 / (length(IU_with_p) - rows.u + nu*T))^(1 / (rows.u - nu*T));

                extra_rows = length(IU_with_p) - rows.u + nu*T;
                for j = 1:rows.u - nu*T
                    row_elim(j, :) = floor(drop_rate*extra_rows);
                    extra_rows = extra_rows - row_elim(j, :);
                end
                fprintf('\nAt each iteration %.3f%% of the remaining rows are dropped\n\n', cast(100*drop_rate, "single"));
            else
                fprintf('\n0%% of remaining rows are dropped at each iteration\n\n');
            end
                
        % Not a valid method is chosen
        otherwise
            row_elim = zeros(rows.u - nu*T, 1);
            fprintf('\nNo large residual rows are dropped at each iteration\n\n');
    end
    
    
    
    %% Select most relevant rows of UF with products with scheduling
    for j = rows.z+length(IU_no_p)+1:rows.z+rows.u
        
        tic 
        % Select j-th row
        res_val = Inf;
        res_idx = 0;
        U_best  = zeros(1, N);
        elim_res = zeros(row_elim(j-nu*T-rows.z), 1);
        elim_idx = zeros(row_elim(j-nu*T-rows.z), 1);
            
        for k = 1:length(IU_with_p)            
        	% Construct k-th row and apply Householder transformation
            U = build_UF_row(dim, u_fut, p_fut, IU_with_p(k));
            U = House_apply(U, H(1:j-1, :));

            % Compute and update residual
            [v, b] = House_vec(U(:, j:end));
            res = res_old(j-1) - norm(YF(:, j:end) * (eye(N-j+1,1) - b*(v.'*v(1))))^2;

            % Select current most relevant row
            if res < res_val
            	U_best  = U;                % most relevant row
                res_val = res;              % residual value
                res_idx = k;                % residual index
            end
            
            % Populate worst residual vector
            if ~isempty(elim_res) && res > elim_res(1)
                elim_res(1) = res;
                elim_idx(1) = k;
                [elim_res, order] = sort(elim_res);
                elim_idx = elim_idx(order);
            end
        end
        
        % Update residual, save index and drop chosen index
        res_old(j)          = res_val;
        IU(j - rows.z)      = IU_with_p(res_idx);
        IU_with_p(res_idx)  = NaN;
        IU_with_p(elim_idx) = NaN;
        
        IU_with_p(isnan(IU_with_p)) = [];
        
        % Compute Householder transformation
        [v, b] = House_vec(U_best(:, j:end));
        H(j, :) = U_best * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        if j < N
            H(j, j+1:end) = v(2:end);
        end
        
        % Apply transformation
        YF = YF * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        
        % Print status
        if strcmp(reduction.residual_method_ZP, "exp")
            time_row = toc * polyval([ones(1, rows.z + rows.u - j), 0], 1-drop_rate);
        else
            time_row = toc * (rows.z + rows.u - j);
        end
        fprintf('|%8s / %-8s|', sprintf('%d', cast(j-rows.z, "int32")), sprintf('%d', cast(rows.u, "int32")));
        fprintf('%-4s%.2e%-4s|', '', cast(res_old(j), "single"), '');
        fprintf('%-4s%3d hours%-4s|', '', cast(floor(time_row/3600), "int16"), '');
        fprintf('%-4s%2d minutes%-4s|', '', cast(floor(rem(time_row, 3600)/60), "int8"), '');
        fprintf('%-4s%5.2f seconds%-4s|\n', '', rem(rem(time_row, 3600), 60), '');
    end  
    fprintf('+-------------------+----------------+-----------------+------------------+---------------------+\n\n');
    
    
    %% Compute LQ factorization for YF
    for j = rows.z+rows.u+1:rows.z+rows.u+ny*T
        
        % Compute Householder transformation
        [v, b] = House_vec(YF(j-rows.z-rows.u, j:end));
        YF = YF * blkdiag(eye(j-1), (eye(N-j+1) - b*(v.'*v)));
        
        % Save LQ factorization with Householder vectors
        H(j, :) = YF(j-rows.z-rows.u, :);
        if j < N
            H(j, j+1:end) = v(2:end);
        end        
    end

end



