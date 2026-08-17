function check_assumptions(dim, reduction, optim, bounds)
    %   Checks the all the assumptions to run the algorithm to build
    %   reduced versions of the ZP and UF matrices. Warnings are given to
    %   inform the user when the matrix reduction algorithm might be
    %   useless
    % 
    %   Inputs
    %           dim:        dimensions and horizons
    %           reduction:  struct containing all the information for the
    %                       reduction of the ZP and UF matrices
    %           optim:      struct containing optimization parameters
    %           bounds:     struct containing saturation bounds


	% Get system dimensions
    ny  = dim.ny;
    nu  = dim.nu;
    T   = dim.T;
    N   = dim.N;
    M   = dim.M;
    
    % Extract rows struct
    rows = reduction.rows;
    
    
    %% System dimensions assumptions
	assert(nu == 1, ['This implementation has been design to work specifically for ', ...)
                     'nu = 1. If you want to apply this script to other scenarios ', ...
                     'you probably have to modify part of the original code'])
                 
	assert(ny == 1, ['This implementation has been design to work specifically for ', ...)
                     'ny = 1. If you want to apply this script to other scenarios ', ...
                     'you probably have to modify part of the original code'])
                 
	assert(M > 0 && T >= 0, "Please set both T and M to be greater than zero")
                 
                 
	%% Matrix reduction assumptions
    % Rows dimension conditions
    assert(rows.z_full >= rows.z, 'The chosen number of rows for ZP is larger than the full height of ZP')
    assert(rows.u_full >= rows.u, 'The chosen number of rows for UF is larger than the full height of UF')
    
    assert(rows.z + rows.u + ny*T <= N, ['The chosen number of rows is too large. The reduced ' ...
                                         'matrix [ZP; UF; YF] would still be taller than larger'])
                                     
    assert(rows.z >= 2*M, ['The chosen number of rows for ZP is too small to include ' ...
                           'in the reduced matrix all the I/O data not multiplied by ' ...
                           'the scheduling variable. The minimum admissible number ' ...
                           'of row is ', num2str(2*M)])
                       
    assert(rows.u >= T, ['The chosen number of rows for UF is too small to include ' ...
                         'in the reduced matrix all the input data not multiplied by ' ...
                         'the scheduling variable. The minimum admissible number ' ...
                         'of row is ', num2str(T)])
    
    if N >= rows.z_full + rows.u_full + ny*T
        warning(['The value of N is bigger than the height of the full matrix ' ...
                 '[ZP; UF; YF]. If not for numerical issues, the LQ factorization ' ...
                 'can be computed with the full data matrices'])
    end
    
    if (rows.z_full == rows.z) && (rows.u_full == rows.u)
        warning(['No reduction is applied to ZP and UF'])
    end
    
    % Drop rows with many multiplications
    assert(reduction.mult_max_ZP > 0, 'reduction.mult_max_ZP must be greater than zero')
    assert(reduction.mult_max_UF > 0, 'reduction.mult_max_UF must be greater than zero')
    
    % Check if after scheduling multiplication reduction there are enough rows to select
    [IZ_no_p, IU_no_p, IZ_with_p, IU_with_p] = reduce_indices(dim, reduction);
    assert(length(IZ_no_p) + length(IZ_with_p) >= rows.z, ['Too many rows in ZP should be removed ', ...
                'due to the frequency of scheduling signal multiplications. The number of remaining ', ...
                'rows would be less than the desired number of rows. Please, either reduce ', ...
                'reduction.rows.z at least down to ', num2str(length(IZ_no_p) + length(IZ_with_p)), ...
                ' or increase reduction.mult_max_ZP'])
            
    assert(length(IU_no_p) + length(IU_with_p) >= rows.u, ['Too many rows in UF should be removed ', ...
                'due to the frequency of scheduling signal multiplications. The number of remaining ', ...
                'rows would be less than the desired number of rows. Please, either reduce ', ...
                'reduction.rows.u at least down to ', num2str(length(IU_no_p) + length(IU_with_p)), ...
                ' or increase reduction.mult_max_UF'])
    
    % Residual drop method
    assert(strlength(reduction.residual_method_ZP) == 3, ['Please, select a residual ', ...
                'drop method name for ZP with 3 letters, the .mex file is compiled assuming ', ...
                'this field has length 3'])
            
    assert(strlength(reduction.residual_method_UF) == 3, ['Please, select a residual ', ...
                'drop method name for UF with 3 letters, the .mex file is compiled assuming ', ...
                'this field has length 3'])
            
    if ~any(strcmp(reduction.residual_method_ZP, ["lin", "exp"]))
        warning(['No valid method has been chosen for the residual drop method for ZP. ', ...
                 'No rows will be dropped from ZP at each iteration due to their large residuals'])
    end
    
    if ~any(strcmp(reduction.residual_method_UF, ["lin", "exp"]))
        warning(['No valid method has been chosen for the residual drop method for UF. ', ...
                 'No rows will be dropped from UF at each iteration due to their large residuals'])
    end
    
    
    %% Optimization parameters assumptions
    assert(optim.nsim > 0, "Set the control simulation length to a value greater than zero")
    
    assert(all(size(optim.Q) == [ny, ny]) && all(size(optim.R) == [nu, nu]) && length(optim.beta2) == 1 ...
           && length(optim.beta3) == 1, 'The weight parameters in the optimization cost function have incorrect size')

    assert(issymmetric(optim.Q) && all(eig(optim.Q) >= 0), ['The weight matrix Q in the optimization ' ...
                        'cost function is not positive semidefinite'])
                    
	assert(issymmetric(optim.R) && all(eig(optim.R) >= 0), ['The weight matrix R in the optimization ' ...
                        'cost function is not positive semidefinite'])
    
    assert(optim.beta2 >= 0, 'The normalization parameter beta2 in the optimization cost function is negative')
    assert(optim.beta3 >= 0, 'The normalization parameter beta3 in the optimization cost function is negative')

    assert(all(optim.y_ref <= bounds.y_max) && all(optim.y_ref >= bounds.y_min), ['The reference ' ...
        'to be tracked during control simulation goes beyond the limits [', num2str(bounds.y_min), ...
        ', ', num2str(bounds.y_max), '] imposed on the output of the system during the optimization'])    
            
end


