function reduced_LQ_codegen(dim, YF, u_past, y_past, u_fut, p_past, p_fut, reduction, force_gen)
    %   Generate C code for the select_rows function, only if it has been
    %   modified since the last code generation, or if forced with the 
    %   force_gen parameter
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
    %           force_gen:  force the code generation if different from 0
    
    
    % Get output file paths
    [output_folder, ~, ~]   = fileparts(mfilename('fullpath'));
    output_file             = [output_folder '/' 'reduced_LQ_mex.' mexext];
    output_c_code_folder    = [output_folder '/codegen/mex/reduced_LQ'];

    % Check if the code must be generated
    if (~isfile(output_file) || ~isfolder(output_c_code_folder))
        % Generate C code if if it has never been done before
        flag = 1;
        
    else
        % Get modification date of function and code
        locale = get(0, 'Language');
        
        fun_modification_date_str       = dir([output_folder '/' 'reduced_LQ.m']).date;
        fun_modification_date           = datetime(fun_modification_date_str, 'Locale', locale);
        codegen_modification_date_str   = dir(output_file).date;
        codegen_modification_date       = datetime(codegen_modification_date_str, 'Locale', locale);
        
        % Generate C code if the function file is newer than the generated code
        if (fun_modification_date >= codegen_modification_date)
            flag = 1;
        else
            flag = 0;
        end
    end
    
    % Generate code only if needed
    if (flag || force_gen)
        fprintf('Generating C code for reduced_LQ function...');
        
        % System dimensions for variable declaration
        nu = dim.nu;
        ny = dim.ny;
        np = dim.np;
        
        % Function arguments
        args = {dim, ...
                coder.typeof(YF, [inf, inf], 1), ...
                coder.typeof(u_past, [nu, inf], [0, 1]), ...
                coder.typeof(y_past, [ny, inf], [0, 1]), ...
                coder.typeof(u_fut, [nu, inf], [0, 1]), ...
                coder.typeof(p_past, [np, inf], [0, 1]), ...
                coder.typeof(p_fut, [np, inf], [0, 1]), ...
                reduction};     
            
        % Generate code
        codegen('reduced_LQ', ...
                '-args', args, ...
                '-o', output_file, ...
                '-d', output_c_code_folder);
        
        fprintf('\tDone!\n\n')
        
    else
        fprintf('Function reduced_LQ has not been modified, skipping code generation\n\n');
    end
    
end


