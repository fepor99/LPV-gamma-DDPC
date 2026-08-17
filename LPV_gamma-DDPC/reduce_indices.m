function [IZ_no_p, IU_no_p, IZ_with_p, IU_with_p] = reduce_indices(dim, reduction)
    %   This function computes which indices of the matrices ZP and UF do 
    %   not contain more than a specified number of multiplications with
    %   scheduling signals. The number of maximum multiplications can be
    %   different between ZP and UF. Notice that each eta variable
    %   contains also double products between elements of the same
    %   scheduling variable, i.e. elements with already two multiplications
    % 
    %   Inputs
    %           dim:        struct of dimensions and horizons
    %           reduction:  struct containing all the information for the
    %                       reduction of the ZP and UF matrices
    %           
    %   Outputs
    %           IZ_no_p:    sorted set of indices of ZP with no multiplications with p
    %           IU_no_p:    sorted set of indices of ZP with no multiplications with p
    %           IZ_with_p:  unsorted set of indices of ZP with multiplications with p
    %           IU_with_p:  unsorted set of indices of UF with multiplications with p
    %
    %
    %   NOTES ABOUT AUXILIARY CELLS FOR COMPUTATIONS:
	%   indexEta is a cell of arrays. Each array is used to store the
    %   indices of a specified number of multiplications with the 
    %   scheduling signal. The first array contains the indices with no
    %   multiplications, the second array contains the indices with one
    %   multiplication, and so on.
    %   
    %   At the first iteration, indexEta contains only the indices of the 
    %   vector eta (if no indices correspond to a given number of 
    %   multiplications, then the array is empty). This is used to compute 
    %   the products for the row of u(t-1). At the second iteration it 
    %   contains only the indices of eta (x) eta (used for u(t-2)), at the 
    %   third of eta (x) eta (x) eta (for u(t-3)), and so on. indexP works 
    %   in a similar way but it is used to account for the components of 
    %   the future scheduling.
    %
    %   shiftEta and shiftP are used to shift the indices contained in
    %   indexEta and indexP, so that the indices of the multiplications
    %   refer to the correct rows of ZP and UF. For example, in ZP, when 
    %   indexEta accounts only for eta, the rows to build are the ones 
    %   multiplied by u(t-1), which are not the first rows of ZP. shiftEta
    %   then modifies the computed rows to reach the correct indices.

    
    % Extract system dimensions
    np      = dim.np;
    nu      = dim.nu;
    ny      = dim.ny;
    neta    = np * (np+1) / 2;
    T       = dim.T;
    M       = dim.M;
    
    % Compute the occurrence of each number of multiplications
	% (used only to initialize the size of the vectors)
    [ZP_occ, UF_occ, ETA_occ] = get_mult_frequency(dim);
    
	% Extract maximum multiplication limit
    limZP = reduction.mult_max_ZP;
    limUF = reduction.mult_max_UF;
    
    % Initialize cells for indices
    IZ_cell     = cell(min(2*M+T+1, limZP), 1);
    IU_cell     = cell(min(T+1, limUF), 1);
    Eta_cell    = cell(min(2*M+1, limZP), 1);
    
    % Instantiate memory for MATLAB coder
    IZ       = coder.nullcopy(IZ_cell);         % indices for ZP
    IU       = coder.nullcopy(IU_cell);         % indices for UF
    indexEta = coder.nullcopy(Eta_cell);        % auxiliary cell for computations
    indexP   = coder.nullcopy(IU_cell);         % auxiliary cell for computations
    
	for i = 1:length(IZ)
    	IZ{i} = zeros(ZP_occ(length(ZP_occ)-i+1), 1);
	end
    for i = 1:length(IU)
        IU{i} = zeros(UF_occ(length(UF_occ)-i+1), 1);
        indexP{i} = zeros(UF_occ(length(UF_occ)-i+1), 1);
    end
    for i = 1:length(indexEta)
        indexEta{i} = zeros(ETA_occ(length(ETA_occ)-i+1), 1);
    end
    
    % Initialize variables to count how many elements have been placed in
    % each cell of IZ, IU, indexEta, and indexP
    IZ_count    = zeros(length(IZ), 1);
    IU_count    = zeros(length(IU), 1);
    Eta_count   = zeros(length(indexEta), 1);
    P_count     = zeros(length(indexP), 1);
    
    
    %% Construct IZ   
    % Initialize indexEta
    % In the first iteration, indexEta accounts for eta, which contains 
    % elements with 0 multiplications with p (only 1), elements with 1
    % multiplications, and elements with 2 multiplications
    
    % 0 multiplications
    indexEta{1}(1) = 1;
    Eta_count(1) = 1;
    
    if limZP > 1
        % 1 multiplications
        indexEta{2}(1:np-1) = 2:np;
        Eta_count(2) = np-1;
        
        if limZP > 2
            % 2 multiplications
            indexEta{3}(1:neta-np) = np+1:neta;
            Eta_count(3) = neta-np;
        end
    end    
    
    % Shift indices related to u_past contained in indexEta and store them in IZ
    shiftEta = neta^2 * (neta^(M-1) - 1) / (neta - 1);
    for i = 1:min(3, limZP)
        
        % Get indices as if dim.nu=1
        indices = indexEta{i}(1:Eta_count(i)) + shiftEta;
        
        % Account for the input dimension and store in IZ
        % We are repeating elements and considering subsequent positions,
        % e.g., if indices=3 and nu=4, we want the elements [9, 10, 11, 12] 
        IZ{i}(1:nu*Eta_count(i)) = repelem(indices*nu, nu, 1) - repmat((nu-1:-1:0)', length(indices), 1);
        IZ_count(i) = nu*Eta_count(i);
    end
    
    % Shift indices related to y_past and store them in IZ
    shiftP = neta / (neta - 1) * ((neta^M - 1) * nu + (neta^(M-1) - 1) * np * ny); 
    
    % Account now for the output dimension
    IZ{1}(IZ_count(1)+1:IZ_count(1)+ny) = shiftP + (1:ny)';
    IZ_count(1) = IZ_count(1) + ny;
    if limZP > 1
        for i = 1:np-1
            IZ{2}(IZ_count(2)+1:IZ_count(2)+ny) = shiftP + (1:ny)' + i*ny;
            IZ_count(2) = IZ_count(2) + ny;
        end
    end
    
    % Loop over remaining number of products (increase elements in indexEta)
    for i = 1:M-1
        
        % Update shiftP and store new rows for y_past
        shiftP = shiftP - ny * np * neta^i;
        for q = 1:min(1 + 2*i, limZP)
            for j = 1:Eta_count(q)
                IZ{q}(IZ_count(q)+1:IZ_count(q)+ny) = (indexEta{q}(j) - 1)*np*ny + shiftP + (1:ny)';
                IZ_count(q) = IZ_count(q) + ny;
                for k = 1:np-1
                    if q < limZP
                        IZ{q+1}(IZ_count(q+1)+1:IZ_count(q+1)+ny) = (indexEta{q}(j) - 1)*np*ny + shiftP + (1:ny)' + k*ny;
                        IZ_count(q+1) = IZ_count(q+1) + ny;
                    end
                end
            end
        end
        
        % Increase the product number contained in indexEta
        for j = min(2*i + 3, limZP):-1:2
            for k = 1:np-1
                indexEta{j}(Eta_count(j)+1:Eta_count(j)+Eta_count(j-1)) = ...
                                k * neta^i + indexEta{j-1}(1:Eta_count(j-1));
                Eta_count(j) = Eta_count(j) + Eta_count(j-1);
            end
            
            if j > 2
                for k = 0:neta - np - 1
                    indexEta{j}(Eta_count(j)+1:Eta_count(j)+Eta_count(j-2)) = ...
                                    (k + np) * neta^i + indexEta{j-2}(1:Eta_count(j-2));
                    Eta_count(j) = Eta_count(j) + Eta_count(j-2);
                end
            end
        end
        
        % Update shiftEta and store new rows for u_past
        shiftEta = shiftEta - neta^(i+1);
        for j = 1:min(3 + 2*i, limZP)
            indices = indexEta{j}(1:Eta_count(j)) + shiftEta;
            IZ{j}(IZ_count(j)+1:IZ_count(j)+Eta_count(j)*nu) = ...
                            repelem(indices*nu, nu, 1) - repmat((nu-1:-1:0)', length(indices), 1);
            IZ_count(j) = IZ_count(j) + Eta_count(j)*nu;
        end
    end
    
    % Add effect of future scheduling matrix [pT (x) ... (x) p1]
    dimEtaP = (nu*neta + ny*np) * (neta^M - 1) / (neta - 1);        % height of ZP before multiplication by [pT (x) ... (x) p1]
    for i = 1:T
        for j = min(i + 2*M + 1, limZP):-1:2
            for k = 1:np-1
                IZ{j}(IZ_count(j)+1:IZ_count(j)+IZ_count(j-1)) = ...
                            k * np^(i-1) * dimEtaP + IZ{j-1}(1:IZ_count(j-1));
                IZ_count(j) = IZ_count(j) + IZ_count(j-1);
            end
        end
    end
    
    
    %% Construct IU
    % Initialize indexP
    indexP{1}(1) = 1;
    P_count(1) = 1;
    
    if limUF > 1
        indexP{2}(1:np-1) = 2:np;
        P_count(2) = np-1;
    end  
    
    % Shift indices related to u_fut contained in indexP and store them in UF
    shiftP = np^2 * (np^(T-1) - 1) / (np - 1);
    for i = 1:min(2, limUF)
        
        % Get indices as if dim.nu=1
        indices = indexP{i}(1:P_count(i)) + shiftP;
        
        % Account for the input dimension and store in IU
        % We are repeating elements and considering subsequent positions,
        % e.g., if indices=3 and nu=4, we want the elements [9, 10, 11, 12] 
        IU{i}(1:nu*P_count(i)) = repelem(indices*nu, nu, 1) - repmat((nu-1:-1:0)', length(indices), 1);
        IU_count(i) = nu*P_count(i);
    end
    
    % Loop over remaining number of products (increase elements in indexP)
    for i = 1:T-1
        
        % Increase the product number contain in indexP
        for j = min(i+2, limUF):-1:2
            for k = 1:np-1
                indexP{j}(P_count(j)+1:P_count(j)+P_count(j-1)) = ...
                                k * np^i + indexP{j-1}(1:P_count(j-1));
                P_count(j) = P_count(j) + P_count(j-1);
            end
        end
        
        % Update shiftP and store new rows for u_fut
        shiftP = shiftP - np^(i+1);
        for j = 1:min(i+2, limUF)
            indices = indexP{j}(1:P_count(j)) + shiftP;
            IU{j}(IU_count(j)+1:IU_count(j)+P_count(j)*nu) = ...
                            repelem(indices*nu, nu, 1) - repmat((nu-1:-1:0)', length(indices), 1);
            IU_count(j) = IU_count(j) + P_count(j)*nu;
        end
    end
    
    
    %% Extract data
    % Elements without products with scheduling signal
    IZ_no_p     = sort(IZ{1});
    IU_no_p     = sort(IU{1});
    
    % Elements with products with scheduling signal
    len_IZ = 0;
    len_IU = 0;
    for i = 2:length(IZ)
        len_IZ = len_IZ + length(IZ{i});
    end
    for i = 2:length(IU)
        len_IU = len_IU + length(IU{i});        
    end
    
    % Initialize matrices
    IZ_with_p   = zeros(len_IZ, 1);
    IU_with_p   = zeros(len_IU, 1);
    
    % Store data
    idx = 1;
    for i = 2:length(IZ)
        IZ_with_p(idx:idx+length(IZ{i})-1) = IZ{i};
        idx = idx + length(IZ{i});
    end
    
    idx = 1;
    for i = 2:length(IU)
        IU_with_p(idx:idx+length(IU{i})-1) = IU{i};
        idx = idx + length(IU{i});
    end
    
    
end


