function A_house = House_apply(A, H)
    %   Given a matrix of saved Householder vectors of a matrix Z, i.e. 
    %
    %       Z * Householder matrix = [lower triangular matrix, zeros] 
    %
    %   this function computes the product of the input matrix A times the 
    %   corresponding Householder matrix without explicitly computing the 
    %   Householder matrix
    % 
    %   Inputs
    %           A:          input matrix
    %           H:          matrix containing the LQ triangular matrix in
    %                       its lower triangular part and the Householder
    %                       vectors (without the unitary first element) in
    %                       the upper triangular part (excluding the
    %                       diagonal)
    %
    %   Outputs
    %           A_house:    A * Householder matrix
    
    
    % Initialize solution
    [m, n]  = size(H);
    A_house = A;
    
    % Compute matrix product
    for j = 1:m
        
        % Get Householder vector v and its corresposing scalar b
        v           = zeros(1, n);
        v(j:end)    = [1, H(j, j+1:end)];
        beta        = 2 / (1 + norm(H(j, j+1:end))^2);
        
        % Update A_house
        A_house(:,j:end) = A_house(:,j:end) - (A_house(:,j:end) * v(j:end).') * (beta * v(j:end));
    end

end


