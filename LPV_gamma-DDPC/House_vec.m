function [v, b] = House_vec(x)
    %   Given a vector x, this function computes the Householder vector of
    %   x, i.e. it finds a scalar b and a vector v with v(1) = 1 such that
    %   P = I - b*v'*v is orthogonal and x*P = [+-norm(x), 0, ..., 0].
    %   The implementation is based upon the the procedure described in 
    %   G.H. Golub and C.F. Van Loan (1996), "Matrix computations" 
    %   (3rd ed.), Chapter 5.1.3, with some corrections to the code, also 
    %   to implement the algorithm for row vectors instead of column 
    %   vectors.
    % 
    %   Inputs
    %           x:          initial vector
    %
    %   Outputs
    %           v:          Householder vector
    %           b:          scalar multiplier
    
    
    % Define variables
    v   = [1, x(2:end)];
	sig = x(2:end) * x(2:end).';
    
    % Compute Householder vector
    if (length(x) <= 1) || (sig == 0)
        b = 2;
    else
        mu = sqrt(x(1)^2 + sig);
        if x(1) <= 0
            v(1) = x(1) - mu;
        else
            v(1) = -sig / (x(1) + mu);
        end
        b = 2 * v(1)^2 / (sig + v(1)^2);
        v = v / v(1);
    end

end


