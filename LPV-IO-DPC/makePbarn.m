function Pbarn = makePbarn(p,n)
%function Pbarn = makePbarn(p,n)
% Creates the block-diagonal kronecker matrix of future scheduling signals 
% Pbar^n, i.e., 
%       \bar{P}^n := (p (o) I_n) = diag( p_1 (x) I_n, ... , p_N (x) I_n )
% p is of dimension np x N, and n is an integer (nu or ny in the LPV-DPC
% paper)
    [np, N] = size(p);
    Pbarn = zeros(np*n*N,N*n);    
    for ii = 1:N
        Pbarn((1+(ii-1)*n*np):(ii*n*np),(1+(ii-1)*n):(ii*n)) = kron(p(:,ii),eye(n));
    end
end