function H = makeHankel(x,L)
% H = makeHankel(x,L)
% Determines the hankel matrix for a dataset $\{x_i\}_{i=1}^N$, where
% $x_i\in\mathbb{R}^n$. L is the length of the trajectory we consider.
% Hence, the hankel matrix will have the form
%
%      [  x1   x2   ..  x_N-L+1 ]
%      [  x2   x3   ..  x_N-L+2 ]
%  H = [  :    :     \ \   :    ]
%      [ x_L  x_L+1 ...   x_N   ]

[n,N] = size(x);
Hs = hankel(1:L,L:N);
Hr = reshape(Hs,1,[]);
xr = x(:,Hr);
H = reshape(xr, n*L,[]);
end

