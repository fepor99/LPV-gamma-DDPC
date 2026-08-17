function wp = make_wp_signal(w, p)
% Makes the w^p signal, i.e., the vector
%       [ p(1) (x) w(1) ]
%       [ p(1) (x) w(1) ]
%       [  :    :   :   ]
%       [p(Nd) (x) w(Nd)]
% input of function:
%  - w of size (nw, Nd)
%  - p of size (np, Nd)

[nw,Ndw] = size(w);
[np,Ndp] = size(p);

if nw > Ndw || np > Ndp
    disp('!! WARNING !!')
    disp('Number of datapoints larger than signal dimension...')
    disp('Ensure dim({w,p})= [{nw,np},Nd]...')
end
assert(Ndw==Ndp,'Number of datasamples not equal'); 
Nd = Ndw;
wp = reshape(w'.*reshape(p',Nd,1,[]),Nd,[])';

end