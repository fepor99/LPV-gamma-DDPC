function [u0k,costopt,Yall,Uall,diagnostics] = lpviodpc_multi_noise(Hy, Hyp, Hu, Hup,uini, pini, yini,ur,yr,Ycon,Ucon,Q,R,Nc,pN, lamb_sig, lamb_g)
% LPV-IO-DPC with terminal equality constraints given p over N
%
% X1: Shifted state matrix from the data-dictionary with a proper length
% Gcal: A data-driven representation using the data-dictionary s.t. rank(Gcal) < (1+np)*(nx+nu)
% Q, R: tuning parameters
% X,U: state and input constraints as polyhdra in H-representation
% Xf: maximum RPI set as a polyhdron in H-representation
% P : Lyapunov matrix used to compute the terminal cost
% Diagnostics: related to the optimization solution
%
% 9.4.2023
%


ny = size(Ycon.A,2);
nu = size(Ucon.A,2);
np = size(pini,1);
tau = size(uini, 2);

Hty = Hy(1:(tau*ny),:);
HNy = Hy((1+tau*ny):end,:);
Htu = Hu(1:(tau*nu),:);
HNu = Hu((1+tau*nu):end,:);

Htyp = Hyp(1:(tau*ny*np),:);
HNyp = Hyp((1+tau*ny*np):end,:);
Htup = Hup(1:(tau*nu*np),:);
HNup = Hup((1+tau*nu*np):end,:);

Pbary = makePbarn(pini, ny);
Pbaru = makePbarn(pini, nu);
Pky = makePbarn(pN, ny);
Pku = makePbarn(pN, nu);

% 
yalmip('clear');
options   = sdpsettings('solver', 'gurobi', 'verbose', 0); 
g         = sdpvar(size(Hy,2),1);
sigm_tau  = sdpvar(size(yini,1)*size(yini,2), 1);
sigm_ptau = sdpvar(size(Htyp,1),1);
sigm_pN   = sdpvar(size(HNyp,1),1);
sigm_N    = sdpvar(size(HNy,1),1);
sigm_term = sdpvar(size(yini,1)*size(yini,2), 1);

sigm = [sigm_tau; sigm_ptau; sigm_pN; sigm_N;sigm_term];

%% build objective constraints
constraints = [];
objective   = lamb_g*(g'*g) + lamb_sig*(sigm'*sigm);
% initial trajectory:
constraints = [constraints;
    Hty*g == reshape(yini,[],1)+sigm_tau;
    Htu*g == reshape(uini,[],1);
    (Htyp-Pbary*Hty)*g == sigm_ptau;
    (Htup-Pbaru*Htu)*g == zeros(size(Htup,1),1)
    (HNyp-Pky*HNy)*g == sigm_pN;
    (HNup-Pku*HNu)*g == zeros(size(HNup,1),1)];
ypred = HNy*g-sigm_N;
upred = HNu*g;
% objective & x,u constraints
for i = 1:Nc
    yrange = (1+ny*(i-1)):(ny*i); 
    urange = (1+nu*(i-1)):(nu*i);
    objective = objective + (upred(urange)-ur(i))'*R*(upred(urange)-ur(i)) + (ypred(yrange)-yr(i))'*Q*(ypred(yrange)-yr(i));
        
    if i > Nc-tau % terminal equality constraints
        constraints = [constraints; upred(urange) == ur(i); ypred(yrange) == yr(i)+sigm_term((1+ny*(i-Nc+tau-1)):(ny*(i-Nc+tau)))];
    else
        constraints = [constraints; Ucon.A*upred(urange) <= Ucon.b; Ycon.A*ypred(yrange) <= Ycon.b];
    end
end

%% solve
diagnostics = optimize(constraints,objective ,options);
gval = value(g);
Uall = reshape(HNu*gval, nu,[]);
u0k = Uall(:,1);
Yall=reshape(HNy*gval, ny,[]);
costopt=double(objective);
