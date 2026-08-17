function ydpcf = LPV_IO_DPC(params, dim, bounds, sys, plt, optim, lambda_sigma, lambda_g, e)
    %% load System and simulation parameters
    ts = params.Ts;
    m = params.m;
    J = params.J;
    Km = params.Km;
    l = params.l;
    g = params.g;
    taum = params.taum;
    pmaxo = params.pmaxo;
    pmino = params.pmino;
    nx = dim.nx;
    nu = dim.nu;
    np = dim.np-1;
    Inx = eye(nx);
    Inp = eye(np);
    Inu = eye(nu);


    % An LPV representation
    pmax  = bounds.p_max;
    pmin  = bounds.p_min;
    ny = size(sys.C{1},1);
    nv     = 2^np;

    % sch. space and  vertices
    pv     = pvec('box',[pmin pmax]);
    pvert  = polydec(pv);

    % output and control constraints and Polyhedral sets
    ymax = bounds.y_max;
    ymin = bounds.y_min;

    Gy  = [eye(ny); -eye(ny)];
    hy  = [ymax  ; -ymin  ];
    Gu  = [eye(nu); -eye(nu)];
    hu  = [bounds.u_max ; bounds.u_max];
    Ymb = Polyhedron(Gy, hy);    % output constraints
    Umb = Polyhedron(Gu, hu);    % input   "  "


    %% generate data
    Nc = dim.T;
    tau = dim.M;
    Nd = dim.T + dim.M + dim.N - 1;

    u = plt.collection.u;
    y = plt.collection.y;
    p = plt.collection.p;

    % make data-structs:
    data.u = u;
    data.up = make_wp_signal(u,p);
    data.y = y;
    data.yp = make_wp_signal(y,p);

    Hddu  =  makeHankel(data.u, Nc+tau);
    Hddup =  makeHankel(data.up,Nc+tau);
    Hddy  =  makeHankel(data.y, Nc+tau);
    Hddyp =  makeHankel(data.yp,Nc+tau);



    %% For the terminal equality constraints:
    yr = 0;
    ur = -(m*g*l/J*sin(yr))*(taum/Km);
    yr = repmat(yr, dim.T, 1);
    ur = repmat(ur, dim.T, 1);
    u_ref = -(m*g*l/J*sin(optim.y_ref))*(taum/Km);

    Q = 16;
    R = 0.01*eye(nu);


    %% settings for simulation

    y0 = -pi/2*ones(1,tau);%[-pi   -pi]; % k-2 k-1
    u0 = -(m*g*l/J*sin(y0))*(taum/Km); % k-2 k-1
    p0 = sinc(y0/pi);% k-2 k-1

    nsim = optim.nsim;


    %% regularization parameters
    PSI = kron(eye(dim.T), R);
    OMEGA = kron(eye(dim.T), Q);
    A = chol(blkdiag(OMEGA, PSI)) * [Hddy((1+tau*ny):end,:); Hddu((1+tau*nu):end,:)];
    b = chol(blkdiag(OMEGA, PSI)) * [yr; ur];


    %% DPC
    % frozen p over N
    yini = y0;
    pini = p0;
    uini = u0;
    yini_real = y0;
    pini_real = p0;

    udpcf   = zeros(nu, nsim+Nc);
    y_real  = zeros(ny, nsim+Nc);
    ydpcf   = zeros(ny, nsim+Nc);
    p_real  = zeros(np, nsim+Nc);
    pdpcf   = zeros(np, nsim+Nc);
    Jdpcf = zeros(nsim+Nc,1);
    
    fprintf('\nStarting simulation with Chris'' controller\n')
    for k = 1:nsim
        y_real(:,k) = (ts^2*Km/taum)*uini(:,end-1) - (1-ts^2*m*g*l/J*pini_real(:,end-1)-ts/taum)*yini_real(:,end-1) - (ts/taum-2)*yini_real(:,end);
        ydpcf(:,k) = y_real(:,k) + e(k);
        p_real(:,k) = sinc(y_real(:,k)/pi);
        pdpcf(:,k) = sinc(ydpcf(:,k)/pi);
        Pk =  kron(ones(1,Nc),pdpcf(:,k));
        
        fprintf('Solving optimization at time step %-10d', k)
        [udpcf(:,k),Jdpcf(k,:),~,~,diagnostics] = lpviodpc_multi_noise(Hddy, Hddyp, Hddu, Hddup, ...
                                                    uini, pini, yini, ur, ...
                                                    yr, Ymb, Umb, Q, R, ...
                                                    Nc, Pk, lambda_sigma, lambda_g);
        fprintf('Done!\n')
        if diagnostics.problem ~= 0
            error('infeasible initial condition check the domain of attraction');
        end
        yini = [yini(:,2:end), ydpcf(:,k)];
        pini = [pini(:,2:end), pdpcf(:,k)];
        uini = [uini(:,2:end), udpcf(:,k)];
        yini_real = [yini_real(:,2:end), y_real(:,k)];
        pini_real = [pini_real(:,2:end), p_real(:,k)];
    end

end


