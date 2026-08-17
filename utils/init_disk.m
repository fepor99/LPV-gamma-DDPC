function [params, bounds, sys, dim] = init_disk(noise)
    %   Initialize all the parameters of the disk system
    %
    %   Inputs
    %           noise:      struct containing noise parameters
    % 
    %   Outputs
    %           params:     struct containing system parameters
    %           bounds:     struct containing saturation bounds
    %           sys:        struct containing state-space matrices
    %           dim:        struct containing system dimensions

    
    % System parameters
    m       = 0.076;        % lumped mass [kg]
    J       = 2.4e-4;       % disk inertia [kg m^2]
    Km      = 11;           % motor constant [-]
    l       = 0.041;        % lumped mass - disk center distance [m]
    g       = 9.81;         % gravitational acceleration [m/s^2]
    taum    = 0.4;          % lumped back EMF constant [-]
    pmaxo   = 1;            % maximum scheduling value [-]
    pmino   = -0.22;        % minimum scheduling value [-]
    Ts      = 0.01;         % sampling time [s]

    % Save to struct
    params.m        = m;
    params.J        = J;
    params.Km       = Km;
    params.l        = l;
    params.g        = g;
    params.taum     = taum;
    params.pmaxo    = pmaxo;
    params.pmino    = pmino;
    params.Ts       = Ts;

    % Input, output, and (normalized) scheduling bounds
    bounds.u_max    = 10;           % maximum input voltage [V]
    bounds.u_min    = -10;          % minimum input voltage [V]
    bounds.y_max    = pi;           % maximum angular position [rad]
    bounds.y_min    = -pi;          % minimum angular position [rad]
    bounds.p_max    = 1;            % scaled maximum scheduling value [-]
    bounds.p_min    = -1;           % scaled minimum scheduling value [-]

    % State space matrices
    sys.A{1} = eye(2) + [0, 1; m*g*l / J * (pmaxo + pmino) / 2, -1/taum] * Ts;
    sys.A{2} = [0, 0; m*g*l / J * (pmaxo - pmino) / 2, 0] * Ts;
    sys.B{1} = [0, Km / taum]' * Ts;
    sys.B{2} = [0, 0]';
    sys.C{1} = [1, 0];
    sys.C{2} = [0, 0];
    sys.D{1} = 0; 
    sys.D{2} = 0;

    % Process noise matrices
    if noise.process
        sys.K{1}    = noise.K{1};
        sys.K{2}    = noise.K{2};
    else
        sys.K{1}    = [0, 0]';
        sys.K{2}    = [0, 0]';
    end
    
    % System dimensions
    [dim.nx, dim.nu] = size(sys.B{1});
    dim.ny = height(sys.C{1});
    dim.np = length(sys.A);
    
end


