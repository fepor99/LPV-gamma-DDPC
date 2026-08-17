function plot_data_mean_std(plt)
    %   Generates all the plot related to the script. It includes the
    %   representation of the collected data and the simulation results of 
    %   the closed-loop system
    % 
    %   Inputs
    %           plt:        struct containing data to plot
    
    

    %% Extract data
    y_ref = plt.y_ref;

    % Control simulation signals
    y_MPC_mean          = mean(plt.MPC.x(1,:,:), 3);
    y_MPC_std           = std(plt.MPC.x(1,:,:), 0, 3);
    y_gamma_DPC_mean    = mean(plt.gamma_DPC.x(1,:,setdiff(1:end,plt.gamma_DPC.rem)), 3);
    y_gamma_DPC_std     = std(plt.gamma_DPC.x(1,:,setdiff(1:end,plt.gamma_DPC.rem)), 0, 3);
    y_IO_DPC_mean       = mean(plt.IO_DPC.y(1,:,setdiff(1:end,plt.IO_DPC.rem)), 3);
    y_IO_DPC_std        = std(plt.IO_DPC.y(1,:,setdiff(1:end,plt.IO_DPC.rem)), 0, 3);

    % Define plot lines
    y_MPC_low           = y_MPC_mean - y_MPC_std;
    y_MPC_high          = y_MPC_mean + y_MPC_std;
    y_gamma_DPC_low     = y_gamma_DPC_mean - y_gamma_DPC_std;
    y_gamma_DPC_high    = y_gamma_DPC_mean + y_gamma_DPC_std;
    y_IO_DPC_low        = y_IO_DPC_mean - y_IO_DPC_std;
    y_IO_DPC_high       = y_IO_DPC_mean + y_IO_DPC_std;

    % Define colored regions
    nsim            = plt.optim.nsim;
    x_range         = [1:nsim, fliplr(1:nsim)] - 1;
    MPC_range       = [y_MPC_high(1,1:nsim), fliplr(y_MPC_low(1,1:nsim))];
    gamma_DPC_range = [y_gamma_DPC_high(1,1:nsim), fliplr(y_gamma_DPC_low(1,1:nsim))];
    IO_DPC_range    = [y_IO_DPC_high(1,1:nsim), fliplr(y_IO_DPC_low(1,1:nsim))];


    %% Plot metadata
    % Set Latex interpreter and font
    set(0, 'defaultTextInterpreter', 'latex')
    set(0, 'defaultAxesTickLabelInterpreter', 'latex');
    set(0, 'defaultLegendInterpreter', 'latex');
    set(0, 'defaultTextFontName', 'CMU Serif')
    set(0, 'defaultAxesFontName', 'CMU Serif')
    set(0, 'defaultLegendFontName', 'CMU Serif');

    % Line width
    set(0, 'defaultLineLineWidth', 1);
    set(0, 'defaultAxesLineWidth', 1.5);

    % Font size
    set(0, 'defaultAxesFontSize', plt.font);

    % Define colors and other properties
    red         = [216,  27,  96]/255;
    blue        = [ 30, 136, 229]/255;
    yellow      = [255, 165,   0]/255;
    line_width  = 3;
    limit_width = 0.1;


    %% Plot control simulation
    figure('Position', [303 294 1310 420]);
    tiledlayout(1, 2, 'padding', 'compact', 'tilespacing', 'loose')

    % MPC+KF vs LPV gamma-DDPC
    nexttile
    grid on; hold on;

    % Shaded areas
    fill(x_range, MPC_range(1,:), [0.6, 0.733, 1], 'FaceAlpha', 0.3, 'FaceColor', red, ...
         'EdgeColor', red, 'LineWidth', limit_width, 'LineStyle', '--', 'HandleVisibility', 'off');
    fill(x_range, gamma_DPC_range(1,:), [0.6, 0.733, 1], 'FaceAlpha', 0.3, 'FaceColor', blue, ...
         'EdgeColor', blue, 'LineWidth', limit_width, 'LineStyle', '--', 'HandleVisibility', 'off');

    % Mean values
    plot(0:nsim-1, y_MPC_mean(1,:), 'color', red, 'linewidth', line_width, 'DisplayName', 'Oracle MPC+KF')
    plot(0:nsim-1, y_gamma_DPC_mean(1,:), 'color', blue, 'linewidth', line_width, 'DisplayName', 'LPV $\gamma$-DDPC')

    % Labels etc...
    plot(0:nsim-1, y_ref(1:nsim), 'k-.', 'linewidth', line_width-1, 'DisplayName', 'Reference')
    legend('Location', 'best')
    ylabel('$y$ [rad]')
    xlabel('Time step')
    axis padded
    xlim([0, nsim-1])
    box on

    % MPC+KF vs LPV-IO-DPC
    nexttile
    grid on; hold on;

    % Shaded areas
    fill(x_range, MPC_range(1,:), [0.6, 0.733, 1], 'FaceAlpha', 0.3, 'FaceColor', red, ...
         'EdgeColor', red, 'LineWidth', limit_width, 'LineStyle', '--', 'HandleVisibility', 'off');
    fill(x_range, IO_DPC_range(1,:), [0.6, 0.733, 1], 'FaceAlpha', 0.3, 'FaceColor', yellow, ...
         'EdgeColor', yellow, 'LineWidth', limit_width, 'LineStyle', '--', 'HandleVisibility', 'off');

    % Mean values
    plot(0:nsim-1, y_MPC_mean(1,:), 'color', red, 'linewidth', line_width, 'DisplayName', 'Oracle MPC+KF')
    plot(0:nsim-1, y_IO_DPC_mean(1,1:nsim), 'color', yellow, 'linewidth', line_width, 'DisplayName', 'LPV-IO-DPC')

    % Labels etc...
    plot(0:nsim-1, y_ref(1:nsim), 'k-.', 'linewidth', line_width-1, 'DisplayName', 'Reference')
    legend('Location', 'best')
    ylabel('$y$ [rad]')
    xlabel('Time step')
    axis padded
    xlim([0, nsim-1])
    box on

    
end


