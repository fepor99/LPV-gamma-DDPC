function plot_residuals(res, IZ, dim, reduction, plt)
    %   Plots the residuals of the row selection algorithm
    %   
    %   Inputs
    %           res:        residuals
    %           IZ:         selected rows of the ZP matrix
    %           dim:        struct containing system horizons
    %           reduction:  struct containing all the information for the
    %                       reduction of the ZP and UF matrices
    %           plt:        struct containing plot metadata
    

    % Change sign to past input dynamics to distinguish them from output 
    % dynamics during plotting
    idx = find(IZ - IZ(dim.M*dim.nu+1) < 0) + 1;
    res(idx) = -res(idx);

    % Set scatter plot settings
    wid = 50;           % scatter plot element size
    linew = 1.3;        % scatter plot edge linewidth
    alphae = 1;         % scatter plot edge alpha
    alphaf = 0.2;       % scatter plot face alpha

    % Initialize flag for each legend entry
    flag1 = 'on';
    flag2 = 'on';
    flag3 = 'on';
    flag4 = 'on';
    flag5 = 'on';
    flag6 = 'on';
    
    
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
    
    
    %% Generate figure
    figure;
    tiledlayout(1, 1, 'padding', 'compact', 'tilespacing', 'compact')

    nexttile;
    for i = 0:reduction.rows.u + reduction.rows.z
        
        % No rows selected
        if i == 0
            scatter(i, res(i+1), wid, 'filled', 'o', 'MarkerFaceColor', [0, 0, 0], ...
                    'MarkerEdgeColor', [0, 0, 0], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'DisplayName', '$|| Y_F ||_F^2$')
                
        % Past input LTI rows
        elseif i <= dim.nu*dim.M
            scatter(i, -res(i+1), wid, 'filled', 'd', 'MarkerFaceColor', [0.92, 0.71, 0.24], ...
                    'MarkerEdgeColor', [0.92, 0.71, 0.24], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag1, 'DisplayName', '$U_P$ without scheduling')
            if length(flag1) == 2
                flag1 = 'off';
            end
            
        % Past output LTI rows
        elseif i <= (dim.nu+dim.ny)*dim.M
            scatter(i, res(i+1), wid, 'filled', 'o', 'MarkerFaceColor', [0, 0.58, 0], ...
                    'MarkerEdgeColor', [0, 0.58, 0], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag2, 'DisplayName', '$Y_P$ without scheduling')
            if length(flag2) == 2
                flag2 = 'off';
            end
            
        % Past input LPV rows
        elseif res(i+1) < 0
            scatter(i, -res(i+1), wid+20, 'filled', 's', 'MarkerFaceColor', [0, 0, 1], ...
                    'MarkerEdgeColor', [0, 0, 1], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag3, 'DisplayName', '$U_P$ with scheduling')
            if length(flag3) == 2
                flag3 = 'off';
            end
            
        % Past output LPV rows
        elseif i <= reduction.rows.z
            scatter(i, res(i+1), wid, 'filled', '^', 'MarkerFaceColor', [1, 0, 0], ...
                    'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag4, 'DisplayName', '$Y_P$ with scheduling')
            if length(flag4) == 2
                flag4 = 'off';
            end
            
        % Future input LTI rows
        elseif i <= reduction.rows.z + dim.T*dim.nu
            scatter(i, res(i+1), wid, 'filled', 'd', 'MarkerFaceColor', [0.15, 0.7, 0.6], ...
                    'MarkerEdgeColor', [0.15, 0.7, 0.6], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag5, 'DisplayName', '$U_F$ without scheduling')
            if length(flag5) == 2
                flag5 = 'off';
            end
            
        % Future input LPV rows
        else
            scatter(i, res(i+1), wid+20, 'filled', 's', 'MarkerFaceColor', [0.53, 0.32, 0.78], ...
                    'MarkerEdgeColor', [0.53, 0.32, 0.78], 'MarkerFaceAlpha', alphaf, ...
                    'MarkerEdgeAlpha', alphae, 'LineWidth', linew, ...
                    'HandleVisibility', flag6, 'DisplayName', '$U_F$ with scheduling')
            if length(flag6) == 2
                flag6 = 'off';
            end
        end
        hold on
    end

    % Axes limits and legend
    box on
    set(gca, 'YScale', 'log', 'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
             'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.1)
    xlim([-5, reduction.rows.z + reduction.rows.u + 5])
    ylim([abs(res(end))/3, abs(res(1))*3])
    xline(reduction.rows.z+0.5, 'k-.', 'linewidth', 1.5, 'HandleVisibility', 'off')
    ylabel('Residual')
    xlabel('Rows selected')
    legend()

end





