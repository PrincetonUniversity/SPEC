function plot_spec_Ivolume(data, cumul, newfig)

% 
% PLOT_SPEC_IVOLUME( DATA, CUMUL, NEWFIG )
% ========================================
%
% Plots volume current
%
% INPUT
% -----
%    data:     Obtained via read_spec(filename)
%    cumul:    Cumulative (=1) or not cumulative (=0) quantity
%    newfig  : Plots on an existing figure (=0), a new one (=1) or
%    overwrite an existing one (=2)
%
% Written by A. Baillod (2019)
%

    if ~any(cumul==[0,1])
        error('InputError: invalid cumul')
    end

    switch newfig
        case 0
            hold on
        case 1
            figure
            hold on
        case 2
            hold off
    end


    [~, I_vol] = get_spec_volume_current(data, cumul);


    bar(I_vol);
    xlab = xlabel('Volume label');
    ylab = ylabel('$\mu_0I_\mathcal{V}$[A]');
    set(gca, 'FontSize', 14)
    set(xlab,'Interpreter','latex');
    set(ylab,'Interpreter','latex');
    grid on;
end
