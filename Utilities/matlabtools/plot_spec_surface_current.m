function plot_spec_surface_current(data, nt, newfig)

%
% PLOT_SPEC_SURFACE_CURRENT( DATA, NS, NT, NEWFIG )
% ===================================================
%
% Plot pressure-driven currents located at each volume interface
%
% INPUT
% -----
%   -data       : data obtained via read_spec(filename)
%   -nt         : number of poloidal interpolation points
%   -newfig     : plots on an existing figure (=0), a new figure (=1) or
%   overwrites last figure (=2)
%
% written by A.Baillod (2019)
%

    % Test input
    if nt<1
        error('InputError: Invalid nt')
    end


    switch newfig
        case 0
            hold on
        case 1
            figure
            hold on
        case 2
            hold off
        otherwise
            error('InputError: Invalid newfig')
    end

    % Evaluate toroidal current
    Itor = get_spec_torcurr_kam_net(data, nt);
    Mvol = data.output.Mvol;


    %plot(tflux(1:end-1), IPDt, '*')
    bar(Itor, 'BarWidth', 0.3);
    grid on
    %xl = xlabel('$\psi_t / \psi_{edge}$');
    xl = xlabel('Surface label');
    yl = ylabel('$\mu_0 I_\mathcal{S}$[A]');
    xlim([0, Mvol])

    set(xl, 'Interpreter', 'latex');
    set(yl, 'Interpreter', 'latex');

    set(gca, 'FontSize', 14)


end
