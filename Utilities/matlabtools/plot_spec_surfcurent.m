function plot_spec_surfcurent(data, ns, nt, zeta, newfig)

%
% PLOT_SPEC_SURFCURRENT( DATA, NS, NT, ZETA, NEWFIG )
% ===================================================
%
% Plot pressure-driven currents located at each volume interface
%
% INPUT
% -----
%   -data       : data obtained via read_spec(filename)
%   -ns         : number of radial interpolation points
%   -nt         : number of poloidal interpolation points
%   -zeta       : toroidal angle 
%   -newfig     : plots on an existing figure (=0), a new figure (=1) or
%   overwrites last figure (=2)
%
% written by A.Baillod (2019)
%

[tflux, IPDt] = get_spec_surface_current(data, ns, nt, zeta);

Nvol = data.input.physics.Nvol + data.input.physics.Lfreebound;

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

%plot(tflux(1:end-1), IPDt, '*')
bar(IPDt, 'BarWidth', 0.3);
grid on
%xl = xlabel('$\psi_t / \psi_{edge}$');
xl = xlabel('Surface label');
yl = ylabel('$\mu_0 I_\mathcal{S}$[A]');
xlim([0, Nvol])

set(xl, 'Interpreter', 'latex');
set(yl, 'Interpreter', 'latex');

set(gca, 'FontSize', 14)


end
