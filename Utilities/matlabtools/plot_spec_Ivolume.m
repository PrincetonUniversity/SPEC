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


[psi_coord, I_vol] = get_spec_volume_current(data, cumul);


% some plots

switch newfig
    case 0
        hold on
    case 1
        figure
        hold on
    case 2
        hold off
end

%plot(psi_coord, I_vol, '*', 'DisplayName', '$I^{vol}_\phi$')
bar(I_vol);
%leg = legend('Location','northwest');
ylab = ylabel('$I_\mathcal{V}$[A]');
%xlab = xlabel('$\psi_t / \psi_{edge}$');
xlab = xlabel('Volume label');
set(gca, 'FontSize', 14)
%set(leg,'Interpreter','latex');
set(xlab,'Interpreter','latex');
set(ylab,'Interpreter','latex');
grid on;
