function plot_spec_Ivolume(data, newfig)

% plot_spec_Ivolume(filename, newfig)
%
% Plots volume current
%
% INPUT
% -----
%    filename: SPEC output filename (.sp.h5)
%    newfig  : Plots on an existing figure (=0), a new one (=1) or
%    overwrite an existing one (=2)
%
% Written by A. Baillod (2019)
%

% % Constant definition
% mu0 = 4*pi*1E-7;
% 
% % Data loading
% fdata = read_spec_field(filename);      % Read data
% Nvol = fdata.Nvol;                      % Total number of volumes
% 
% % Data processing
% 
% % First, get the current in each volume
% psi_coord = zeros(1, Nvol);             % Allocate memory
% I_vol = zeros(1, Nvol);
% 
% mu = fdata.mu;
% tflux = fdata.tflux;
% sumI = 0;
% phiedge = fdata.phiedge;
%     
% for ivol=1:Nvol
% 
%     if ivol==1    
%         I_vol(ivol) = mu(ivol) / mu0 * tflux(ivol) * phiedge;
%     else
%         % Add previous current volumes (sumI) since we use a cumulative 
%         % representation
%         I_vol(ivol) = mu(ivol) / mu0 * (tflux(ivol) - tflux(ivol-1)) * phiedge + sumI;
%     end
%     
%     psi_coord(ivol) = tflux(ivol);    
%     
%     sumI = I_vol(ivol);
% end

[psi_coord, I_vol] = get_spec_volume_current(data);


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
