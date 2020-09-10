function plot_spec_Bfield(data, component, theta, phi, nr, newfig)
%
% Plot SPEC magnetic field solution
%
% INPUT
% -----
%   data:   	data obtained from read_spec(filename)
%   component:  ='psi' to plot r-component, 'theta' to plot theta component
%               and 'phi' to plot phi component, ='all' for all components
%   theta:      Angle theta at which the field is plotted
%   phi:        Angle phi at which the field is plotted
%   nr:         Number of radial points
%   newfig:     open (=1) a new figure or use the current figure and hold
%               on (=0) or old off (=2)
%
%
% Written by A.Baillod (2019)

    switch newfig
        case 0
            hold on;
        case 1
            figure
            hold on;
        case 2
            hold off;
    end
    
    fdata = fdata_from_data(data);
    
    [r_end, z_end] = get_spec_radius(data, theta, phi, fdata.Nvol);
    [r_start, z_start] = get_spec_radius(data, theta, phi, 0);
    a = sqrt((r_end-r_start)^2 + (z_end-z_start)^2);
    
    r = linspace(0, a, nr);
    B = get_full_field(data, r, theta, phi, nr);
    
    switch component
        case 'psi'
            plot(r, B(1,:))
            ylab = ylabel('$B_\psi$ [T]');
        case 'theta'
            plot(r, B(2,:))
            ylab = ylabel('$B_\theta$ [T]');
        case 'phi'
            plot(r, B(3,:))
            ylab = ylabel('$B_\zeta$ [T]');
        case 'all'
            plot(r, B(1,:))
            hold on;
            plot(r, B(2,:))
            plot(r, B(3,:))
            ylab = ylabel('B [T]');
            leg = legend('$B_\psi$', '$B_\theta$','$B_\zeta$' );
            set(leg,'Interpreter','latex');
    end
    
    xlab = xlabel('Distance to magnetic axis [m]');
    
    set(gca, 'FontSize', 14)
    set(xlab,'Interpreter','latex');
    set(ylab,'Interpreter','latex');
    
end