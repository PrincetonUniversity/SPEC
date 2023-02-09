%% plot_spec_fieldlines( DATA, NPOINTS, NPERIODS, NEWFIG )
% =======================================================
%
% Traces magnetic field lines in the (phi, theta) plane
% and gives an output containing the coordinates 
% (used e.g with the code get_spec_straightfieldlines)
% 
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -Npoints  : number of points along the field line 
%   -Nperiods : number of toroidal periods over which the field line is
%               traced
%   -Newfig   : opens (=1) or not (=0) a new figure, or overwrites (=2)
%   last plot
%
% ------------------------------------%
% Written by S.Guinchard (03/01/22)   %
% Last modified (05/15/22)            %
% ------------------------------------%

function coord = plot_spec_fieldlines(d,Npoints,Nperiods,newfig)

 phi   = linspace(0,2*Nperiods*pi,Npoints);
 dphi  = phi(2)-phi(1);
 s     = 1;
 Fs    = 18;
 theta_temp(1) = 0;

 for i = 1:length(phi)-1
        Bfield_temp = get_spec_magfield(d,1,s,theta_temp(i),phi(i));
        theta_temp(i+1)  = theta_temp(i)+dphi*(cell2mat(Bfield_temp(2))/cell2mat(Bfield_temp(3)));
 end

 coord.theta = wrapTo2Pi(theta_temp);
 coord.phi   = wrapTo2Pi(phi);
 switch newfig
   
   case 0
      hold on 
        scatter(coord.phi, coord.theta, 'b.')
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        xlim([0 2*pi])
        ylim([0 2*pi])
        set (gca, 'fontsize', Fs)
        xlabel('$\phi$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
        ylabel('$\theta$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
        
   case 1
        figure
        scatter(coord.phi, coord.theta, 'b.')
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        xlim([0 2*pi])
        ylim([0 2*pi])
        set (gca, 'fontsize', Fs)
        xlabel('$\phi$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
        ylabel('$\theta$', 'FontSize', Fs+10 , 'Interpreter', 'latex')

  case 2
     hold off
        scatter(coord.phi, coord.theta, 'b.')
        ax = gca;
        ax.TickLabelInterpreter = 'latex';
        xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
        yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
        xlim([0 2*pi])
        ylim([0 2*pi])
        set (gca, 'fontsize', Fs)
        xlabel('$\phi$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
        ylabel('$\theta$', 'FontSize', Fs+10 , 'Interpreter', 'latex')

 end

end