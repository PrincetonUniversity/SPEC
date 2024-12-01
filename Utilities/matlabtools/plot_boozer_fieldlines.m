%% plot_boozer_fieldlines( COORDB, NEWFIG )
% =============================================
% 
% Plots a magnetic fieldline in Boozer coordinates
% starting from the point (thetab, phib) = (0,0)
%
% INPUT
% -----
%   
%   -coordb   : must be produced using get_boozer_coordinates
%   -Newfig   : opens (=1) or not (=0) a new figure, or overwrites (=2)
%               last plot
%
% ------------------------------------%
% Written by S.Guinchard (05/15/22)   %
% ------------------------------------%

function plot_boozer_fieldlines(coordb, Newfig)

    phiB   = coordb.phib;
    thetaB = coordb.thetab;
    Fs = 18;
    switch Newfig
       
       case 0
          hold on 
            scatter(phiB, thetaB, 'k.')
            ax = gca;
            ax.TickLabelInterpreter = 'latex';
            xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            xlim([0 2*pi])
            ylim([0 2*pi])
            set (gca, 'fontsize', Fs)
            xlabel('$\phi_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
            ylabel('$\theta_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
            
       case 1
            figure
            scatter(phiB, thetaB, 'k.')
            ax = gca;
            ax.TickLabelInterpreter = 'latex';
            xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            xlim([0 2*pi])
            ylim([0 2*pi])
            set (gca, 'fontsize', Fs)
            xlabel('$\phi_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
            ylabel('$\theta_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
    
      case 2
         hold off
            scatter(phiB, thetaB, 'k.')
            ax = gca;
            ax.TickLabelInterpreter = 'latex';
            xticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            xticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            yticks([0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3 2*pi])
            yticklabels({'$0$','$\frac{\pi}{3}$','$\frac{2\pi}{3}$','$\pi$', '$\frac{4\pi}{3}$','$\frac{5\pi}{3}$','$2\pi$'});
            xlim([0 2*pi])
            ylim([0 2*pi])
            set (gca, 'fontsize', Fs)
            xlabel('$\phi_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
            ylabel('$\theta_B$', 'FontSize', Fs+10 , 'Interpreter', 'latex')
    
     end

end