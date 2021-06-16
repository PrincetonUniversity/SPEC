function plot_spec_current_profile( data, newfig )
%
% PLOT_SPEC_CURRENT_PROFILE( DATA, NEWFIG )
% =========================================
%
% Plot the total enclosed toroidal current as a function of the minor
% radius for a given SPEC equilibrium.
%
% INPUTS
% ------
%   * data: SPEC output data obtained with read_spec(filename)
%   * newfig: (0) plot on current figure
%             (1) open a new figure
%             (2) Erase current figure and use it
%
% Written by A. Baillod (2020)
%


  switch newfig
      case 0 
          hold on
      case 1 
          figure('Position', [200 200 900 700])
          hold on
      case 2
          hold off
  end
  
  nsucctrj = length(data.poincare.R(:,1,1)); % number of successfully followed trajectories
  sval    = data.transform.fiota(1:nsucctrj,1);
  nvol    = data.input.physics.Nvol;
  mvol    = data.output.Mvol;
  nptrj   = zeros(1,mvol);
  count   = 1;
  
  ind = find(sval==-2); %Remove wrongly written data
  sval(ind) = [];
  
  for is=1:length(sval)-1
   if(sval(is)>0 && sval(is+1)<0) %Reached the end of a volume
    if(count==1) 
    nptrj(count) = is;
    else
    nptrj(count) = is-sum(nptrj);
    end
    count        = count+1;
   end
  end
  nptrj(mvol)    = length(sval)-sum(nptrj);


  ns       = 32;
  nt       = 32;
  cumflux  = 0;
  cumcur   = 0;
  kstart   = 1;
  psitor   = zeros(1,length(sval));

  Iphi = zeros(1, length(sval) );
  mu = data.output.mu;

  for lvol=1:mvol
      
      for k=kstart:kstart-1+nptrj(lvol)
          tflux = get_spec_torflux(data,lvol,0,-1,sval(k),ns,nt);
          psitor(k) = cumflux + tflux;
          
          Iphi(k) = cumcur + mu(lvol) * tflux;
      end
      cumflux = psitor(k);
      kstart  = kstart+nptrj(lvol);
      %cumcur  = Iphi(k);
      
      cumcur = Iphi(k) + data.output.IPDt(lvol); % add surface current
  end
  
  phiedge=psitor(end);
  
  
  plot( sqrt(psitor / phiedge), Iphi, 'LineWidth', 3 )
  xlabel('$(\Psi_t / \Psi_{edge})^{1/2}$', 'Interpreter', 'latex')
  ylabel('$\mu_0I_\phi$[Tm]', 'Interpreter', 'latex')
  set(gca, 'FontSize', 18)
  hold on;
  
end