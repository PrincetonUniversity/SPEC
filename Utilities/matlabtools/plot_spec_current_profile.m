function plot_spec_current_profile( data, iflag, newfig, varargin )
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
%   * iflag: (0) plots Ivolume
%            (1) plots Isurf
%            (2) plots both
%   * newfig: (0) plot on current figure
%             (1) open a new figure
%             (2) Erase current figure and use it
%   * Optional input: any combination of
%           -'LineWidth', value (default 2)
%           -'Color', value (default 'r')
%           -'Marker', value (default 'none')
%           -'MarkerSize', value (default 8)
%           -'LineStyle', value (default '-')
%
% Written by A. Baillod (2020)
%

    % Check inputs
    if ~any(iflag==[0,1,2])
        error('InputError: Invalid iflag')
    end

    switch newfig
      case 0 
          hold on
      case 1 
          figure('Position', [200 200 900 700],'Color','w')
          hold on
      case 2
          hold off
      otherwise
          error('InputError: invalid newfig')
    end


    l = length(varargin);
    if mod(l,2)~=0
        error('Invalid number of argument')
    end

    opt.LineWidth = 2;
    opt.Color = 'r';
    opt.Marker = 'none';
    opt.MarkerSize = 8;
    opt.LineStyle = '-';
    for ii=1:l/2
        field = varargin{2*ii-1};
        value = varargin{2*ii  };

        opt.(field) = value;
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
        start = -1;
        if lvol==1
            start=-0.999;
        end

        for k=kstart:kstart-1+nptrj(lvol)
            tflux = get_spec_torflux(data,lvol,0,start,sval(k),ns,nt);
            psitor(k) = cumflux + tflux;

            if iflag==0 || iflag==2
                Iphi(k) = cumcur + mu(lvol) * tflux;
            else
                Iphi(k) = cumcur;
            end
        end
        cumflux = psitor(k);
        kstart  = kstart+nptrj(lvol);
        %cumcur  = Iphi(k);

        if iflag==1 || iflag==2
            cumcur = Iphi(k) + data.output.IPDt(lvol); % add surface current
        else
            cumcur = Iphi(k);
        end
    end

    phiedge=psitor(end);


    plot( sqrt(psitor / phiedge), Iphi, 'LineWidth', opt.LineWidth, 'Color', opt.Color, 'Marker', opt.Marker, 'LineStyle', opt.LineStyle )
    xlabel('$(\Psi_t / \Psi_{edge})^{1/2}$', 'Interpreter', 'latex')
    ylabel('$\mu_0I_\phi$[Tm]', 'Interpreter', 'latex')
    set(gca, 'FontSize', 18)
    hold on;
  
end