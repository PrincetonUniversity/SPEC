function plot_SFL_modB_boundary(data,interface,innout,nt,nz,dimension)

%
% PLOT_SPEC_MODB_BOUNDARY( BDATA, DATA, VOL, NT, NZ )
% ============================================
%
% Produces plot of |B|_b on the full boundary surface in toroidal geometry.
%
% INPUT
% -----
%   -bdata   : data obtained via read_boozer(filename, root)
%   -data    : data obtained via read_spec(filename)
%   -interface    : Volume on which modB should be plotted
%   -innout: (0) - inner side of interface
%            (1) - outer side of interface
%   -nt     : poloidal resolution for the plotting (e.g. nt=64)
%   -nz     : toroidal resolution for the plotting (e.g. nz=64)
%   -dimension: (2) - plot a colorplot on a theta-phi grid
%               (3) - plot a colorplot on a 3d surface
%
% ------------------------------------%
% Written by S.Guinchard (05/30/22)   %
% ------------------------------------%
%
% TODO: adapt for more than 1 volume 
%       when booz_xform has been adapted
innout =0;
interface = 1;
Mvol = data.output.Mvol;

if interface<1 || interface>Mvol
    error('Invalid interface')
end

switch innout
    case 0        
        vol  = interface;
        sarr = 1;
    case 1
        if interface==Mvol
            error('Cannot plot on the outer side of last interface!')
        end
        vol  = interface+1;
        sarr = -1;
    otherwise
        error('Interface should be 0 or 1')
end

tarr = linspace(0,2*pi,nt);
zarr = linspace(0,2*pi,nz);

if(vol>Mvol || vol<1)
  error('vol not valid')
end


switch dimension
    case 2

        plot_SFL_modB(data,nt,nz,1,1)
        xlabel('$\phi_b$', 'Interpreter', 'latex')
        ylabel('$\theta_b$', 'Interpreter', 'latex')  
        set(gca,'FontSize',20)
        set(gcf,'Position',[200 200 900 700])
        
    
    case 3
        
        % Compute |B|
        modB   = get_SFL_fourier_modB(data,nt,nz);

        % Compute function (R,Z)(s,theta,zeta)
        sarr = 1;
        rzdata = get_spec_rzarr(data,vol,sarr,tarr,zarr);

        R = squeeze(rzdata{1});   
        Z = squeeze(rzdata{2});


        % Construct cartesian corrdinates 

        X = zeros(nt,nz);
        Y = zeros(nt,nz);

        for it=1:nt
         for iz=1:nz
          X(it,iz) = R(it,iz)*cos(zarr(iz));
          Y(it,iz) = R(it,iz)*sin(zarr(iz));
         end
        end


        % Plot

        figure
        h=surf(X,Y,Z,squeeze(modB.modB));
        axis equal
        shading interp
        colorbar
        title('$| B_b |$', 'interpreter', 'latex')
        
    otherwise
        error('Invalid dimension')
end
