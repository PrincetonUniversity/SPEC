function plot_spec_modB_boundary(data,interface,innout,nt,nz,dimension)

%
% PLOT_SPEC_MODB_BOUNDARY( DATA, INTERFACE, INNOUT, NT, NZ, DIMENSION )
% =====================================================================
%
% Produces plot of |B| on the full boundary surface in toroidal geometry.
%
% INPUT
% -----
%   -data   : data obtained via read_spec(filename)
%   -interface    : Interface on which modB should be plotted
%   -innout: (0) - inner side of interface
%            (1) - outer side of interface
%   -nt     : poloidal resolution for the plotting (e.g. nt=64)
%   -nz     : toroidal resolution for the plotting (e.g. nz=64)
%   -dimension: (2) - plot a colorplot on a theta-phi grid
%               (3) - plot a colorplot on a 3d surface
%
% written by J.Loizu (2016)
% modified by J.Loizu (01.2017)

    Mvol = data.output.Mvol;
    if interface<1 || interface>Mvol
        error('InputError: invalid interface')
    end

    switch innout
        case 0        
            vol = interface;
            sarr = 1;
        case 1
            if interface==Mvol
                error('InputError: Cannot plot on the outer side of last interface!')
            end
            vol = interface+1;
            sarr = -1;
        otherwise
            error('InputError: interface should be 0 or 1')
    end

    if nt<1
        error('InputError: invalid nt')
    end
    if nz<1
        error('InputError: invalid nz')
    end
    if ~any(dimension==[2,3])
        error('InputError: invalid dimension')
    end

    tarr = linspace(0,2*pi,nt);
    zarr = linspace(0,2*pi,nz);

    if(vol>Mvol || vol<1)
      error('vol not valid')
    end

if(vol>Mvol)
    % Compute |B|
end

    modB   = get_spec_modB(data,vol,sarr,tarr,zarr);

    switch dimension
        case 2
            figure

            [t,z] = meshgrid(tarr, zarr);

            modB = reshape(modB, nt, nz);

            pcolor( z, t, modB' )
            shading interp

            colorbar
            xlabel('$\phi$', 'Interpreter', 'latex')
            ylabel('$\theta$', 'Interpreter', 'latex')

            ax = gca;
            ax.YTick = [pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
            ax.YTickLabel = {'$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$', '$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
            ax.XTick = [pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4];
            ax.XTickLabel = {'$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$', '$\frac{3\pi}{2}$','$\frac{7\pi}{4}$'};
            ax.TickLabelInterpreter = 'latex';

            set(gca,'FontSize',18)
            set(gcf,'Position',[200 200 900 700])


        case 3
            % Compute function (R,Z)(s,theta,zeta)
            R = get_spec_R_derivatives(data,vol,sarr,tarr,zarr,'R');
            Z = get_spec_R_derivatives(data,vol,sarr,tarr,zarr,'Z');

            R = squeeze(R{1});   
            Z = squeeze(Z{1});


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

            h=surf(X,Y,Z,squeeze(modB(1,:,:)));

            axis equal
            shading interp
            colorbar
            title('| B |')

        otherwise
            error('Invalid dimension')
    end
end
