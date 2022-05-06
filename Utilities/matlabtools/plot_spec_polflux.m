function plot_spec_polflux(data, theta, cumulative, newfig)
%
%
% PLOT_SPEC_POLFLUX( FILENAME, ZETA, CUMULATIVE, NEWFIG )
% -----------------------------------------
%
% Plots the poloidal flux from the output file filename, in a cumulative or
% non-cumulative way.
%
% INPUTS
% ------
%   data:       data obtained from read_spec(filename);
%   theta:      Poloidal angle
%   cumulative: True to get cumulative plot (\psi_a = \int_0^a B_\phi dS)
%               or False to get a non-cumulative plot (\psi_a = 
%               \int_{a-1}^a B_\phi dS)
%   newfig:     0:  plots on an existing figure without erasing previous
%                   plot
%               1:  plots on a new figure
%               2:  plots on an existing figure and erase previous plot
%
%
% Written by A. Baillod (2019)
%
%

    Nvol = data.input.physics.Nvol;

    polflux = zeros(1,Nvol);

    % Start at one only if slab geometry (otherwise it is zero by definition)
    sarr = linspace(-1, 1, 64);
    if data.input.physics.Igeometry==1
        polflux(1) = get_spec_polflux(data,1,theta,sarr,64);
    end

    for lvol=2:Nvol
        tmp = get_spec_polflux(data,lvol,theta,sarr,64);

        if cumulative
            polflux(lvol) = polflux(lvol-1) + tmp;
        else
            polflux(lvol)=tmp;
        end

    end


    switch newfig
        case 0
            hold on;
        case 1
            figure
            hold on;
        case 2
            hold off;
        otherwise
            error('InputError: invalid newfig')
    end

    bar(polflux)
    xlabel('Volume label')
    ylabel('Poloidal flux')
    set(gca, 'FontSize', 14)
    xticks(1:1:Nvol)
    grid on;

end
