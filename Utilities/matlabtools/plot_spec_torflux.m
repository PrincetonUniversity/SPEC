function plot_spec_torflux(data, zeta, cumulative, newfig)

%
% PLOT_SPEC_TORFLUX( DATA, ZETA, CUMULATIVE, NEWFIG )
% -----------------------------------------
%
% Integrates the magnetic field on a constant toroidal surface,
% in a cumulative or non-cumulative way.
%
% INPUTS
% ------
%   data:       data obtained from read_spec(filename)
%   zeta:       Toroidal angle
%   cumulative: true to get cumulative plot (\psi_a = \int_0^a B_\phi dS)
%               or false to get a non-cumulative plot (\psi_a = 
%               \int_{a-1}^a B_\phi dS)
%   newfig:     0:  plots on an existing figure without erasing previous
%                   plot
%               1:  plots on a new figure
%               2:  plots on an existing figure and erase previous plot
%
%
% Written by A. Baillod (2019)
%
    Mvol = data.output.Mvol;
    torflux = zeros(1,Mvol);
    for lvol=1:Mvol
        start=-1;
        if(lvol==1)
        start=-0.999;
        end
        tmp = get_spec_torflux(data,lvol,zeta,start,1,64,64);

Nvol = data.input.physics.Nvol;
        if cumulative
            if lvol==1
                torflux(lvol)=tmp;
            else
                torflux(lvol) = torflux(lvol-1) + tmp;
            end
        else
            torflux(lvol)=tmp;
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
            error('InputError: Invalid newfig')
    end

    bar(torflux)
    xlabel('Volume label')
    ylabel('Toroidal flux')
    set(gca, 'FontSize', 14)
    xticks(1:1:Nvol)
    grid on;

end
