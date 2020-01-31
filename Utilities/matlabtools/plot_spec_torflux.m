function plot_spec_torflux(data, zeta, cumulative, newfig)
%
%
% PLOT_SPEC_TORFLUX( FILENAME, CUMULATIVE )
% -----------------------------------------
%
% Plots the toroidal flux from the output file filename, in a cumulative or
% non-cumulative way.
%
% INPUTS
% ------
%   data:       data obtained from read_spec(filename)
%   zeta:       Toroidal angle
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


fdata = fdata_from_data(data);
Nvol = fdata.Nvol;


torflux = zeros(1,Nvol);

for lvol=1:Nvol
    tmp = get_spec_torflux(fdata,lvol,zeta,-1,1,64,64);
    
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
        error('Unsupported newfig value')
end

bar(torflux)
xlabel('Volume label')
ylabel('Toroidal flux')
set(gca, 'FontSize', 14)
xticks(1:1:Nvol)
grid on;

end