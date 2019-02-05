function [r_space, jp, jz, Bp, Bz] = ScrewPinch_IotaSolver(L, a, r_space, p, iota, Bz0, fig)

% Solve the screw pinch equlibrium for a given pressure and rotational
% transform profile
%
% INPUT
% -----
%   L:      Length of the screw pinch [m]
%   a:      Radius of the screw pinch [m]
%   r_space Radial coordinate
%   p:      Pressure profile as a function of r_space (not normalized by
%           mu0)
%   iota:   Rotational transform profile, multiplied by 2pi (as defined in
%           Freidberg's book)
%   Bz0:    Toroidal field at r=0
%   fig:    Plot (=1) or don't plot (=0) additional figures
%
% OUTPUT
% ------
%   r_space Radial coordinate
%   jp:     Poloidal current density
%   jz:     Toroidal current density
%   Bp:     Poloidal magnetic field
%   Bz:     Toroidal magnetic field
%
% Written by A.Baillod (2019)

    [r_space, jp, jz, Bp, Bz] = forward(L, a, r_space, p, iota, Bz0);


    if fig % Some figures to debug
        figure;

        subplot(2, 2, 1);
        plot(r_space, p);
        title('Pressure')

        subplot(2, 2, 2);
        plot(r_space, iota);
        title('rotational transform')

        subplot(2, 2, 3);
        plot(r_space, jp)
        yyaxis right
        plot(r_space, jz)
        legend('J_\theta', 'J_z')
        title('Current density')

        subplot(2, 2, 4);
        plot(r_space, Bp)
        yyaxis right
        plot(r_space, Bz)
        legend('B_\theta', 'B_z')
        title('Magnetic field')
    end
end



function [r_space, jp, jz, Bp, Bz] = forward(L, a, r_space, p, iota, Bz0)

    % Solve the ode
    r_span = [0, a];    % Space on which the ode is solved
    [r_Bz, Bz] = ode45(@(r, Bz) forward_ode(r, Bz, p, iota, r_space, L),...
                                            r_span, Bz0); 

    Bz = interp1(r_Bz, Bz, r_space, 'pchip');
    r_Bz = r_space;
    Bp = iota .* r_space .* Bz / L;

    jp = -diff(Bz);
    r_jp = r_Bz(1:end-1) + diff(r_Bz) / 2.0;
    jp = interp1(r_jp, jp, r_space, 'pchip');

    dBp = diff(Bp) ./ diff(r_space);
    drspace = r_space(1:end-1) + diff(r_space);
    dBp = interp1(drspace, dBp, r_space, 'pchip');
    jz = dBp + Bp ./ r_space;


end

function dBzdr = forward_ode(r, Bz, p, iota, r_space, L)

mu0 = 4*pi*1E-7;

% Derivatives evaluation
h = r_space(2) - r_space(1);
dpdr = diff(p) / h;
didr = diff(iota) / h;
dr_space = r_space(1:end-1) + h/2;
    
% Interpolation
int_dpdr = interp1(dr_space, dpdr, r, 'spline');
int_iota = interp1( r_space, iota, r, 'spline');
int_didr = interp1(dr_space, didr, r, 'spline');

% Evaluation
num = int_dpdr + r * int_iota^2 / (mu0 * L^2) * (2 + ...
                                           r * int_didr / int_iota) * Bz^2;
denum = (1 + (r * int_iota / L)^2) * Bz / mu0;
dBzdr = - num / denum;

end
