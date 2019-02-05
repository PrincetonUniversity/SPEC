function phi = tflux(r, B, r_kam, total)
% Compute the poloidal flux in cylindrical coordinate with toroidal and
% poloidal symmetry
%
% INPUT
% -----
%   r:      radial coordinate of B
%   B:      magnetic field
%   r_kam:  radial coordinate between which the toroidal flux is computed
%   total:  annulus (=0) or cylinder (=1) toroidal flux
%
% OUPUT
% -----
%   phi:    Toroidal flux as a function of r_flux

    Nvol = length(r_kam);
    phi = zeros(1, Nvol);
    

    for jj=1:Nvol

	  if total || jj==1
              lbound = 0;
          else
              lbound = r_kam(jj-1);
          end
            ubound = r_kam(jj);
        
        r_temp = linspace(lbound, ubound, 1E3);
        B_temp = interp1(r, B, r_temp, 'spline');
        
        if false
           figure
           plot(r_temp, B_temp, '*')
           hold on;
           plot(r, B)
        end

        phi(jj) = 2 * pi * trapz(r_temp, r_temp .* B_temp);
    end
end
