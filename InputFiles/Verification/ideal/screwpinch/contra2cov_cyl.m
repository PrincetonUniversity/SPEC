function vcov = contra2cov_cyl(filename, vol, s, vcontrav, ns, theta, phi, norm)

%
% Transform a contravariant vector to a covariant one in cylindrical
% coordinates
%
% INPUT
% -----
%   filename:   spec hdf5 output file
%   vol:        volume
%   s:          radial coordinate
%   vcontrav:   contravariant vector as a function of r, size = 3xlength(r)
%   ns:         number of point in each volume for interpolation
%   theta:      theta angle
%   phi:        phi angle
%   norm:       use the unitary basis (=1) or the general basis(=0)
%
% OUTPUT
% ------
%   vcov:       covariant vector as a function v


[sarr, Rarr] = getR_derivatives(filename, vol, ns, theta, phi);


% transform in covariant basis
vcov = zeros(3, length(s));

R = interp1(sarr, Rarr(1,:), s);
Rs = interp1(sarr, Rarr(2,:), s);
Rtheta = interp1(sarr, Rarr(3,:), s);
Rphi = interp1(sarr, Rarr(4,:), s);
    
vcov(1,:) = Rs.^2 .* vcontrav(1,:) + Rs .* Rtheta .* vcontrav(2,:) ...
            + Rs .* Rphi .* vcontrav(3,:);
vcov(2,:) = Rs .* Rtheta .* vcontrav(1,:) + (Rtheta.^2 + R.^2) .* vcontrav(2,:)...
            + Rtheta .* Rphi .* vcontrav(3,:);
vcov(3,:) = Rs.*Rphi.*vcontrav(1,:) + Rtheta.*Rphi.*vcontrav(3,:)...
            + (1+Rphi.^2).*vcontrav(3,:);

if norm %Normalize with respect to basis vector norm.
   vcov(1,:) = vcov(1,:) ./ sqrt((Rtheta.^2 + R.^2 + (R.*Rphi ).^2) ./ (Rs.*R).^2);
   vcov(2,:) = vcov(2,:) ./ R;
end







