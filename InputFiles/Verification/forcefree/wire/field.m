function [Br, Bz] = field( R0, I, mu0, R, Z )
%Calculate the field by a current carrying wire
% R0    - major radius of the wire
% I     - current of the wire
% mu0   - mu0
% R     - R of the measurement point
% Z     - Z of the measurement point

nint = 10000; % number of integration points

t = linspace(0, 2*pi, nint+1);

ind1 = cos(t)./((R-R0*cos(t)).^2 + (R0 * sin(t)).^2 + Z^2).^(3/2);
ind2 = 1 ./((R-R0*cos(t)).^2 + (R0 * sin(t)).^2 + Z^2).^(3/2);

int1 = trapz(t, ind1);
int2 = trapz(t, ind2);

Br = R0 * Z * int1 * I * mu0 / 4 / pi;
Bz = R0 * (R0 * int2 - R * int1) * I * mu0 / 4 / pi;

end

