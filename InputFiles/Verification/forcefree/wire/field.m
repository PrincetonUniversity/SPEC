function [Br, Bz] = field( R0, mu0I, R, Z, nint )
%Calculate the field by a current carrying wire
% R0    - major radius of the wire
% mu0I  - mu0 * current of the wire
% R     - R of the measurement point
% Z     - Z of the measurement point
% nint  - number of integration points

t = linspace(0, 2*pi, nint+1);

ind1 = cos(t)./((R-R0*cos(t)).^2 + (R0 * sin(t)).^2 + Z^2).^(3/2);
ind2 = 1 ./((R-R0*cos(t)).^2 + (R0 * sin(t)).^2 + Z^2).^(3/2);

int1 = trapz(t, ind1);
int2 = trapz(t, ind2);

Br = R0 * Z * int1 * mu0I / 4 / pi;
Bz = R0 * (R0 * int2 - R * int1) * mu0I / 4 / pi;

end

