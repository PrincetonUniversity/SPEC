function Avol = get_spec_area(data,lvol,smax,ns,nt,phi0)
 
%
% GET_SPEC_AREA( DATA, LVOL, NS, NT, PHI0 )
% =========================================
%
% Calculates cross-sectional area of a given volume at fixed phi
%
% INPUT
% -----
%   -data    : must be produced by calling read_spec(filename)
%   -lvol    : volume number
%   -smax    : max s
%   -ns      : is the resolution in the s-coordinate     (e.g. 64)
%   -nt      : is the resolution in the theta-coordinate (e.g. 64)
%   -phi0    : toroidal angle defining a toroidal plane
%
% OUTPUT
% ------
%   -Avol    : area in m^2 if geometrical dimensions (R,Z) are interpreted in meters.
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2018)
%

smin     = -0.999; %avoids singular inversion of the metric matrix

% Define arrays
sarr     = linspace(smin,smax,ns);
tarr     = linspace(0,2*pi,nt);

% Evaluate Jacobian and metric
jacobian = get_spec_jacobian(data,lvol,sarr,tarr,phi0);
gcontrav = get_spec_metric_contrav(data,lvol,sarr,tarr,phi0);
sqrtgradphi = sqrt(gcontrav{3}{3});
  
% Evaluate area
Avol = trapz(tarr, trapz(sarr, jacobian.*sqrtgradphi, 1) );

end
