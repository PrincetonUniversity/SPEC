function Lvol = get_spec_area(data,lvol,nt,phi0)
 
%
% GET_SPEC_AREA( DATA, LVOL, NT, PHI0 )
% =====================================
%
% Calculates length of a curve defined by a volume boundary at fixed phi
%
% INPUT
% -----
%   -data    : must be produced by calling e.g. read_spec(filename)
%   -lvol    : volume number
%   -nt      : is the resolution in the theta-coordinate (e.g. 64)
%   -phi0    : toroidal angle defining a toroidal plane
%
% OUTPUT
% ------
%   -Lvol    : length in meters if geometrical dimensions (R,Z) are interpreted in meters.
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2018)

tarr        = linspace(0,2*pi,nt);

gcov        = get_spec_metric(data,lvol,1,tarr,phi0);
sqrtetheta  = sqrt(gcov{2}{2});
  
Lvol  = sum(sqrtetheta(:))*(2*pi)/nt;

end
