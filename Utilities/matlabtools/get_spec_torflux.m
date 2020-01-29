function psitor = get_spec_torflux(fdata,lvol,zeta,start,send,ns,nt)

% Computes total enclosed toroidal flux in the poloidal cross-section defined by zeta,
% inside the volume number lvol and across the radial extension defined by start and send
%
% INPUT
%   -fdata   : must be produced by calling read_spec_field(filename)
%   -lvol    : volume number
%   -zeta    : toroidal angle at which the flux is calculated
%   -start   : first point in the radial direction (e.g. start=-1)
%   -send    : last point in the radial direction (e.g. send=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nt      : poloidal resolution (e.g. 64)
%
% OUPUT
%   -psitor  : total enclosed toroidal flux
%
%   written by J.Loizu (2016)
%   modified by J.Loizu (01.2017)
%   modified by J.Loizu (06.2017)
%   modified by A.Baillod (06.2019) - added switch for geometry


sarr = linspace(start,send,ns);

tarr = linspace(0,2*pi,nt);

ds   = sarr(2)-sarr(1);

dth  = tarr(2)-tarr(1);

if(ds==0 || dth==0)

 psitor = 0;
 
else

Bcontrav = get_spec_magfield(fdata,lvol,sarr,tarr,zeta);
jac      = get_spec_jacobian(fdata,lvol,sarr,tarr,zeta);
  

 % Compute surface integral

 Bzeta    = Bcontrav{3};
 psitor   = sum(sum( jac(2:end,:).*Bzeta(2:end,:) ))*ds*dth;

end
