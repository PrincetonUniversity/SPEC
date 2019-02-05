function psitor = get_spec_torflux_cyl(fdata,lvol,zeta,s_start,s_end,ns,nt)

% Computes total enclosed toroidal flux in the poloidal cross-section defined by zeta,
% inside the volume number lvol and across the radial extension defined by start and send
% in cylindrical geometry
%
% INPUT
% -----
%   -fdata   : must be produced by calling read_spec_field(filename)
%   -lvol    : volume number
%   -zeta    : toroidal angle at which the flux is calculated
%   -s_start   : first point in the radial direction (e.g. start=-1)
%   -s_end    : last point in the radial direction (e.g. send=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nt      : poloidal resolution (e.g. 64)
%
% OUPUT
% -----
%   -psitor  : total enclosed toroidal flux
%
%   written by A.Baillod (2019)


sarr = linspace(s_start,s_end,ns);

tarr = linspace(0,2*pi,nt);

ds   = sarr(2)-sarr(1);

dth  = tarr(2)-tarr(1);

if(ds==0 || dth==0)

 psitor = 0;
 
else

 % Get B^{zeta}

 Bcontrav = get_spec_magfield_cyl(fdata,lvol,sarr,tarr,zeta);

 Bzeta    = Bcontrav{3};


 % Get Jacobian of the coordinates

 jac      = get_spec_jacobian_cyl(fdata,lvol,sarr,tarr,zeta);


 % Compute surface integral

 psitor   = sum(sum( jac(2:end,:).*Bzeta(2:end,:) ))*ds*dth;

end
