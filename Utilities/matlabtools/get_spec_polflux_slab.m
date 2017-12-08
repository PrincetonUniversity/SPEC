function psipol = get_spec_polflux_slab(fdata,lvol,theta,start,send,ns,nz)

% Computes total enclosed poloidal flux in the surface defined by theta, in slab geometry
% inside the volume number lvol and across the radial extension defined by start and send
%
% INPUT
%   -fdata   : must be produced by calling read_spec_field(filename)
%   -lvol    : volume number
%   -theta   : poloidal angle at which the flux is calculated
%   -start   : first point in the radial direction (e.g. start=-1)
%   -send    : last point in the radial direction (e.g. send=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nz      : toroidal resolution (e.g. 64)
%
% OUPUT
%   -psipol  : total enclosed poloidal flux 
%
%   written by J.Loizu (2017)

sarr = linspace(start,send,ns);

zarr = linspace(0,2*pi,nz);

ds   = sarr(2)-sarr(1);

dz   = zarr(2)-zarr(1);


% Get B^{theta}

Bcontrav = get_spec_magfield_slab(fdata,lvol,sarr,theta,zarr);

Btheta    = Bcontrav{2};


% Get Jacobian of the coordinates

jac      = squeeze(get_spec_jacobian_slab(fdata,lvol,sarr,theta,zarr));


% Compute surface integral

psipol   = sum(sum( jac(2:end,:).*Btheta(2:end,:) ))*ds*dz;
