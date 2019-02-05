function psipol = get_spec_polflux_cyl(fdata,lvol,theta,s_start,s_end,ns,nz)

% Computes total enclosed poloidal flux in the surface defined by theta
% inside the volume number lvol and across the radial extension defined by 
% s_start and s_end. Cylindrical geometry.
%
% INPUT
% -----
%   -fdata   : must be produced by calling read_spec_field(filename)
%   -lvol    : volume number
%   -theta   : poloidal angle at which the flux is calculated
%   -s_start : first point in the radial direction (e.g. s_start=-1)
%   -s_end   : last point in the radial direction (e.g. s_end=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nz      : toroidal resolution (e.g. 64)
% OUPUT
% -----
%   -psipol  : total enclosed poloidal flux 
%
%   written by J.Loizu (2016)

sarr = linspace(s_start,s_end,ns);

zarr = linspace(0,2*pi,nz);

ds   = sarr(2)-sarr(1);

dz   = zarr(2)-zarr(1);


% Get B^{theta}

Bcontrav = get_spec_magfield_cyl(fdata,lvol,sarr,theta,zarr);

Btheta    = Bcontrav{2};


% Get Jacobian of the coordinates

jac      = squeeze(get_spec_jacobian_cyl(fdata,lvol,sarr,theta,zarr));


% Compute surface integral

psipol   = sum(sum( jac(2:end,:).*Btheta(2:end,:) ))*ds*dz;

end
