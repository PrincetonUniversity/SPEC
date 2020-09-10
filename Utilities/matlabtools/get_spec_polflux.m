function psipol = get_spec_polflux(data,lvol,theta,start,send,ns,nz)

%
% GET_SPEC_POLFLUX( DATA, LVOL, THETA, START, SEND, NS, NZ )
% ==========================================================
%
% Computes total enclosed poloidal flux in the surface defined by theta
% inside the volume number lvol and across the radial extension defined by start and send
%
% INPUT
% -----
%   -data   : must be produced by calling read_spec(filename)
%   -lvol    : volume number
%   -theta   : poloidal angle at which the flux is calculated
%   -start   : first point in the radial direction (e.g. start=-1)
%   -send    : last point in the radial direction (e.g. send=+1)
%   -ns      : radial resolution   (e.g. 64)
%   -nz      : toroidal resolution (e.g. 64)
%
% OUPUT
% -----
%   -psipol  : total enclosed poloidal flux 
%
%   written by J.Loizu (2016)
%   modified by A. Baillod (2019) - Added switch for geometry

sarr = linspace(start,send,ns);

zarr = linspace(0,2*pi,nz);

ds   = sarr(2)-sarr(1);

dz   = zarr(2)-zarr(1);


Bcontrav = get_spec_magfield(data,lvol,sarr,theta,zarr);
jac      = squeeze(get_spec_jacobian(data,lvol,sarr,theta,zarr));

        
        
Btheta    = Bcontrav{2};

% Compute surface integral

psipol   = sum(sum( jac(2:end,:).*Btheta(2:end,:) ))*ds*dz;
