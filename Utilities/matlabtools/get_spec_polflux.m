function psipol = get_spec_polflux(data,lvol,theta,sarr,nz)

%
% GET_SPEC_POLFLUX( DATA, LVOL, THETA, SARR, NZ )
% ===============================================
%
% Computes total enclosed poloidal flux in the surface defined by theta
% inside the volume number lvol and across the radial extension defined by 
% start and send
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

% Check input
Igeometry = data.input.physics.Igeometry;
if (sarr(1)==-1) && (lvol==1) && (Igeometry~=1)
    error('Singularity in first volume for s=-1. Set sarr to start from >-1')
end

% Build zeta array
zarr = linspace(0,2*pi,nz);

% Get B theta contravariant and the jacobian
Bcontrav = get_spec_magfield(data,lvol,sarr,theta,zarr);
Btheta   = squeeze(Bcontrav{2});
jac      = squeeze(get_spec_jacobian(data,lvol,sarr,theta,zarr));

      
% Compute surface integral
psipol = trapz(zarr, trapz(sarr, jac.*Btheta, 1 ) );
