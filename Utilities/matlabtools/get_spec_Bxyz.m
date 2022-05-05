function xyz = get_spec_Bxyz( data, lvol, sarr, tarr, zarr )
%
% GET_SPEC_BXYZ( DATA, LVOL, SARR, TARR, ZARR )
% =============================================
%
% Write the cartesian components of the magnetic field on a 
% (x,y,z) grid on the outer side of the plasma-vacuum interface
% given a SPEC equilibrium.
%
% INPUTS
% ------
%   * DATA: data read from spec output via read_spec( fname )
%   * LVOL: volume
%   * SARR: Array of s-coordinate
%   * TARR: Array of poloidal angles
%   * ZARR: Array of toroidal angles
%
% Written by A. Baillod (2021)
%

% Check input
Igeometry = data.input.physics.Igeometry;
if Igeometry~=3
    error('Invalide geometry')
end

% Read some input from data
Mvol = data.output.Mvol;

if( lvol<1 || lvol>Mvol )
  error('invalid volume')
end

% Reshape input
ns = length(sarr);
sarr = reshape( sarr, 1, ns);

nt = length(tarr);
theta = reshape( tarr, 1, nt );

nz = length(zarr);
phi = reshape( zarr, 1, nz );

Rarr = get_spec_R_derivatives( data, Mvol, sarr, theta, phi, 'R' );
Zarr = get_spec_R_derivatives( data, Mvol, sarr, theta, phi, 'Z' );


% Transform to cartesian coordinates (x,y,z)
x = zeros(ns,nt,nz);
y = zeros(ns,nt,nz);
for is = 1:ns
  for it = 1:nt
    Rtmp = reshape(Rarr{1}(is,it,:), 1, nz);

    x(is,it,:) = Rtmp.*cos(phi);
    y(is,it,:) = Rtmp.*sin(phi);
  end
end
z = Zarr{1};
z = reshape(z, ns, nt, nz);

% Now evaluate field on each grid point
Bcontrav = get_spec_magfield( data, Mvol, sarr, theta, phi );

Bx = zeros(ns,nt,nz);
By = zeros(ns,nt,nz);
Bz = zeros(ns,nt,nz);

for is = 1:ns
  for it=1:nt
    Bs   = reshape( Bcontrav{1}(is,it,:), 1, nz);
    Bt   = reshape( Bcontrav{2}(is,it,:), 1, nz);
    Bphi = reshape( Bcontrav{3}(is,it,:), 1, nz);

    R  = reshape( Rarr{1}(is,it,:), 1, nz );
    Rs = reshape( Rarr{2}(is,it,:), 1, nz );
    Rt = reshape( Rarr{3}(is,it,:), 1, nz );
    Rz = reshape( Rarr{4}(is,it,:), 1, nz );

    Zs = reshape( Zarr{2}(is,it,:), 1, nz );
    Zt = reshape( Zarr{3}(is,it,:), 1, nz );
    Zz = reshape( Zarr{4}(is,it,:), 1, nz );

    Bx(is,it,:) = Bs.*(Rs.*cos(phi)              ) ...
                + Bt.*(Rt.*cos(phi)              ) ...
                + Bphi.*(Rz.*cos(phi) - R.*sin(phi));
    By(is,it,:) = Bs.*(Rs.*sin(phi)              ) ...
                + Bt.*(Rt.*sin(phi)              ) ...
                + Bphi.*(Rz.*sin(phi) + R.*cos(phi));
    Bz(is,it,:) = Bs.*Zs ...
                + Bt.*Zt ...
                + Bphi.*Zz;
  end

xyz.theta = theta;
xyz.phi   = phi;
xyz.x  = x;
xyz.y  = y;
xyz.z  = z;
xyz.Bx = Bx;
xyz.By = By;
xyz.Bz = Bz;

end
