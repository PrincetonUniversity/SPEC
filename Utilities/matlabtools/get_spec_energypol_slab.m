function W = get_spec_energypol_slab(data,lv,ns,nt,nz)
 
 
% Calculates Plasma Poloidal Magnetic Energy in volume lv in slab geometry 
%
% INPUT
%   -data      : must be produced by calling read_spec_field(filename)
%   -lv        : volume in which to calculate the energy (total energy for lvol=0)
%   -ns        : is the s-coordinate resolution 
%   -nt        : is the theta-coordinate resolution
%   -nz        : is the zeta-coordinate resolution
%
% OUTPUT
%   -W         : total energy 
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2019)

Nvol = data.Nvol;

W    = 0;

sarr = linspace(-1,1,ns);
tarr = linspace(0,2*pi,nt);
zarr = linspace(0,2*pi,nz);

if(lv==0)

 for lvol=1:Nvol

  modBp = get_spec_modBpol(data,lvol,sarr,tarr,zarr);

  jac   = get_spec_jacobian_slab(data,lvol,sarr,tarr,zarr);

  F     = jac.*modBp.^2;

  dW   = trapz(sarr,trapz(tarr,trapz(zarr,F,3),2));

  W    = W + 0.5*dW;

 end

else

  lvol = lv;

  modBp = get_spec_modBpol(data,lvol,sarr,tarr,zarr);

  jac   = get_spec_jacobian_slab(data,lvol,sarr,tarr,zarr);

  F     = jac.*modBp.^2;

  dW   = trapz(sarr,trapz(tarr,trapz(zarr,F,3),2));

  W    = W + 0.5*dW;

end


