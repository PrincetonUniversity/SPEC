function modB = get_spec_modB(data,lvol,sarr,tarr,zarr)
 
% 
% GET_SPEC_MODB( DATA, LVOL, SARR, TARR, ZARR )
% =============================================
%
% Calculates mod(B) in volume lvol 
%
% INPUT
% -----
%   -data  : must be produced by calling read_spec(filename)
%   -lvol   : volume number
%   -sarr   : is the array of values for the s-coordinate
%   -tarr   : is the array of values for the theta-coordinate
%   -zarr   : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -modB   : array of |B| with size ns*nt*nz where ns=size(sarr), nt=size(tarr), nz=size(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2016)

ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

modB    = zeros(ns,nt,nz); % allocate data for mod(B)


% Calculate the magnetic field contravariant components

bvec = get_spec_magfield(data,lvol,sarr,tarr,zarr); 
gmat = get_spec_metric(data,lvol,sarr,tarr,zarr);


% Calculate modB = sqrt(Bcontrav*gmat*Bcontrav)

for is=1:ns
 for it=1:nt
  for iz=1:nz
   a0 = gmat{1,1}(is,it,iz);
   b0 = gmat{1,2}(is,it,iz);
   c0 = gmat{1,3}(is,it,iz);
   d0 = gmat{2,1}(is,it,iz);
   e0 = gmat{2,2}(is,it,iz);
   f0 = gmat{2,3}(is,it,iz);
   g0 = gmat{3,1}(is,it,iz);
   h0 = gmat{3,2}(is,it,iz);
   i0 = gmat{3,3}(is,it,iz);
 
   bs = bvec{1}(is,it,iz);
   bt = bvec{2}(is,it,iz);
   bz = bvec{3}(is,it,iz);

   modB(is,it,iz) =  sqrt(( bs*(a0*bs+b0*bt+c0*bz)  ...
                           +bt*(d0*bs+e0*bt+f0*bz)  ...
                           +bz*(g0*bs+h0*bt+i0*bz)) ...
                          );
  end
 end
end


