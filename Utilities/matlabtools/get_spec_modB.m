function modB = get_spec_modB(fdata,lvol,sarr,tarr,zarr)
 
 
% Calculates mod(B) in volume lvol 
%
% INPUT
%   -fdata  : must be produced by calling read_spec_field(filename)
%   -lvol   : volume number
%   -sarr   : is the array of values for the s-coordinate
%   -tarr   : is the array of values for the theta-coordinate
%   -zarr   : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -modB   : array of |B| with size ns*nt*nz where ns=size(sarr), nt=size(tarr), nz=size(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2016)


ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

J       = zeros(ns,nt,nz); % allocate data for the Jacobian
modB    = zeros(ns,nt,nz); % allocate data for mod(B)


% Calculate the magnetic field contravariant components

bvec = get_spec_magfield(fdata,lvol,sarr,tarr,zarr); 


% Calculate the metric matrix 

gmat = get_spec_metric(fdata,lvol,sarr,tarr,zarr); 


% Calculate modB = sqrt(Bcontrav*gmat*Bcontrav)

for is=1:ns
 for it=1:nt
  for iz=1:nz
   a = gmat{1}{1}(is,it,iz);
   b = gmat{1}{2}(is,it,iz);
   c = gmat{1}{3}(is,it,iz);
   d = gmat{2}{1}(is,it,iz);
   e = gmat{2}{2}(is,it,iz);
   f = gmat{2}{3}(is,it,iz);
   g = gmat{3}{1}(is,it,iz);
   h = gmat{3}{2}(is,it,iz);
   i = gmat{3}{3}(is,it,iz);
 
   bs = bvec{1}(is,it,iz);
   bt = bvec{2}(is,it,iz);
   bz = bvec{3}(is,it,iz);

   modB(is,it,iz) =  sqrt(( bs*(a*bs+b*bt+c*bz)  ...
                           +bt*(d*bs+e*bt+f*bz)  ...
			   +bz*(g*bs+h*bt+i*bz)) ...
                        );
  end
 end
end


