function modBp = get_spec_modBpol(fdata,lvol,sarr,tarr,zarr)
 
 
% Calculates mod(Bpol) in volume lvol 
%
% INPUT
%   -fdata  : must be produced by calling read_spec_field(filename)
%   -lvol   : volume number
%   -sarr   : is the array of values for the s-coordinate
%   -tarr   : is the array of values for the theta-coordinate
%   -zarr   : is the array of values for the zeta-coordinate
%
% OUTPUT
%   -modBp  : array of |Bpol| with size ns*nt*nz where ns=size(sarr), nt=size(tarr), nz=size(zarr)
%
% Note: Stellarator symmetry is assumed
%
% written by J.Loizu (2019)


ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);

modBp   = zeros(ns,nt,nz); % allocate data for mod(Bpol)


% Calculate the magnetic field contravariant components

if(fdata.Igeometry == 1)
bvec = get_spec_magfield_slab(fdata,lvol,sarr,tarr,zarr); 
else
bvec = get_spec_magfield(fdata,lvol,sarr,tarr,zarr); 
end

% Calculate the metric matrix 

if(fdata.Igeometry == 1)
gmat = get_spec_metric_slab(fdata,lvol,sarr,tarr,zarr); 
else
gmat = get_spec_metric(fdata,lvol,sarr,tarr,zarr); 
end

% Calculate modB = sqrt(Bcontrav*gmat*Bcontrav)

for is=1:ns
 for it=1:nt
  for iz=1:nz
   a0 = gmat{1}{1}(is,it,iz);
   b0 = gmat{1}{2}(is,it,iz);
   d0 = gmat{2}{1}(is,it,iz);
   e0 = gmat{2}{2}(is,it,iz);
 
   bs = bvec{1}(is,it,iz);
   bt = bvec{2}(is,it,iz);

   modBp(is,it,iz) =  sqrt(( bs*(a0*bs+b0*bt)  ...
                           +bt*(d0*bs+e0*bt)  ...
                        ));
  end
 end
end


