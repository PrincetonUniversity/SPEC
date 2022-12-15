function Acov = get_spec_vecpot_slab(data,lvol,sarr,tarr,zarr)
 
 
% Computes covariant components of A in volume lvol in slab geometry
%   -fdata  must be produced by calling read_spec_field(filename)
%   -lvol   is the volume number (from 1 to Nvol)
%   -sarr   is the array of values for the s-coordinate
%   -tarr   is the array of values for the theta-coordinate
%   -zarr   is the array of values for the zeta-coordinate
%   written by J.Loizu (2018)

fdata   = data.vector_potential;

Ate     = fdata.Ate{lvol};
Aze     = fdata.Aze{lvol};
Ato     = fdata.Ato{lvol};
Azo     = fdata.Azo{lvol};

Lrad    = data.input.physics.Lrad(lvol);

sarr    = transpose(sarr);
ns      = length(sarr);
nt      = length(tarr);
nz      = length(zarr);
sbar    = (sarr+1)/2;

mn      = data.output.mn;
im      = double(data.output.im);
in      = double(data.output.in);

At      = zeros(ns,nt,nz); % allocate data for vector potential along theta
Az      = zeros(ns,nt,nz); % allocate data for vector potential along zeta
dAt     = zeros(ns,nt,nz); % allocate data for vector potential radial derivative
dAz     = zeros(ns,nt,nz); % allocate data for vector potential radial derivative

T       = cell(Lrad+1,2);  % allocate data for Chebyshev polynomials and their derivatives
fac     = cell(mn,2);      % allocate data for regularization factors and their derivatives

T{1}{1} = ones(ns,1);
T{1}{2} = zeros(ns,1);

T{2}{1} = sarr;
T{2}{2} = ones(ns,1);


% Construct Chebyshev polynomials and their derivatives

for l=3:Lrad+1
  T{l}{1} = 2*sarr.*T{l-1}{1} - T{l-2}{1};
  T{l}{2} = 2*T{l-1}{1} + 2*sarr.*T{l-1}{2} - T{l-2}{2};
end


% Construct regularization factors and their derivatives

for j=1:mn
  fac{j}{1}  = ones(ns,1);
  fac{j}{2}  = zeros(ns,1);
end


% Construct magnetic field contravariant components

for l=1:Lrad+1
  for j=1:mn
    for it=1:nt
      for iz=1:nz
       cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
       sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
       At(:,it,iz) = At(:,it,iz) + fac{j}{1}.*T{l}{1}.*( Ate(l,j)*cosa + Ato(l,j)*sina );
       Az(:,it,iz) = Az(:,it,iz) + fac{j}{1}.*T{l}{1}.*( Aze(l,j)*cosa + Azo(l,j)*sina );
       dAt(:,it,iz) = dAt(:,it,iz) + fac{j}{1}.*T{l}{2}.*( Ate(l,j)*cosa + Ato(l,j)*sina );
       dAz(:,it,iz) = dAz(:,it,iz) + fac{j}{1}.*T{l}{2}.*( Aze(l,j)*cosa + Azo(l,j)*sina );
      end
    end
  end
end


Acov{1} = At; 
Acov{2} = Az; 
Acov{3} = dAt; 
Acov{4} = dAz; 


