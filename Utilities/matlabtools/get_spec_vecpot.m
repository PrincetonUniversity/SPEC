function Acov = get_spec_vecpot(data,lvol,sarr,tarr,zarr)
 
%
% GET_SPEC_VECPOT( DATA, LVOL, SARR, TARR, ZARR )
% ===============================================
%
% Computes covariant components of A in volume lvol
%
% INPUT
% -----
%   -data     : must be produced by calling read_spec(filename)
%   -lvol     : is the volume number 
%   -sarr     : is the array of values for the s-coordinate
%   -tarr     : is the array of values for the theta-coordinate
%   -zarr     : is the array of values for the zeta-coordinate
%
% OUTPUT
% ------
%   -Acov     : cell structure with 2 arrays: A_theta, A_zeta each with size length(sarr)*length(tarr)*length(zarr)
%
% written by J.Loizu (2018)


Ate     = data.vector_potential.Ate{lvol};
Aze     = data.vector_potential.Aze{lvol};
Ato     = data.vector_potential.Ato{lvol};
Azo     = data.vector_potential.Azo{lvol};

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

T       = cell(Lrad+1,1);  % allocate data for Chebyshev polynomials 
fac     = cell(mn,1);      % allocate data for regularization factors 

T{1}    = ones(ns,1);

T{2}    = sarr;



% Construct Chebyshev polynomials 

for l=3:Lrad+1
  T{l} = 2*sarr.*T{l-1} - T{l-2};
end


% Construct regularization factors

for j=1:mn
  if(lvol>1 || im(j)==0) 
   fac{j}  = ones(ns,1);
  elseif(im(j)==2)
   fac{j}  = sbar;
  else
   fac{j}  = sbar.^(im(j)/2);
  end
end


% Construct vector potential covariant components

for l=1:Lrad+1
  for j=1:mn
    for it=1:nt
      for iz=1:nz
       cosa = cos(im(j)*tarr(it)-in(j)*zarr(iz));
       sina = sin(im(j)*tarr(it)-in(j)*zarr(iz));
       At(:,it,iz) = At(:,it,iz) + fac{j}.*T{l}.*( Ate(l,j)*cosa + Ato(l,j)*sina );
       Az(:,it,iz) = Az(:,it,iz) + fac{j}.*T{l}.*( Aze(l,j)*cosa + Azo(l,j)*sina );
      end
    end
  end
end


Acov{1} = At;
Acov{2} = Az; 
